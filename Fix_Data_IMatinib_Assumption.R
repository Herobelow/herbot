# =============================================================================
#  DATA CORRECTION SCRIPT - IMATINIB 1L ASSUMPTION + FOLLOW-UP FLOW AUDIT
#  KNH CML TKI STUDY
#
#  WHAT IT DOES
#    1. Finds EVERY cell containing
#         "Not applicable (patient not treated with imatinib)"
#       (and wording variants) and corrects it under the study assumption
#       that EVERY patient started IMATINIB as 1L and then progressed to
#       later-line TKIs:
#         - exposure-type columns   -> Yes / Imatinib (assumed 1L)
#         - duration-type columns   -> left BLANK (dates were never
#                                      captured - a number cannot be
#                                      invented) + logged
#         - response-type columns   -> "Not documented"
#         - anything else           -> "Not documented (REVIEW)"
#    2. Audits the whole dataset for follow-up flow consistency
#       (survival status <-> OS time <-> cause of death; PFS event <->
#       PFS time; Grade 3+ ADR vs Any ADR; relapse vs TKI modification;
#       phase vs TKI; age range; duplicate IDs; T315I vs TKI).
#       Contradictions are FLAGGED (your call), never silently changed.
#    3. Writes two files next to your data file (original untouched):
#         <yourfile>_CORRECTED.csv   - the corrected dataset
#         DATA_CORRECTION_AUDIT.csv  - every change + every flagged issue
#                                      (your paper trail for the panel)
#
#  HOW TO RUN: RStudio -> File > Open File -> select this script ->
#              Run All -> pick your CSV (e.g. "am at this point 2.csv")
#              in the window that opens.
# =============================================================================

# ---------------------------------------------------------------------------
# 0. PACKAGES
# ---------------------------------------------------------------------------
for (pkg in c("readr", "dplyr", "stringr", "tibble")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(">>> Installing missing package: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# ---------------------------------------------------------------------------
# 1. FILE WINDOW
# ---------------------------------------------------------------------------
csv_path <- tryCatch(file.choose("Select your CML data file (.csv)"),
                     error = function(e) NA_character_)
if (is.na(csv_path) || !nzchar(csv_path) || !file.exists(csv_path)) {
  typed <- tryCatch(input("Full path to .csv file: "), error = function(e) "")
  typed <- gsub('"', "", trimws(as.character(typed)))
  if (nzchar(typed)) csv_path <- typed
}
if (!file.exists(csv_path)) {
  csv_path <- "am at this point 2.csv"
  if (!file.exists(csv_path)) stop("No data file found - select it in the window.", call. = FALSE)
}
if (tolower(tools::file_ext(csv_path)) %in% c("xlsx", "xls")) {
  stop("Selected file is an Excel workbook. Save it as CSV (comma delimited) first.", call. = FALSE)
}

message("\n>>> Loading: ", normalizePath(csv_path))
df <- tryCatch(
  readr::read_csv(csv_path, show_col_types = FALSE,
                  locale = readr::locale(encoding = "UTF-8")),
  error = function(e) readr::read_csv(csv_path, show_col_types = FALSE,
                                      locale = readr::locale(encoding = "Latin1"))
)
names(df) <- trimws(names(df))
message(">>> Loaded ", nrow(df), " rows x ", ncol(df), " columns.")
if (ncol(df) < 5 || nrow(df) < 5) {
  stop("Loaded file looks wrong (", nrow(df), " x ", ncol(df), ").", call. = FALSE)
}

corrected <- df

# Output location: next to the data file
out_dir <- dirname(normalizePath(csv_path))
base    <- tools::file_path_sans_ext(tools::file_path_sans_ext(basename(csv_path)))
fixed_f <- file.path(out_dir, paste0(base, "_CORRECTED.csv"))
audit_f <- file.path(out_dir, "DATA_CORRECTION_AUDIT.csv")

# ---------------------------------------------------------------------------
# 2. COLUMN MAP (automatic)
# ---------------------------------------------------------------------------
find_col <- function(df, exact = NULL, patterns = NULL) {
  nms <- names(df)
  if (!is.null(exact)) {
    hit <- intersect(exact, nms)
    if (length(hit) > 0) return(hit[1])
  }
  if (!is.null(patterns)) {
    for (pat in patterns) {
      hit <- nms[grepl(pat, nms, ignore.case = TRUE, perl = TRUE)]
      if (length(hit) > 0) return(hit[1])
    }
  }
  NULL
}

colmap <- list(
  study_id  = find_col(df, c("Study_ID", "Study ID"), patterns = c("Study[_ ]?ID")),
  tki       = find_col(df, c("TKI Used", "TKI"), patterns = c("TKI[_ ]?Used", "^TKI$")),
  survival  = find_col(df, c("Survival_Status", "Survival Status"), patterns = c("Survival[_ ]?Status")),
  os_time   = find_col(df, c("OS_Time_Months"), patterns = c("OS[_ ]?Time")),
  pfs_event = find_col(df, c("PFS_Event"), patterns = c("PFS[_ ]?Event")),
  pfs_time  = find_col(df, c("PFS_Time_Months"), patterns = c("PFS[_ ]?Time")),
  age       = find_col(df, c("Age"), patterns = c("^Age$")),
  phase     = find_col(df, c("Phase of CML", "CML Phase"), patterns = c("Phase")),
  t315i     = find_col(df, c("T315I_Mutation", "T315I"), patterns = c("T315I")),
  adr_any   = find_col(df, c("Any ADRs", "Any ADR"), patterns = c("Any[_ ]?ADR")),
  adr_sev   = find_col(df, c("Grade_3_ADR", "Grade 3 ADR"), patterns = c("Grade[_ ]?[33]")),
  relapse   = find_col(df, c("Disease_Relapse", "Disease Relapse"), patterns = c("Relapse")),
  mod       = find_col(df, c("TKI_Treatment_Modification"), patterns = c("Modification")),
  cause     = find_col(df, c("Cause of Death", "Cause_Of_Death"), patterns = c("cause[_ ]?of[_ ]?death"))
)

study_id_at <- function(row) {
  if (is.null(colmap$study_id)) return(NA_character_)
  as.character(corrected[[colmap$study_id]][row])
}

# ---------------------------------------------------------------------------
# 3. AUDIT LOG
# ---------------------------------------------------------------------------
audit <- tibble(Row = integer(), Study_ID = character(), Column = character(),
                Old_Value = character(), New_Value = character(),
                Action = character(), Reason = character())

add_audit <- function(row, col, oldv, newv, action, reason) {
  audit <<- bind_rows(audit, tibble(
    Row      = as.integer(row),
    Study_ID = as.character(study_id_at(row)),
    Column   = col,
    Old_Value = if (is.null(oldv)) "" else paste(as.character(oldv), collapse = " | "),
    New_Value = if (is.null(newv)) "" else paste(as.character(newv), collapse = " | "),
    Action   = action,
    Reason   = reason))
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# 4. STEP A - FIND EVERY "NOT TREATED WITH IMATINIB" CELL
# ---------------------------------------------------------------------------
na_pat <- "not\\s+(treated|on|exposed to)\\s+with\\s*imatinib|never\\s+(treated|on|exposed to)\\s+imatinib"

message("\n=============================================================")
message(" STEP 1: 'not treated with imatinib' cells")
message("=============================================================")

im_hits <- list()   # column -> rows
for (col in names(df)) {
  v <- as.character(df[[col]])
  m <- !is.na(v) & grepl(na_pat, tolower(v), perl = TRUE)
  if (any(m)) im_hits[[col]] <- which(m)
}

if (length(im_hits) == 0) {
  message(">>> No 'not treated with imatinib' cells found in this file.")
} else {
  message(">>> Found in ", length(im_hits), " column(s):\n")
  for (col in names(im_hits)) {
    message("    - ", col, "  (", length(im_hits[[col]]), " row(s))")
  }
}

# ---------------------------------------------------------------------------
# 5. STEP B - CORRECT EACH HIT (column-type rules)
# ---------------------------------------------------------------------------
# Returns a class label for the column.
classify_col <- function(col) {
  cl <- tolower(col)
  dur_kw  <- c("duration", "time", "months", "free", "interval", "follow")
  resp_kw <- c("response", "mmr", "cmr", "molecular", "haematolog", "haematologic", "depth")
  exp_kw  <- c("treated", "therapy", "treatment", "line", "use")
  if (grepl("imatinib", cl)) {
    if (any(vapply(dur_kw, function(k) grepl(k, cl), logical(1))))  return("duration")
    if (any(vapply(resp_kw, function(k) grepl(k, cl), logical(1)))) return("response")
    if (any(vapply(exp_kw, function(k) grepl(k, cl), logical(1))))  return("exposure")
    return("other")
  }
  # Column mentions the label but not the word imatinib (e.g. "1L therapy")
  if (any(vapply(exp_kw, function(k) grepl(k, cl), logical(1)))) return("exposure")
  "other"
}

# Is the column essentially a Yes/No column?
looks_binary <- function(col) {
  vals <- tolower(trimws(as.character(corrected[[col]])))
  vals <- vals[!is.na(vals) & vals != ""]
  vals <- vals[!grepl(na_pat, vals, perl = TRUE)]
  if (length(vals) == 0) return(TRUE)
  mean(vals %in% c("yes", "no", "y", "n", "1", "0", "true", "false")) >= 0.8
}

n_fixed <- n_missing <- n_review <- 0
for (col in names(im_hits)) {
  rows  <- im_hits[[col]]
  ctype <- classify_col(col)
  is_num <- is.numeric(corrected[[col]])

  for (r in rows) {
    oldv <- as.character(corrected[[col]][r])
    if (ctype == "exposure") {
      newv <- if (looks_binary(col)) "Yes" else "Imatinib (assumed 1L)"
      corrected[[col]][r] <- newv
      n_fixed <<- n_fixed + 1
      add_audit(r, col, oldv, newv, "CORRECTED",
                "Study assumption: every patient started imatinib as 1L")
    } else if (ctype == "duration") {
      if (is_num) {
        corrected[[col]][r] <- NA_real_
        n_missing <<- n_missing + 1
        add_audit(r, col, oldv, NA, "SET_MISSING",
                  "Imatinib exposure assumed (1L) but start/relapse dates were never captured - a duration cannot be computed; left blank")
      } else {
        newv <- "Not documented (imatinib start/relapse dates not captured)"
        corrected[[col]][r] <- newv
        n_missing <<- n_missing + 1
        add_audit(r, col, oldv, newv, "SET_MISSING",
                  "Imatinib exposure assumed (1L) but start/relapse dates were never captured")
      }
    } else if (ctype == "response") {
      newv <- "Not documented"
      corrected[[col]][r] <- newv
      n_missing <<- n_missing + 1
      add_audit(r, col, oldv, newv, "SET_MISSING",
                "Imatinib exposure assumed (1L); response outcome not captured in dataset")
    } else {
      newv <- "Not documented (REVIEW)"
      corrected[[col]][r] <- newv
      n_review <<- n_review + 1
      add_audit(r, col, oldv, newv, "REVIEW",
                "Imatinib-related 'not treated' label found - replaced with 'Not documented (REVIEW)'; check this column")
    }
  }
}
message("\n>>> Corrections: ", n_fixed, " exposure cells set to Yes/Imatinib; ",
        n_missing, " cells set to missing/Not documented; ", n_review, " flagged for review.")

# ---------------------------------------------------------------------------
# 6. STEP C - FOLLOW-UP FLOW CONSISTENCY AUDIT (flags only)
# ---------------------------------------------------------------------------
message("\n=============================================================")
message(" STEP 2: follow-up flow consistency audit")
message("=============================================================")

flag_count <- 0
blank_like <- function(x) {
  v <- trimws(as.character(x))
  is.na(v) | v %in% c("", "na", "n/a", "not documented", "not applicable", "none", "-")
}
val_low <- function(x) tolower(trimws(as.character(x)))

for (r in seq_len(nrow(corrected))) {
  sid <- study_id_at(r)

  if (!is.null(colmap$survival)) {
    surv <- val_low(corrected[[colmap$survival]][r])
    os_t <- suppressWarnings(as.numeric(corrected[[colmap$os_time]][r]))
    caus <- if (!is.null(colmap$cause)) corrected[[colmap$cause]][r] else NA
    is_dead  <- surv %in% c("deceased", "dead", "died")
    is_alive <- surv %in% c("alive", "living", "survived", "on treatment")
    if (is_dead) {
      if (is.null(colmap$os_time) || is.na(os_t) || blank_like(corrected[[colmap$os_time]][r])) {
        add_audit(r, "OS_Time_Months", surv, NA, "FLAGGED",
                  "Recorded deceased but OS time missing"); flag_count <<- flag_count + 1
      }
      if (!is.null(colmap$cause) && (blank_like(caus))) {
        add_audit(r, "Cause of Death", surv, NA, "FLAGGED",
                  "Recorded deceased but cause of death missing"); flag_count <<- flag_count + 1
      }
    }
    if (is_alive && !is.null(colmap$cause) && !blank_like(caus)) {
      add_audit(r, "Cause of Death", surv, caus, "FLAGGED",
                "Recorded alive but cause of death is filled"); flag_count <<- flag_count + 1
    }
    if (!is_dead && !is_alive && !blank_like(surv)) {
      add_audit(r, colmap$survival, surv, NA, "FLAGGED",
                "Unrecognised survival status value - check spelling"); flag_count <<- flag_count + 1
    }
  }

  if (!is.null(colmap$os_time)) {
    os_t <- suppressWarnings(as.numeric(corrected[[colmap$os_time]][r]))
    if (!is.na(os_t) && os_t < 0) {
      add_audit(r, colmap$os_time, os_t, NA, "FLAGGED",
                "Negative OS time"); flag_count <<- flag_count + 1
    }
  }

  if (!is.null(colmap$pfs_event) && !is.null(colmap$pfs_time)) {
    pe <- val_low(corrected[[colmap$pfs_event]][r])
    pt <- suppressWarnings(as.numeric(corrected[[colmap$pfs_time]][r]))
    if (pe %in% c("1", "true", "yes") && (is.na(pt) || blank_like(corrected[[colmap$pfs_time]][r]))) {
      add_audit(r, "PFS_Time_Months", pe, NA, "FLAGGED",
                "PFS event = 1 but PFS time missing"); flag_count <<- flag_count + 1
    }
    if (!is.na(pt) && pt < 0) {
      add_audit(r, colmap$pfs_time, pt, NA, "FLAGGED",
                "Negative PFS time"); flag_count <<- flag_count + 1
    }
    os_t <- if (!is.null(colmap$os_time)) suppressWarnings(as.numeric(corrected[[colmap$os_time]][r])) else NA
    if (!is.na(pt) && !is.na(os_t) && pt > os_t && os_t >= 0) {
      add_audit(r, "PFS_Time_Months", pt, os_t, "FLAGGED",
                "PFS time later than OS (death) time"); flag_count <<- flag_count + 1
    }
  }

  if (!is.null(colmap$adr_any) && !is.null(colmap$adr_sev)) {
    sev <- val_low(corrected[[colmap$adr_sev]][r])
    any <- val_low(corrected[[colmap$adr_any]][r])
    if (sev %in% c("yes", "y", "1", "true") && !(any %in% c("yes", "y", "1", "true"))) {
      add_audit(r, colmap$adr_any, any, sev, "FLAGGED",
                "Grade 3+ ADR = Yes but 'Any ADRs' is not Yes (should it be Yes?)")
      flag_count <<- flag_count + 1
    }
  }

  if (!is.null(colmap$relapse) && !is.null(colmap$mod)) {
    rel <- val_low(corrected[[colmap$relapse]][r])
    mod_v <- as.character(corrected[[colmap$mod]][r])
    if (rel %in% c("yes", "y", "1", "true") && blank_like(mod_v)) {
      add_audit(r, colmap$mod, "relapse=Yes", NA, "FLAGGED",
                "Relapse = Yes but no TKI modification recorded (was a switch missed?)")
      flag_count <<- flag_count + 1
    }
  }

  if (!is.null(colmap$phase) && !is.null(colmap$tki)) {
    ph <- val_low(corrected[[colmap$phase]][r])
    tk <- val_low(corrected[[colmap$tki]][r])
    if (grepl("blast", ph, fixed = TRUE) && tk %in% c("imatinib", "")) {
      add_audit(r, colmap$phase, paste(ph, tk), NA, "FLAGGED",
                "Blast phase recorded with imatinib (monotherapy) - verify")
      flag_count <<- flag_count + 1
    }
  }

  if (!is.null(colmap$age)) {
    ag <- suppressWarnings(as.numeric(corrected[[colmap$age]][r]))
    if (is.na(ag)) {
      add_audit(r, colmap$age, corrected[[colmap$age]][r], NA, "FLAGGED",
                "Age missing or not numeric"); flag_count <<- flag_count + 1
    } else if (ag < 18 || ag > 95) {
      add_audit(r, colmap$age, ag, NA, "FLAGGED",
                "Age outside plausible range 18-95"); flag_count <<- flag_count + 1
    }
  }

  if (!is.null(colmap$t315i) && !is.null(colmap$tki)) {
    t3 <- val_low(corrected[[colmap$t315i]][r])
    tk <- val_low(corrected[[colmap$tki]][r])
    if (t3 %in% c("yes", "y", "1", "positive") && tk != "ponatinib") {
      add_audit(r, colmap$tki, paste("T315I=Yes", tk), NA, "FLAGGED",
                "T315I+ but recorded TKI is not ponatinib (see T315I report)")
      flag_count <<- flag_count + 1
    }
  }

  if (!is.null(colmap$tki) && blank_like(as.character(corrected[[colmap$tki]][r]))) {
    add_audit(r, colmap$tki, NA, NA, "FLAGGED",
              "No TKI recorded"); flag_count <<- flag_count + 1
  }
}

if (!is.null(colmap$study_id)) {
  sids <- as.character(corrected[[colmap$study_id]])
  dups <- sids[!is.na(sids) & duplicated(sids)]
  if (length(dups) > 0) {
    for (r in which(sids %in% dups)) {
      add_audit(r, colmap$study_id, sids[r], NA, "FLAGGED", "Duplicate Study_ID")
      flag_count <<- flag_count + 1
    }
  }
}

message(">>> Consistency flags raised: ", flag_count)

# ---------------------------------------------------------------------------
# 7. WRITE OUTPUTS
# ---------------------------------------------------------------------------
write.csv(corrected, fixed_f, row.names = FALSE)
write.csv(audit, audit_f, row.names = FALSE)

message("\n=============================================================")
message(" DONE - files written next to your data file")
message("=============================================================")
message("  Corrected dataset: ", fixed_f)
message("  Audit log:         ", audit_f)
message("\nSummary of audit log:")
if (nrow(audit) > 0) {
  print(audit %>% count(Action, Column, name = "n", sort = TRUE))
}
message("\n>>> Open DATA_CORRECTION_AUDIT.csv to review every change and every")
message(">>> flagged inconsistency, then send it (or paste its summary) for review.")
