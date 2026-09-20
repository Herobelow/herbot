# =============================================================================
#  T315I MUTATION REPORT  -  KNH CML TKI STUDY  (panel comment i)
#
#  Standalone script. It answers the panel's two questions:
#    (1) Which mutation was most frequently detected?
#    (2) What therapy was used in the T315I-positive patients - in
#        particular ponatinib, the TKI that is active against T315I?
#
#  HOW TO RUN IN RSTUDIO
#    1. File > Open File...  -> select this script  (or paste it in).
#    2. Press "Run All" (or Ctrl+Shift+S).
#    3. A file-selection WINDOW will pop up - choose your data file
#       (e.g. "am at this point 2.csv").
#       If you close the window, you will be asked to TYPE the full path.
#    4. Every table is printed in the Console AND saved as a CSV in an
#       "output" folder next to your data file.
#
#  NOTES
#    - Reads .csv files only. If you select an Excel (.xlsx) file the
#      script stops and tells you how to save it as CSV.
#    - Column names are matched automatically - no editing needed.
#    - The main analysis script (CML_TKI_KNH_Main_Analysis.R) is unchanged;
#      this file only produces the T315I mutation report.
# =============================================================================


# ---------------------------------------------------------------------------
# 0. PACKAGES
# ---------------------------------------------------------------------------
required_pkgs <- c("readr", "dplyr", "stringr", "tibble")
for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(">>> Installing missing package: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# survival is OPTIONAL - only needed for the bonus T315I+ vs T315I- survival
have_survival <- requireNamespace("survival", quietly = TRUE)
if (!have_survival) {
  install.packages("survival", repos = "https://cloud.r-project.org")
  have_survival <- requireNamespace("survival", quietly = TRUE)
}
if (have_survival) suppressPackageStartupMessages(library(survival))


# ---------------------------------------------------------------------------
# 1. FILE UPLOAD WINDOW
# ---------------------------------------------------------------------------
cat("\n=============================================================\n")
cat(" T315I MUTATION REPORT - file selection\n")
cat("=============================================================\n")

csv_path <- tryCatch(file.choose("Select your CML data file (.csv)"),
                     error = function(e) NA_character_)

if (is.na(csv_path) || !nzchar(csv_path) || !file.exists(csv_path)) {
  cat("No file was selected in the window.\n")
  cat("Type the FULL path to your .csv file (folder + filename) and press Enter:\n")
  typed <- tryCatch(input("Full path to .csv file: "), error = function(e) "")
  typed <- gsub('"', "", trimws(as.character(typed)))
  if (nzchar(typed)) csv_path <- typed
}

if (!file.exists(csv_path)) {
  cat("Falling back to the default filename...\n")
  csv_path <- "am at this point 2.csv"
  if (!file.exists(csv_path)) {
    stop("No data file found. Run the script again and select your .csv file in the window.", call. = FALSE)
  }
}

# Guard: CSV only - Excel workbooks are refused with clear instructions.
if (tolower(tools::file_ext(csv_path)) %in% c("xlsx", "xls")) {
  stop("The selected file is an Excel workbook (", basename(csv_path), ").\n",
       "This script reads CSV files. In Excel use File > Save As > CSV (comma\n",
       "delimited) (.csv), then select that .csv file in the window.", call. = FALSE)
}

message("\n>>> Loading data from: ", normalizePath(csv_path))
df <- tryCatch(
  readr::read_csv(csv_path, show_col_types = FALSE,
                  locale = readr::locale(encoding = "UTF-8")),
  error = function(e) {
    message("UTF-8 read failed (", conditionMessage(e), "); retrying with Latin-1.")
    readr::read_csv(csv_path, show_col_types = FALSE,
                    locale = readr::locale(encoding = "Latin1"))
  }
)
names(df) <- trimws(names(df))
message(">>> Loaded ", nrow(df), " rows x ", ncol(df), " columns.")

if (ncol(df) < 5 || nrow(df) < 5) {
  stop("Loaded file has ", nrow(df), " row(s) x ", ncol(df),
       " column(s) - it does not look like the CML dataset.\n",
       "Check that you selected the correct .csv file with a header row.", call. = FALSE)
}

# Output folder (next to the data file)
out_dir <- file.path(dirname(normalizePath(csv_path)), "output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
save_csv <- function(x, fname) {
  f <- file.path(out_dir, fname)
  write.csv(x, f, row.names = FALSE)
  message("    Saved: ", f)
  invisible(x)
}


# ---------------------------------------------------------------------------
# 2. COLUMN MAPPING (automatic)
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
  study_id  = find_col(df, c("Study_ID", "Study ID"),
                       patterns = c("Study[_ ]?ID")),
  t315i     = find_col(df, c("T315I_Mutation", "T315I Mutation", "T315I"),
                       patterns = c("T315I")),
  tki       = find_col(df, c("TKI Used", "TKI", "Current TKI",
                             "TKI (Current)", "TKI (current)"),
                       patterns = c("TKI[_ ]?Used", "Current[_ ]?TKI", "^TKI$")),
  phase     = find_col(df, c("Phase of CML", "CML Phase", "Phase"),
                       patterns = c("Phase")),
  survival  = find_col(df, c("Survival_Status", "Survival Status"),
                       patterns = c("Survival[_ ]?Status")),
  os_time   = find_col(df, c("OS_Time_Months", "OS Time (Months)", "OS_Time"),
                       patterns = c("OS[_ ]?Time")),
  cause     = find_col(df, c("Cause of Death", "Cause_Of_Death", "CauseOfDeath"),
                       patterns = c("cause[_ ]?of[_ ]?death")),
  age       = find_col(df, c("Age"), patterns = c("^Age$")),
  sex       = find_col(df, c("Gender", "Sex"), patterns = c("^(Gender|Sex)$"))
)

cat("\n>>> Columns matched:\n")
for (nm in names(colmap)) {
  v <- colmap[[nm]]
  cat("    ", nm, ": ", if (is.null(v)) "NOT FOUND" else v, "\n", sep = "")
}

stopifnot("No T315I column found in the file (looked for names containing 'T315I')." =
            !is.null(colmap$t315i))
stopifnot("No TKI column found in the file (looked for 'TKI Used' / 'Current TKI')." =
            !is.null(colmap$tki))


# ---------------------------------------------------------------------------
# 3. PREPARE VARIABLES
# ---------------------------------------------------------------------------
# --- TKI used + TKI generation (1G/2G/3G) ------------------------------------
tki <- tolower(str_trim(as.character(df[[colmap$tki]])))
tki[tki %in% c("", "na", "n/a", "not documented")] <- NA_character_
df$tki <- tki
df$gen <- factor(
  case_when(
    tki == "imatinib" ~ "1G",
    tki %in% c("dasatinib", "nilotinib", "bosutinib") ~ "2G",
    tki == "ponatinib" ~ "3G",
    TRUE ~ NA_character_
  ),
  levels = c("1G", "2G", "3G")
)
gen_desc <- c("1G" = "imatinib",
              "2G" = "dasatinib / nilotinib / bosutinib",
              "3G" = "ponatinib")

# --- T315I status (4-level text) + positive flag -----------------------------
t3 <- tolower(str_trim(as.character(df[[colmap$t315i]])))
t3[t3 %in% c("", "nd", "not documented")] <- "Not documented"
t3[t3 %in% c("na", "n/a", "not applicable")] <- "Not applicable"
t3[t3 %in% c("yes", "y", "1", "positive")] <- "Yes"
t3[t3 %in% c("no", "n", "0", "negative")] <- "No"
df$T315I_cat <- factor(t3, levels = c("Yes", "No", "Not applicable", "Not documented"))
df$flag_T315I <- df$T315I_cat == "Yes"

# --- survival status + OS time ----------------------------------------------
if (!is.null(colmap$survival)) {
  surv_raw <- tolower(str_trim(as.character(df[[colmap$survival]])))
  df$died <- surv_raw %in% c("deceased", "dead", "died")
} else {
  df$died <- NA
}
if (!is.null(colmap$os_time)) {
  df$OS_time <- suppressWarnings(as.numeric(as.character(df[[colmap$os_time]])))
} else {
  df$OS_time <- NA_real_
}

# --- helper columns for the patient list -------------------------------------
if (!is.null(colmap$study_id)) df$Study_ID <- df[[colmap$study_id]] else df$Study_ID <- NA
if (!is.null(colmap$age))      df$Age      <- suppressWarnings(as.numeric(as.character(df[[colmap$age]]))) else df$Age <- NA
if (!is.null(colmap$sex))      df$Sex      <- str_trim(as.character(df[[colmap$sex]])) else df$Sex <- NA
if (!is.null(colmap$phase))    df$Phase    <- str_trim(as.character(df[[colmap$phase]])) else df$Phase <- NA
if (!is.null(colmap$cause))    df$Cause_of_Death <- str_trim(as.character(df[[colmap$cause]])) else df$Cause_of_Death <- NA

message(">>> Preparation complete.")


# ---------------------------------------------------------------------------
# 4. REPORT A - T315I TESTING AND DETECTION
# ---------------------------------------------------------------------------
cat("\n=============================================================\n")
cat(" REPORT A: T315I TESTING AND DETECTION\n")
cat("=============================================================\n")

t3_status <- df %>%
  mutate(T315I_cat = ifelse(is.na(T315I_cat), "Not documented", as.character(T315I_cat))) %>%
  count(T315I_cat, name = "n", sort = TRUE) %>%
  mutate(Pct_of_cohort = round(100 * n / nrow(df), 1))

n_tested <- sum(df$T315I_cat %in% c("Yes", "No"))
n_pos    <- sum(df$T315I_cat == "Yes")
t3_status <- bind_rows(
  t3_status,
  tibble(T315I_cat = "TESTED (Yes or No documented)", n = n_tested,
         Pct_of_cohort = round(100 * n_tested / nrow(df), 1))
)
save_csv(t3_status, "T315I_01_Testing_and_Detection.csv")
print(t3_status)

message("\n>>> ", n_tested, " of ", nrow(df), " patients (",
        round(100 * n_tested / nrow(df), 1), "%) had T315I testing documented.")
message(">>> T315I was DETECTED in ", n_pos, " patients (",
        round(100 * n_pos / nrow(df), 1), "% of the cohort",
        if (n_tested > 0) paste0("; ", round(100 * n_pos / n_tested, 1), "% of those tested") else "",
        ").")


# ---------------------------------------------------------------------------
# 5. REPORT B - WHICH MUTATION WAS MOST FREQUENTLY DETECTED?
# ---------------------------------------------------------------------------
cat("\n=============================================================\n")
cat(" REPORT B: MUTATION FREQUENCIES\n")
cat("=============================================================\n")

# All columns that look like mutation columns (the T315I column is handled
# separately above; anything else is examined for additional mutations).
mut_cols <- setdiff(grep("mutat", names(df), ignore.case = TRUE, value = TRUE),
                   colmap$t315i)

neg_vals <- c("na", "n/a", "none", "nd", "negative", "not documented",
              "not applicable", "no mutation", "no mutation detected")

mut_freq <- NULL
if (length(mut_cols) > 0) {
  message(">>> Additional mutation column(s) found: ", paste(mut_cols, collapse = ", "))
  rows <- list()
  for (col in mut_cols) {
    vals <- tolower(str_trim(as.character(df[[col]])))
    vals[is.na(vals) | vals %in% neg_vals] <- NA_character_
    tab <- table(vals, useNA = "no")
    if (length(tab) > 0) {
      rows[[col]] <- tibble(Mutation_column = col,
                            Mutation = names(tab),
                            n = as.integer(tab),
                            Pct_of_cohort = round(100 * as.integer(tab) / nrow(df), 1))
    }
  }
  if (length(rows) > 0) mut_freq <- bind_rows(rows)
} else {
  message(">>> No additional mutation columns found in the file.")
  message(">>> The dataset records ONLY the T315I mutation, so it is the only")
  message(">>> mutation for which frequency can be reported.")
}

if (!is.null(mut_freq) && nrow(mut_freq) > 0) {
  save_csv(mut_freq, "T315I_02_Mutation_Frequencies.csv")
  print(mut_freq)
  top <- mut_freq[which.max(mut_freq$n), ]
  message("\n>>> Most frequently detected mutation: ", top$Mutation,
          " (", top$n, " patient(s)).")
} else {
  save_csv(tibble(Mutation = "T315I (BCR-ABL1 p.251)", n = n_pos,
                  Pct_of_cohort = round(100 * n_pos / nrow(df), 1),
                  Pct_of_tested = if (n_tested > 0) round(100 * n_pos / n_tested, 1) else NA_real_),
           "T315I_02_Mutation_Frequencies.csv")
  message("\n>>> Most frequently (and only) detected mutation: T315I - ",
          n_pos, " of ", nrow(df), " patients (",
          round(100 * n_pos / nrow(df), 1), "% of the cohort",
          if (n_tested > 0) paste0(", ", round(100 * n_pos / n_tested, 1), "% of those tested") else "",
          ").")
}


# ---------------------------------------------------------------------------
# 6. REPORT C - THE T315I-POSITIVE PATIENTS (per-patient list)
# ---------------------------------------------------------------------------
cat("\n=============================================================\n")
cat(" REPORT C: T315I-POSITIVE PATIENTS - PER-PATIENT LIST\n")
cat("=============================================================\n")

t3_pos <- df %>%
  filter(flag_T315I) %>%
  select(Study_ID, Age, Sex, Phase, tki, gen, died, OS_time, Cause_of_Death) %>%
  mutate(Survival = ifelse(died, "Deceased",
                           ifelse(is.na(died), NA_character_, "Alive")))

if (nrow(t3_pos) > 0) {
  save_csv(t3_pos, "T315I_03_Positive_Patients.csv")
  print(t3_pos)
} else {
  message(">>> No T315I-positive patients found.")
}


# ---------------------------------------------------------------------------
# 7. REPORT D - THERAPY USED IN T315I+ PATIENTS (the key panel question)
# ---------------------------------------------------------------------------
cat("\n=============================================================\n")
cat(" REPORT D: THERAPY USED IN T315I+ PATIENTS\n")
cat("=============================================================\n")

if (nrow(t3_pos) > 0) {

  # 7.1 TKI distribution among T315I+ patients (within group = 100%) --------
  therapy_dist <- t3_pos %>%
    count(tki, gen, name = "n", sort = TRUE) %>%
    mutate(Generation = ifelse(is.na(gen), NA_character_,
                               paste0(gen, " (", gen_desc[gen], ")")),
           Pct_of_T315I_positive = round(100 * n / nrow(t3_pos), 1)) %>%
    select(tki, Generation, n, Pct_of_T315I_positive)

  # 7.2 Concordance: is the recorded TKI active against T315I? --------------
  therapy_conc <- t3_pos %>%
    mutate(therapy_match = ifelse(tki == "ponatinib",
                                  "Appropriate: ponatinib (active vs T315I)",
                                  "TKI without proven activity vs T315I")) %>%
    select(Study_ID, tki, gen, Phase, OS_time, Survival, therapy_match)

  save_csv(therapy_dist, "T315I_04a_Therapy_Distribution.csv")
  save_csv(therapy_conc, "T315I_04b_Therapy_Concordance.csv")

  cat("\nTKI distribution among T315I-positive patients (group total = 100%):\n\n")
  print(therapy_dist)
  cat("\nPer-patient therapy concordance:\n\n")
  print(therapy_conc)

  n_pono  <- sum(t3_pos$tki == "ponatinib", na.rm = TRUE)
  n_other <- sum(!t3_pos$tki %in% c("ponatinib", NA), na.rm = TRUE)
  n_unk   <- sum(is.na(t3_pos$tki), na.rm = TRUE)

  message("\n>>> ", n_pono, " of ", nrow(t3_pos), " T315I+ patients (",
          round(100 * n_pono / nrow(t3_pos), 1), "%) were recorded on PONATINIB -",
          " the TKI effective against T315I.")
  if (n_other > 0) {
    other_tkis <- paste(sort(unique(t3_pos$tki[!t3_pos$tki %in% c("ponatinib", NA)])),
                        collapse = ", ")
    message(">>> ", n_other, " patient(s) were recorded on ", other_tkis,
            " - TKIs that lack activity against T315I.")
    message(">>> Interpretation: because the dataset records only the CURRENT TKI,",
            " these records may reflect T315I detected before a planned switch,")
    message(">>> or a documentation limitation. This must be stated as a study limitation.")
  }
  if (n_unk > 0) message(">>> ", n_unk, " patient(s) had no TKI recorded.")

} else {
  message(">>> No T315I-positive patients - therapy report not applicable.")
}


# ---------------------------------------------------------------------------
# 8. REPORT E (BONUS) - T315I+ vs T315I- OVERALL SURVIVAL
# ---------------------------------------------------------------------------
cat("\n=============================================================\n")
cat(" REPORT E (BONUS): T315I+ vs T315I- OVERALL SURVIVAL\n")
cat("=============================================================\n")

if (!have_survival) {
  message(">>> 'survival' package not available - survival comparison skipped.")
} else if (is.null(colmap$os_time) || is.null(colmap$survival)) {
  message(">>> OS time or survival-status column not found - comparison skipped.")
} else {
  d <- df %>%
    filter(!is.na(OS_time), !is.na(died),
           T315I_cat %in% c("Yes", "No")) %>%
    mutate(grp = factor(T315I_cat, levels = c("Yes", "No")))
  if (nrow(d) < 10 || sum(d$died) < 2) {
    message(">>> Too few T315I-tested patients with complete OS data (n = ", nrow(d),
            ", events = ", sum(d$died, na.rm = TRUE),
            ") - comparison skipped.")
  } else {
    tryCatch({
      logr <- survival::survdiff(survival::Surv(OS_time, died) ~ grp, data = d)
      p_lr <- 1 - pchisq(logr$chisq, df = length(logr$nevent) - 1)

      # KM estimate at time t (last estimate if follow-up is shorter)
      km_estimate_at <- function(time, event, t) {
        f <- survival::survfit(survival::Surv(time, event))
        tt <- f$time; s <- f$surv
        if (length(tt) == 0) return(100)
        idx <- which(tt <= t)
        if (length(idx) == 0) return(100)
        100 * s[max(idx)]
      }

      report_t <- if (max(d$OS_time) >= 60) 60 else max(d$OS_time)
      surv_out <- tibble(
        Group = levels(d$grp),
        n = t(d$grp),
        events = t(d$grp[d$died]),
        estimate_pct_at_report_time = sapply(
          levels(d$grp),
          function(g) round(km_estimate_at(d$OS_time[d$grp == g],
                                           d$died[d$grp == g], report_t), 1)
        )
      )
      save_csv(surv_out, "T315I_05_Survival_T315Iplus_vs_minus.csv")
      message("\n>>> Log-rank test: chi-sq = ", round(logr$chisq, 3),
              ", p = ", format.pval(p_lr, digits = 3), "\n")
      message(">>> Survival estimates refer to ", round(report_t, 1),
              " months (60 months where follow-up allowed).\n")
      print(surv_out)
    }, error = function(e) {
      message(">>> Survival comparison could not be computed: ", conditionMessage(e))
    })
  }
}


# ---------------------------------------------------------------------------
# 9. NARRATIVE SUMMARY
# ---------------------------------------------------------------------------
cat("\n=============================================================\n")
cat(" SUMMARY (for the thesis / panel response)\n")
cat("=============================================================\n")
cat(paste0(
  "\nT315I mutation: documented in ", n_tested, " of ", nrow(df),
  " patients (", round(100 * n_tested / nrow(df), 1), "%). ",
  if (n_tested > 0) paste0("Among those tested, T315I was detected in ", n_pos,
                           " patients (", round(100 * n_pos / n_tested, 1), "%).") else "",
  " Because the dataset records only T315I, it is the only mutation whose",
  " frequency can be reported: detected in ", n_pos, " of ", nrow(df),
  " patients (", round(100 * n_pos / nrow(df), 1), "% of the cohort).\n\n",
  "Therapy of T315I-positive patients: of ", nrow(t3_pos), " T315I+ patient(s), ",
  sum(t3_pos$tki == "ponatinib", na.rm = TRUE),
  " were recorded on ponatinib (the TKI active against T315I) and ",
  sum(!t3_pos$tki %in% c("ponatinib", NA), na.rm = TRUE),
  " on TKIs without proven activity against T315I. The per-patient list with",
  " phase, survival status and cause of death is in T315I_03_Positive_Patients.csv.",
  " Limitation: only the CURRENT TKI is recorded, so recorded therapy may not",
  " reflect the TKI in use at the time the mutation was detected.\n"
))

cat("\n>>> Done. All tables are in: ", out_dir, "\n", sep = "")
