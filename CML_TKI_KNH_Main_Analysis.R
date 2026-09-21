# =============================================================================
#  MAIN ANALYSIS SCRIPT (REVISED - PANEL COMMENT VERSION)
#  Chronic Myeloid Leukemia (CML) patients managed with Tyrosine Kinase
#  Inhibitors (TKIs) at Kenyatta National Hospital (KNH)
#  Thesis project - KNH CML/TKI cohort
# =============================================================================
#
#  RESEARCH OBJECTIVES (reordered per panel comment vii:
#  treatment patterns -> survival -> predictors -> ADRs)
#
#   Objective 1: To describe TKI treatment patterns among CML patients
#                managed with TKIs at KNH, including drug selection, treatment
#                modifications, mutation testing, and evidence of sequential
#                TKI use.
#   Objective 2: To assess five-year overall survival (OS) and
#                progression-free survival (PFS) among CML patients managed
#                with TKIs at KNH.
#   Objective 3: To investigate predictors of survival outcomes (OS and PFS)
#                among CML patients managed with TKIs at KNH.
#   Objective 4: To determine the prevalence of adverse drug reactions (ADRs)
#                associated with TKIs, including severe (grade 3+) ADRs and
#                drug-specific ADR profiles.
#
#  PANEL COMMENT -> WHERE ADDRESSED
#   (i)    T315I mutation & therapy used           -> Section 4.7  (O1_11..13)
#   (ii)   Generation-outcome discrepancy          -> Sections 4.2 + 5.5
#                                                    (phase x group cross-tab,
#                                                    SUPPLEMENTARY 3-group KM
#                                                    slide table; the PRIMARY
#                                                    survival comparison is the
#                                                    2-group first-line vs
#                                                    subsequent-line LOT)
#   (iii)  Results as tables, not prose            -> all O*_0* tables are
#                                                     formatted for slides
#   (iv)   Grade >=3 ADR prevalence & distribution -> Section 7.2  (O4_04)
#   (v)    Every HR reported with 95% CI and p     -> Section 6 (Report column
#                                                     in all Cox tables)
#   (vi)   Full Cox variable table + EPV control   -> Section 6 (O3_02, EPV
#                                                     guard, Schoenfeld test)
#   (vii)  Reordering of objectives                -> overall file structure
#   (viii) Sequential TKI use                      -> Section 4.5  (O1_08/09)
#   (ix)   ADR distribution by individual drug     -> Section 7.3  (O4_05)
#   (x)    Line of therapy (first-line vs          -> Section 4.5b (O1_14) +
#            subsequent-line)
#                                                    PRIMARY survival
#                                                    comparison 5.1-5.4
#                                                    (O2_02..O2_08) + Table 1
#                                                    + Cox (2-group exposure)
#   (xi)   Grade 3+ ADRs -> OS/PFS: univariable
#                            + forced into the multivariable model
#                                                    -> Section 6 (O3_08)
#
#  IMPORTANT - DATA LIMITATIONS (state in the thesis):
#   * The dataset records only the CURRENT/primary TKI per patient (a single
#     "TKI Used" field), not a dated longitudinal treatment history.
#     The 1G/2G/3G labels are RECORDED TKI GROUPS, not verified treatment
#     lines: a 2G TKI is not automatically second-line, and ponatinib is not
#     automatically third-line. They must NOT be presented as a confirmed
#     treatment flow.
#   * LINE OF THERAPY (LOT) VARIABLE - PRIMARY SURVIVAL GROUPING (agreed
#     2-group design): First-line TKI = the first TKI the patient received for
#     CML; Subsequent-line TKI = any TKI given after the first TKI (2nd, 3rd,
#     ...). It is derived per patient from a line-of-therapy or first/initial-
#     TKI column if present, otherwise from the imatinib exposure column,
#     otherwise from the study assumption that every patient started imatinib
#     (see O1_14c_LOT_Derivation_Patients.csv). With the current data:
#     imatinib -> First-line; dasatinib/nilotinib/bosutinib/ponatinib ->
#     Subsequent-line. The ponatinib group has only n = 5 - far too few for any
#     meaningful separate survival analysis - so all 2G/3G TKIs are pooled as
#     "subsequent line" and the 2G TKIs are NOT analysed separately either
#     (their per-drug numbers are also too small). The 1G/2G/3G KM breakdown
#     (O2_10..O2_15) is retained as SUPPLEMENTARY descriptive material only.
#   * A BLANK modification field means no modification was recorded - it does
#     NOT mean the patient never switched. "Dose switch" should be checked
#     against the study codebook before being interpreted as a TKI change;
#     it is treated here as a documented TKI modification only.
#   * Sequential TKI use is inferred indirectly (cross-referencing "TKI Used"
#     with the blast-phase treatment field); no ordered, dated three-drug
#     sequence can be established from this dataset.
#   * Only T315I was captured; other BCR-ABL1 mutations (E255K, Y253H, F317L,
#     ...) were not assessed.
#   * ADRs are recorded for the current TKI (cross-sectional); grading is
#     absent in a substantial fraction of the cohort, so severe-ADR
#     prevalence is likely underestimated.
#   * Five-year estimates use 60 months. If maximum observed follow-up is
#     < 60 months, estimates are reported at the longest follow-up and
#     flagged Five_year_estimable = FALSE.
#   * Percentages for treatment-flow tables are computed WITHIN each recorded
#     TKI group (each group = 100%), per thesis reporting rules.
#
#  HOW TO RUN
#   1. Open this file in RStudio (R >= 4.0).
#   2. Source the file (or run it top to bottom).
#   3. Select your data CSV in the file dialog
#      (default: "am at this point 2.csv").
#   4. All tables (.csv) and figures (.png) are written to an "output" folder
#      next to the data file.
#
#  EXPECTED COLUMNS (exact names or close variants - auto-detected):
#    TKI_Used / TKI Used, Survival_Status, Event,
#    Duration_of_Survival_Months / OS_Time_Months, PFS_Event, PFS_Time_Months,
#    Age, Gender / Sex, Phase_of_CML, T315I_Mutation,
#    Any_ADRs / Any ADRs, Grade_3_ADR, Grade_3_ADR_Specified,
#    Name_of_ADR, Toxicity_1 ... Toxicity_4,
#    TKI_Treatment_Modification, TKI_Modification_Reason,
#    Disease_Progression, Disease_Relapse, Any_Comorbidity,
#    BCR_ABL_Monitoring, Blast_Phase_Treatment, Cause_of_Death,
#    Study_ID / Study ID
#
# =============================================================================

# -----------------------------------------------------------------------------
# 0. PACKAGES AND OPTIONS
# -----------------------------------------------------------------------------
required_packages <- c("dplyr", "readr", "tidyr", "stringr", "tibble",
                       "purrr", "survival", "survminer", "ggplot2", "broom")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Package '", pkg, "' is not installed. Installing now...")
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

options(stringsAsFactors = FALSE, width = 110)

T_FIVE_YEARS <- 60   # months; "five-year" = 60 months
VAR_ENTRY_P  <- 0.20 # univariable entry threshold for the multivariable model
EPV_MIN      <- 10   # recommended minimum events per covariate

# Registry of everything saved to disk (printed at the end).
saved_files <- character(0)

save_csv <- function(x, name) {
  path <- file.path(out_dir, name)
  write.csv(x, path, row.names = FALSE)
  saved_files <<- c(saved_files, path)
  message(">>> saved table:  ", name)
  invisible(x)
}

save_fig <- function(p, name, w = 8, h = 6) {
  path <- file.path(out_dir, name)
  ok <- tryCatch({
    ggsave(path, plot = p, width = w, height = h, dpi = 300)
    TRUE
  }, error = function(e) {
    message(">>> could not save figure '", name, "': ", conditionMessage(e))
    FALSE
  })
  if (ok) {
    saved_files <<- c(saved_files, path)
    message(">>> saved figure: ", name)
  }
  invisible(ok)
}

# -----------------------------------------------------------------------------
# HELPER FUNCTIONS
# -----------------------------------------------------------------------------

# Find a column in df: exact names first, then regex patterns (ignore.case).
find_col <- function(df, exact = character(), patterns = character()) {
  nms <- names(df)
  hit <- nms[nms %in% exact]
  if (length(hit) == 0L && length(patterns) > 0L) {
    hit <- nms[grepl(paste0(patterns, collapse = "|"), nms, ignore.case = TRUE)]
  }
  if (length(hit) >= 1L) hit[1L] else NULL
}

as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# Convert a status column to logical: TRUE = event (death), FALSE = alive.
status_01 <- function(x, event_words = c("deceased", "dead", "died", "death", "1")) {
  v <- tolower(trimws(as.character(x)))
  v[is.na(v) | v == ""] <- NA_character_
  is_event <- v %in% event_words
  is_alive <- v %in% c("alive", "living", "censored", "surviving", "0", "no")
  out <- is_event & !is_alive
  out[is.na(v)] <- NA
  out
}

# Convert a Yes/No-style column to a logical flag (NA preserved).
to_flag <- function(x, positives = c("yes", "y", "1", "true", "present")) {
  v <- tolower(trimws(as.character(x)))
  v[is.na(v) | v == ""] <- NA_character_
  v %in% positives
}

# 95% CI for a proportion (Wald; rule-of-three at the boundaries).
prop_ci <- function(x, n) {
  if (is.na(x) || is.na(n) || n < 1) return(c(NA_real_, NA_real_))
  if (x == 0)  return(c(0, min(1, 3 / n)))
  if (x >= n)  return(c(max(0, 1 - 3 / n), 1))
  ci <- tryCatch(
    suppressWarnings(prop.test(x, n, correct = FALSE)$conf.int),
    warning = function(w) c(NA_real_, NA_real_),
    error   = function(e) c(NA_real_, NA_real_)
  )
  if (length(ci) != 2 || all(is.na(ci))) c(NA_real_, NA_real_) else ci
}

# p-value for a contingency table (Fisher when small, else chi-square).
cat_pvalue <- function(tab) {
  if (is.null(dim(tab)) || nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  tab <- tab[rowSums(tab) > 0, , drop = FALSE]
  tab <- tab[, colSums(tab) > 0, drop = FALSE]
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  if (nrow(tab) == 2 && ncol(tab) == 2) {
    ft <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
    if (!is.na(ft)) return(ft)
  }
  if (any(tab < 5)) {
    ft <- tryCatch(
      suppressWarnings(fisher.test(tab, simulate.p.value = TRUE, B = 2000)$p.value),
      error = function(e) NA_real_
    )
    if (!is.na(ft)) return(ft)
  }
  tryCatch(suppressWarnings(chisq.test(tab, correct = FALSE)$p.value),
           error = function(e) NA_real_)
}

# Numeric summary: "median (Q1 - Q3)".
fmt_median_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  q <- quantile(x, c(0.25, 0.75), names = FALSE)
  sprintf("%.1f (%.1f - %.1f)", median(x), q[1], q[2])
}

# Factor summary: "level: n (%) | level: n (%) | ...".
fmt_factor_dist <- function(x) {
  if (all(is.na(x))) return("NA")
  tab <- table(factor(x))
  n_known <- sum(tab)
  if (n_known == 0) return("NA")
  paste(sprintf("%s: %d (%.1f%%)", names(tab), tab, 100 * tab / n_known),
        collapse = " | ")
}

# Logical-flag summary: "n (%)".
fmt_pct <- function(x) {
  if (all(is.na(x))) return("NA")
  n_known <- sum(!is.na(x))
  if (n_known == 0) return("NA")
  n_yes <- sum(x, na.rm = TRUE)
  sprintf("%d (%.1f%%)", n_yes, 100 * n_yes / n_known)
}

# Kruskal-Wallis p-value for a numeric variable across groups.
pval_numeric <- function(x, g) {
  d <- data.frame(x = x, g = g)
  d <- d[complete.cases(d), , drop = FALSE]
  if (nrow(d) < 3 || length(unique(d$g)) < 2) return(NA_real_)
  tryCatch(suppressWarnings(kruskal.test(x ~ g, data = d)$p.value),
           error = function(e) NA_real_)
}

# Categorical p-value (chi-square / Fisher) across groups.
pval_categorical <- function(x, g) {
  d <- data.frame(x = x, g = g)
  d <- d[complete.cases(d), , drop = FALSE]
  if (nrow(d) < 3 || length(unique(d$g)) < 2) return(NA_real_)
  cat_pvalue(table(d$g, d$x))
}

# Restricted mean survival time over [0, tau] from a one-group survfit
# (exact step-function integration of the KM curve).
rmst_step <- function(fit, tau) {
  if (is.null(fit) || length(fit$time) == 0) return(NA_real_)
  t <- c(0, fit$time)
  s <- c(1, fit$surv)
  if (t[length(t)] < tau) {
    t <- c(t, tau)
    s <- c(s, s[length(s)])
  }
  keep <- t <= tau
  t <- t[keep]
  s <- s[keep]
  if (length(t) < 2) return(NA_real_)
  sum(diff(t) * s[-length(s)])
}

# KM estimate at time t from a one-group survfit, with log-scale 95% CI.
km_estimate_at <- function(fit, t) {
  s <- summary(fit, times = t)
  sv <- s$survival[1]
  se <- s$std.err[1]
  lo <- NA_real_
  hi <- NA_real_
  if (!is.na(sv) && sv > 0 && is.finite(se)) {
    lo <- exp(log(sv) - 1.96 * se)
    hi <- exp(log(sv) + 1.96 * se)
  }
  list(estimate = sv, ci_lower = lo, ci_upper = hi)
}

# Median survival time from a one-group survfit; NA if not reached.
km_median <- function(fit) {
  if (is.null(fit) || length(fit$time) == 0) return(NA_real_)
  s <- fit$surv
  if (min(s) > 0.5) return(NA_real_)
  idx <- which(s <= 0.5)[1]
  fit$time[idx]
}

# One-row KM summary (n, events, median, estimate at min(60, max follow-up)).
outcome_summary <- function(fit, d, label) {
  n <- nrow(d)
  ne <- sum(d$Event)
  tmax <- max(d$Time)
  t_rep <- min(T_FIVE_YEARS, tmax)
  est <- km_estimate_at(fit, t_rep)
  tibble(
    Outcome = label,
    N = n,
    Events = ne,
    Event_rate_pct = round(100 * ne / n, 1),
    Median_time_months = round(median(d$Time, na.rm = TRUE), 1),
    Max_followup_months = round(tmax, 1),
    Estimate_at_months = t_rep,
    Estimate_pct = if (is.na(est$estimate)) NA_real_ else round(100 * est$estimate, 1),
    CI95_lower_pct = if (is.na(est$ci_lower)) NA_real_ else round(100 * est$ci_lower, 1),
    CI95_upper_pct = if (is.na(est$ci_upper)) NA_real_ else round(100 * est$ci_upper, 1),
    Five_year_estimable = tmax >= T_FIVE_YEARS,
    RMST_months = round(rmst_step(fit, t_rep), 1)
  )
}

# KM estimate at min(60, max follow-up), one row per level of `group_var`.
# group_var = "LOT" for the PRIMARY 2-group comparison (first-line vs
# subsequent-line); "TKI_Group" for the SUPPLEMENTARY 3-group breakdown.
km_by_group <- function(d, label, group_var = "LOT") {
  d <- d[!is.na(d[[group_var]]), , drop = FALSE]
  if (nrow(d) == 0) return(NULL)
  t_rep <- min(T_FIVE_YEARS, max(d$Time))
  rows <- list()
  for (g in levels(droplevels(d[[group_var]]))) {
    dg <- d[d[[group_var]] == g, , drop = FALSE]
    if (nrow(dg) == 0) next
    est <- km_estimate_at(survfit(Surv(Time, Event) ~ 1, data = dg), t_rep)
    rows[[g]] <- tibble(
      Outcome = label,
      Group = g,
      N = nrow(dg),
      Events = sum(dg$Event),
      Estimate_at_months = t_rep,
      Estimate_pct = if (is.na(est$estimate)) NA_real_ else round(100 * est$estimate, 1),
      CI95_lower_pct = if (is.na(est$ci_lower)) NA_real_ else round(100 * est$ci_lower, 1),
      CI95_upper_pct = if (is.na(est$ci_upper)) NA_real_ else round(100 * est$ci_upper, 1)
    )
  }
  if (length(rows) == 0) NULL else bind_rows(rows)
}

# Log-rank p-value of Time~Event across the levels of a group variable.
logrank_p <- function(d, group_var) {
  d2 <- d[!is.na(d[[group_var]]), , drop = FALSE]
  if (length(unique(d2[[group_var]])) < 2) return(NA_real_)
  f <- as.formula(paste("Surv(Time, Event) ~", group_var))
  sd <- survdiff(f, data = d2)
  dfc <- length(sd$n) - 1
  if (dfc < 1) return(NA_real_)
  pchisq(sd$chisq, df = dfc, lower.tail = FALSE)
}

# p-value formatted for reporting ("<0.001" style).
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

# Panel comment (v): standardized HR reporting string.
hr_report <- function(hr, lo, hi, p) {
  sprintf("HR %.2f (95%% CI %.2f - %.2f, p = %s)", hr, lo, hi, fmt_p(p))
}

# Extract TKI names mentioned in a free-text field (e.g. blast-phase
# treatment). Returns a collapsed string, "" if none, NA if no text.
ALL_TKIS <- c("imatinib", "dasatinib", "nilotinib", "bosutinib", "ponatinib")
tkis_in_text_str <- function(x) {
  v <- tolower(as.character(x))
  vapply(seq_along(v), function(i) {
    if (is.na(v[i]) || trimws(v[i]) == "") return(NA_character_)
    hit <- ALL_TKIS[vapply(ALL_TKIS, function(k) grepl(k, v[i], fixed = TRUE),
                           logical(1))]
    if (length(hit) == 0) "" else paste(hit, collapse = " + ")
  }, character(1))
}

# Light normalization of free-text ADR names for grouping.
norm_adr <- function(x) {
  v <- tolower(trimws(as.character(x)))
  v[is.na(v) | v == ""] <- NA_character_
  v <- ifelse(grepl("diarrhoea", v) & !grepl("diarrhea", v),
              sub("diarrhoea", "diarrhea", v), v)
  v <- ifelse(grepl("^(elevated |elev )?(ast|alt)$", v) |
                grepl("hepatit", v),
              "hepatotoxicity (elevated ast/alt)", v)
  v
}

# Known class-signal label for drug-specific ADR profiles (panel comment ix).
class_signal <- function(adr_type) {
  v <- tolower(adr_type)
  if (grepl("pleural", v))        return("Dasatinib signature (pleural effusion)")
  if (grepl("pericard", v))       return("Dasatinib signature (pericardial effusion)")
  if (grepl("gi bleed|gastrointest.*bleed", v)) return("Dasatinib signal (GI bleeding)")
  if (grepl("hepat|ast|alt", v))  return("Hepatotoxicity (nilotinib/ponatinib signal)")
  if (grepl("bilirubin", v))      return("Nilotinib signal")
  if (grepl("thrombocytopen|myelo|pancyt|neutropen|anaem|anemia", v))
    return("Haematological (class effect)")
  if (grepl("rash|dermat|prurit", v)) return("Dermatological (imatinib signal)")
  if (grepl("diarr|vomit|nausea|gastro", v)) return("GI (bosutinib signal)")
  ""
}

# Univariate + backward-stepwise multivariate Cox models for one endpoint.
# `data` must contain Time, Event and (where available) the candidate
# covariates. Returns univ / multiv / zph tables, the fitted model, and the
# number of events and covariates used (for EPV reporting).
run_cox_models <- function(data, force_vars = character(0)) {
  # TKI exposure = the TWO-GROUP line of therapy (LOT), not the 3-level
  # recorded TKI group: the ponatinib (3G) level has only n = 5 (its hazard
  # ratio is not reliably estimable) and per-drug 2G numbers are also small,
  # so all 2G/3G TKIs are pooled as "subsequent line" (agreed 2-group design).
  candidate_names <- c("Phase_group", "Age_group", "Sex", "LOT",
                       "flag_T315I", "flag_ADR_any", "flag_ADR_sev",
                       "flag_Progression", "flag_Relapse", "flag_Comorbidity",
                       "BCR_monitored", "mod_switch")
  var_labels <- c(
    Phase_group      = "CML phase (recorded)",
    Age_group        = "Age group (<40 / 40-59 / >=60 years)",
    Sex              = "Sex",
    LOT              = "Line of therapy (first TKI received for CML vs any TKI after it)",
    flag_T315I       = "T315I mutation",
    flag_ADR_any     = "Any ADR",
    flag_ADR_sev     = "Grade 3+ ADR",
    flag_Progression = "Disease progression",
    flag_Relapse     = "Disease relapse",
    flag_Comorbidity = "Any comorbidity",
    BCR_monitored    = "BCR-ABL monitoring done",
    mod_switch       = "Documented TKI dose switch"
  )
  candidates <- candidate_names[candidate_names %in% names(data)]

  # Keep only variables with usable variation and enough data.
  keep <- character(0)
  for (v in candidates) {
    x <- data[[v]]
    ok <- sum(!is.na(x)) >= 10
    if (is.logical(x)) {
      ok <- ok && length(unique(x[!is.na(x)])) >= 2
    } else if (is.numeric(x)) {
      ok <- ok && !is.na(stats::sd(x, na.rm = TRUE)) && stats::sd(x, na.rm = TRUE) > 0
    } else {
      ok <- ok && length(unique(x[!is.na(x)])) >= 2
    }
    if (ok) keep <- c(keep, v)
  }

  n_events <- sum(data$Event, na.rm = TRUE)

  # --- univariate ----------------------------------------------------------
  univ_list <- list()
  for (v in keep) {
    d2 <- data[!is.na(data$Time) & !is.na(data$Event) & !is.na(data[[v]]), , drop = FALSE]
    if (nrow(d2) < 10) next
    f <- as.formula(paste("Surv(Time, Event) ~", v))
    fit <- tryCatch(coxph(f, data = d2),
                    error = function(e) {
                      message("    coxph failed for '", v, "': ", conditionMessage(e))
                      NULL
                    })
    if (is.null(fit)) next
    tt <- broom::tidy(fit, conf.int = TRUE, exponent = TRUE)
    if (!"HR" %in% names(tt) && "estimate" %in% names(tt)) tt$HR <- tt$estimate
    univ_list[[v]] <- tt %>%
      filter(term != "(Intercept)") %>%
      mutate(Covariate = v, Label = var_labels[v],
             HR_rounded = round(HR, 2),
             P_report = map_chr(p.value, fmt_p),
             Report = map2_chr(HR, p.value, function(h, p) hr_report(h, 0, 0, p)))
  }
  # fix Report (needs CI): rebuild properly below
  if (length(univ_list) > 0) {
    univ_raw <- bind_rows(univ_list)
    univ <- univ_raw %>%
      mutate(Report = mapply(hr_report, HR, conf.low, conf.high, p.value))
  } else {
    univ <- tibble()
  }

  # --- multivariate (variables with univariate p < VAR_ENTRY_P) -------------
  multiv <- tibble()
  zph_tbl <- tibble()
  n_used <- NA_integer_
  n_cov_final <- 0L
  cox_fit <- NULL
  sel <- unique(univ$Covariate[univ$p.value < VAR_ENTRY_P])
  if (length(sel) == 0) sel <- keep

  # Panel requirement (xi): Grade 3+ ADR status and line of therapy must
  # appear in the multivariable model even if their univariable p >= 0.20.
  sel <- union(sel, intersect(force_vars, keep))

  # EPV guard (panel comment vi): limit covariates to ~events/EPV_MIN.
  # Forced covariates are retained; remaining slots go to the most
  # significant univariable p-values.
  max_cov <- max(2L, as.integer(floor(n_events / EPV_MIN)))
  if (length(sel) > max_cov) {
    forced <- intersect(force_vars, sel)
    pv <- aggregate(p.value ~ Covariate, data = univ, FUN = min)
    ranked <- pv$Covariate[order(pv$p.value, na.last = TRUE)]
    rest <- ranked[ranked %in% setdiff(sel, forced)]
    sel <- c(forced, head(rest, max(0L, max_cov - length(forced))))
    message("    EPV guard: multivariable model limited to ", max_cov,
            " covariates (events = ", n_events, "; ", EPV_MIN,
            " events per covariate); retained forced covariates: ",
            if (length(forced) > 0) paste(forced, collapse = ", ") else "none", ".")
  }

  if (length(sel) >= 1) {
    f_full <- as.formula(paste("Surv(Time, Event) ~", paste(sel, collapse = " + ")))
    d_full <- data[!is.na(data$Time) & !is.na(data$Event), , drop = FALSE]
    full_fit <- tryCatch(coxph(f_full, data = d_full),
                         error = function(e) {
                           message("    multivariate Cox model failed: ", conditionMessage(e))
                           NULL
                         })
    if (!is.null(full_fit)) {
      cox_fit <- tryCatch(step(full_fit, direction = "backward", trace = FALSE),
                          error = function(e) full_fit)
      n_used <- as.integer(summary(cox_fit)$n)
      mt <- broom::tidy(cox_fit, conf.int = TRUE, exponent = TRUE)
      if (!"HR" %in% names(mt) && "estimate" %in% names(mt)) mt$HR <- mt$estimate
      term_var <- vapply(mt$term, function(term) {
        if (term == "(Intercept)") return("Intercept")
        hits <- sel[vapply(sel, function(v) startsWith(term, v), logical(1))]
        if (length(hits) >= 1) hits[1] else NA_character_
      }, character(1))
      mt$Level <- vapply(mt$term, function(term) {
        if (term == "(Intercept)") return("")
        sub(paste0("^(Intercept|", paste(sel, collapse = "|"), ")"), "", term)
      }, character(1))
      multiv <- mt %>%
        mutate(Covariate = term_var,
               Label = ifelse(term_var == "Intercept", "Intercept (baseline hazard)",
                              var_labels[term_var]),
               HR_rounded = round(HR, 2),
               CI_report = paste0(round(conf.low, 2), " - ", round(conf.high, 2)),
               P_report = map_chr(p.value, fmt_p),
               Report = mapply(hr_report, HR, conf.low, conf.high, p.value))
      n_cov_final <- length(setdiff(unique(multiv$Covariate), "Intercept"))
      z <- tryCatch(cox.zph(cox_fit, test = TRUE), error = function(e) NULL)
      if (!is.null(z)) {
        zph_tbl <- bind_rows(
          tibble(Test = "Global", Chi2 = z$global$chisq, df = z$global$df,
                 P_value = z$global[["p value"]]),
          tibble(Test = rownames(z$cov), Chi2 = z$cov[, "chi2"],
                 df = z$cov[, "df"], P_value = z$cov[, "p value"])
        )
      }
    }
  }
  list(univ = univ, multiv = multiv, zph = zph_tbl, n_used = n_used,
       cox_fit = cox_fit, n_events = n_events, n_cov_final = n_cov_final)
}

# =============================================================================
# 1. UPLOAD / SELECT THE DATA FILE
# =============================================================================
csv_path <- NULL

tryCatch(
  {
    message("\n>>> Please select your data CSV file in the window that opens...")
    selected <- file.choose()
    if (!is.na(selected) && file.exists(selected)) csv_path <- selected
  },
  error = function(e) message("File dialog unavailable: ", conditionMessage(e))
)

if (is.null(csv_path)) {
  message("\n>>> File dialog cancelled/unavailable. Type the CSV path and press Enter,")
  message(">>> or press Enter to use the default file.")
  typed <- trimws(readline(prompt = "CSV path: "))
  if (is.na(typed)) typed <- ""
  if (nzchar(typed) && file.exists(typed)) {
    csv_path <- typed
  } else if (nzchar(typed)) {
    stop("File not found at: ", typed)
  }
}

if (is.null(csv_path)) {
  default_path <- "am at this point 2.csv"
  if (file.exists(default_path)) {
    csv_path <- default_path
  } else {
    stop("No file selected and default file not found: ", default_path)
  }
}

# Guard: this script reads CSV files, not Excel workbooks.
if (tolower(tools::file_ext(csv_path)) %in% c("xlsx", "xls")) {
  stop("The selected file is an Excel workbook (", basename(csv_path), ").\n",
       "This script reads CSV files. In Excel use File > Save As > CSV (comma\n",
       "delimited) (.csv), then select that .csv file in the dialog.")
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

# Guard: the CML data file is wide (dozens of variables). A 1- or 2-column
# load almost always means the wrong file was selected (e.g. an exported
# Excel/zip part) or the file has no proper header row.
if (ncol(df) < 5 || nrow(df) < 5) {
  stop("Loaded file has ", nrow(df), " row(s) x ", ncol(df),
       " column(s) - it does not look like the CML dataset.\n",
       "Check that you selected the correct .csv file with a header row.")
}

if (!"Study_ID" %in% names(df) && "Study ID" %in% names(df)) {
  df <- df %>% rename(Study_ID = `Study ID`)
}

out_dir <- file.path(dirname(normalizePath(csv_path)), "output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
message(">>> Results will be saved to: ", out_dir)

# -----------------------------------------------------------------------------
# 2. COLUMN DISCOVERY (auto-detect expected columns)
# -----------------------------------------------------------------------------
toxicity_cols <- names(df)[grepl("^toxicity", names(df), ignore.case = TRUE)]

colmap <- list(
  tki         = find_col(df, c("TKI_Used", "TKI Used", "Current TKI", "TKI"),
                         patterns = c("^TKI[_ ]?Used", "^Current[_ ]?TKI", "^TKI$")),
  os_status   = find_col(df, c("Event", "Survival_Status", "Survival Status", "Status"),
                         patterns = c("Survival[_ ]?Status", "^Status$")),
  os_time     = find_col(df, c("Duration_of_Survival_Months", "Duration of Survival (Months)",
                               "OS_Time_Months", "OS Time (Months)", "OS (Months)",
                               "Overall_Survival_Months", "OS_Time"),
                         patterns = c("Duration[_ ]of[_ ]Survival", "OS[_ ]?Time",
                                      "Overall[_ ]?Survival[_ ]?(Time|Month)")),
  pfs_event   = find_col(df, c("PFS_Event", "PFS Event"),
                         patterns = c("^PFS[_ ]?Event")),
  pfs_time    = find_col(df, c("PFS_Time_Months", "Duration_of_PFS_Months",
                               "PFS Time (Months)", "PFS (Months)"),
                         patterns = c("PFS[_ ]?Time", "PFS[_ ]?Months?")),
  age         = find_col(df, c("Age", "age", "Age (years)", "Age at diagnosis"),
                         patterns = c("^Age")),
  sex         = find_col(df, c("Gender", "Sex", "gender", "sex"),
                         patterns = c("^Sex", "^Gender")),
  phase       = find_col(df, c("Phase_of_CML", "Phase of CML", "CML Phase", "CML_Phase",
                               "Phase"),
                         patterns = c("Phase")),
  t315i       = find_col(df, c("T315I_Mutation", "T315I Mutation", "T315I"),
                         patterns = c("T315I")),
  adr_any     = find_col(df, c("Any_ADRs", "Any ADRs", "Any ADR", "Any_ADR"),
                         patterns = c("Any[_ ]?ADR", "^ADR")),
  adr_sev     = find_col(df, c("Grade_3_ADR", "Grade 3 ADR", "Grade3_ADR", "Grade 3+ ADR"),
                         patterns = c("Grade[_ ]?3")),
  adr_sev_sp  = find_col(df, c("Grade_3_ADR_Specified", "Grade 3 ADR Specified"),
                         patterns = c("Grade[_ ]?3[_ ]?ADR[_ ]?Specif")),
  adr_name    = find_col(df, c("Name_of_ADR", "Name of ADR", "ADR_Name"),
                         patterns = c("Name[_ ]of[_ ]ADR")),
  relapse     = find_col(df, c("Disease_Relapse", "Disease Relapse", "Relapse"),
                         patterns = c("Relapse")),
  progression = find_col(df, c("Disease_Progression", "Disease Progression", "Progression"),
                         patterns = c("Progression")),
  comorbidity = find_col(df, c("Any_Comorbidity", "Any Comorbidity", "Comorbidity"),
                         patterns = c("Comorbid")),
  bcr_mon     = find_col(df, c("BCR_ABL_Monitoring", "BCR-ABL Monitoring",
                               "BCR ABL Monitoring", "BCR_ABL_Monitoring_Done"),
                         patterns = c("BCR[-_ ]?ABL[_ ]?Monitoring|^BCR[_ ]?Monitoring")),
  modif       = find_col(df, c("TKI_Treatment_Modification", "TKI Treatment Modification"),
                         patterns = c("Treatment[_ ]?Modification", "^Modification")),
  reason      = find_col(df, c("TKI_Modification_Reason", "TKI Modification Reason"),
                         patterns = c("Modification[_ ]?Reason")),
  blast_tx    = find_col(df, c("Blast_Phase_Treatment", "Blast Phase Treatment",
                               "Blast treatment"),
                         patterns = c("Blast[_ ]Phase[_ ]Treatment|^Blast[_ ]?Treatment")),
  cause_death = find_col(df, c("Cause of Death", "Cause_of_Death"),
                         patterns = c("Cause[_ ]?of[_ ]?Death")),
  study_id    = find_col(df, c("Study_ID", "Study ID"), patterns = c("^Study[_ ]?ID"))
)
colmap$toxicity_cols <- toxicity_cols

message("\n>>> Column mapping:")
for (nm in names(colmap)) {
  if (nm == "toxicity_cols") {
    if (length(toxicity_cols) > 0) {
      message("      ", strrep("-", 16), "toxicity_cols : ",
              paste(toxicity_cols, collapse = ", "))
    } else {
      message("      ", strrep("-", 16), "toxicity_cols : NOT FOUND")
    }
    next
  }
  if (is.null(colmap[[nm]])) {
    message("      ", strrep("-", 16), nm, " : NOT FOUND")
  } else {
    message("      ", strrep("-", 16), nm, " : '", colmap[[nm]], "'")
  }
}

if (is.null(colmap$tki)) {
  stop("Column 'TKI Used / TKI_Used' not found - cannot describe the TKI groups. ",
       "Check that the CSV has a column listing the current/primary TKI.")
}

# =============================================================================
# 3. DATA PREPARATION
# =============================================================================
message("\n>>> Preparing variables...")

# --- current/primary TKI and recorded TKI group (1G/2G/3G) ------------------
# NOTE: these are RECORDED TKI GROUPS, not verified treatment lines.
df$tki_clean <- tolower(str_trim(as.character(df[[colmap$tki]])))
df$tki_clean[is.na(df$tki_clean) | df$tki_clean == ""] <- NA_character_

df$TKI_Group <- factor(
  case_when(
    df$tki_clean == "imatinib" ~ "1G",
    df$tki_clean %in% c("dasatinib", "nilotinib", "bosutinib") ~ "2G",
    df$tki_clean == "ponatinib" ~ "3G",
    TRUE ~ NA_character_
  ),
  levels = c("1G", "2G", "3G")
)
# NOTE: the per-line data set (df_groups) is built at the END of this section,
# after ALL derived columns (blast_tx_*, seq_inferred, ...) exist, so the
# snapshot used by every later section is complete.

group_desc <- c("1G" = "imatinib",
                "2G" = "dasatinib / nilotinib / bosutinib",
                "3G" = "ponatinib")

# --- time-to-event variables -------------------------------------------------
have_os  <- !is.null(colmap$os_status) && !is.null(colmap$os_time)
have_pfs <- !is.null(colmap$pfs_event) && !is.null(colmap$pfs_time)

if (have_os) {
  df$OS_Status <- as.integer(status_01(df[[colmap$os_status]]))
  df$OS_Time   <- as_num(df[[colmap$os_time]])
} else {
  df$OS_Status <- NA
  df$OS_Time   <- NA
  message("WARNING: OS columns not found; Objective 2 (OS part) will be limited.")
}

if (have_pfs) {
  df$PFS_Event <- as.integer(status_01(df[[colmap$pfs_event]]))
  df$PFS_Time  <- as_num(df[[colmap$pfs_time]])
} else {
  df$PFS_Event <- NA
  df$PFS_Time  <- NA
  message("WARNING: PFS columns not found; Objective 2 (PFS part) will be limited.")
}

# --- death flag (for flow tables) ---------------------------------------------
df$died <- if (!is.null(colmap$os_status)) status_01(df[[colmap$os_status]]) else FALSE

# --- documented treatment modifications -----------------------------------------
if (!is.null(colmap$modif)) {
  modif_txt <- tolower(str_trim(coalesce(as.character(df[[colmap$modif]]), "")))
  df$mod_switch       <- str_detect(modif_txt, "switch")
  df$mod_interruption <- str_detect(modif_txt, "interruption")
  df$mod_dose_reduction <- str_detect(modif_txt, "reduction")
  df$any_modification <- str_detect(modif_txt, "switch|interruption|reduction|discontin")
} else {
  df$mod_switch       <- FALSE
  df$mod_interruption <- FALSE
  df$mod_dose_reduction <- FALSE
  df$any_modification <- FALSE
}

# --- clinical factor flags -------------------------------------------------------
flag_map <- c(flag_T315I = colmap$t315i, flag_ADR_any = colmap$adr_any,
              flag_ADR_sev = colmap$adr_sev, flag_Relapse = colmap$relapse,
              flag_Progression = colmap$progression,
              flag_Comorbidity = colmap$comorbidity)
for (vname in names(flag_map)) {
  src <- flag_map[[vname]]
  df[[vname]] <- if (is.null(src)) NA else to_flag(df[[src]])
}

# --- age / age group / sex / phase / phase group ----------------------------------
df$Age_num <- if (!is.null(colmap$age)) as_num(df[[colmap$age]]) else NA
df$Age_group <- factor(
  cut(df$Age_num, breaks = c(-Inf, 39, 59, Inf),
      labels = c("<40", "40-59", ">=60"))
)

if (!is.null(colmap$sex)) {
  s <- tolower(str_trim(as.character(df[[colmap$sex]])))
  s[s %in% c("m", "male")] <- "Male"
  s[s %in% c("f", "female")] <- "Female"
  df$Sex <- factor(s, levels = c("Female", "Male"))
} else {
  df$Sex <- NULL
}

if (!is.null(colmap$phase)) {
  ph <- tolower(str_trim(as.character(df[[colmap$phase]])))
  df$Phase_raw <- factor(trimws(as.character(df[[colmap$phase]])))
  df$Phase_group <- factor(
    ifelse(grepl("chron", ph), "Chronic",
           ifelse(grepl("accel", ph), "Accelerated",
                  ifelse(grepl("blast", ph), "Blast", NA))),
    levels = c("Chronic", "Accelerated", "Blast")
  )
} else {
  df$Phase_raw <- NULL
  df$Phase_group <- NULL
}

# --- BCR-ABL monitoring (proxy for quality of monitoring) --------------------------
if (!is.null(colmap$bcr_mon)) {
  b <- tolower(str_trim(as.character(df[[colmap$bcr_mon]])))
  b[is.na(b) | b == ""] <- NA_character_
  df$BCR_monitored <- factor(
    ifelse(b %in% c("not done", "not documented", "no", "0"), "No",
           ifelse(b %in% c("yes", "1", "done", "completed") |
                    grepl("mmr|dmr|major|complete|optimal|deep|mr[0-9]", b), "Yes",
                  NA_character_)),
    levels = c("No", "Yes")
  )
} else {
  df$BCR_monitored <- NULL
}

# --- T315I status (4-level text) + positive flag -----------------------------------
if (!is.null(colmap$t315i)) {
  t3 <- tolower(str_trim(as.character(df[[colmap$t315i]])))
  t3[t3 %in% c("", "nd", "not documented")] <- "Not documented"
  t3[t3 %in% c("na", "n/a", "not applicable")] <- "Not applicable"
  t3[t3 %in% c("yes", "y", "1", "positive")] <- "Yes"
  t3[t3 %in% c("no", "n", "0", "negative")] <- "No"
  df$T315I_cat <- factor(t3,
                         levels = c("Yes", "No", "Not applicable", "Not documented"))
  df$flag_T315I <- df$T315I_cat == "Yes"
} else {
  df$T315I_cat <- NULL
}

# --- modification reason: ADR-related? ---------------------------------------------
if (!is.null(colmap$reason)) {
  df$reason_txt <- tolower(str_trim(as.character(df[[colmap$reason]])))
  df$reason_adr_related <- str_detect(
    df$reason_txt,
    "adverse|\\badr\\b|intoler|toxic|side.effect|hypersens"
  )
} else {
  df$reason_txt <- NA
  df$reason_adr_related <- FALSE
}

# --- cause of death (cleaned text) ----------------------------------------------------
df$Cause_of_Death_clean <- if (!is.null(colmap$cause_death))
  str_trim(as.character(df[[colmap$cause_death]])) else NA

# --- blast-phase treatment: TKIs mentioned (sequential-use evidence) -----------------
df$blast_tx_raw <- if (!is.null(colmap$blast_tx))
  str_trim(as.character(df[[colmap$blast_tx]])) else NA
df$blast_tkis <- tkis_in_text_str(df$blast_tx_raw)

# TKIs in the blast-phase field that DIFFER from the recorded TKI.
df$blast_other_tki <- mapply(function(tki, txt) {
  if (is.na(txt) || trimws(txt) == "" || tolower(txt) %in%
      c("not applicable", "not documented", "none")) return(NA_character_)
  pool <- if (is.na(tki)) ALL_TKIS else ALL_TKIS[ALL_TKIS != tki]
  others <- pool[vapply(pool, function(k) grepl(k, tolower(txt), fixed = TRUE),
                        logical(1))]
  if (length(others) == 0) NA_character_ else paste(others, collapse = " + ")
}, df$tki_clean, df$blast_tx_raw)

# Evidence of sequential/different TKI exposure (panel comment viii).
df$seq_evidence <- case_when(
  !is.na(df$blast_other_tki) ~ "Different TKI documented in blast-phase treatment",
  df$TKI_Group == "3G"       ~ "Implied (ponatinib is almost never first-line)",
  TRUE ~ NA_character_
)

df$seq_inferred <- mapply(function(other, tki) {
  if (!is.na(other) && !is.na(tki)) {
    return(paste(unique(c(other, tki)), collapse = " -> "))
  }
  NA_character_
}, df$blast_other_tki, df$tki_clean)
idx3 <- which(!is.na(df$TKI_Group) & df$TKI_Group == "3G" & is.na(df$seq_inferred))
df$seq_inferred[idx3] <- "prior TKI(s) -> ponatinib"

# --- Line of therapy (LOT): First-line vs Subsequent-line TKI ------------------------
# Definition (agreed):
#   First-line TKI      = the first TKI the patient received for CML
#   Subsequent-line TKI = any TKI given AFTER the first TKI (2nd-line,
#                         3rd-line, ...)
#
# Derivation, in priority order (the dataset records only the CURRENT/primary
# TKI per patient):
#   1. a "Line of Therapy" column, if the dataset has one (line labels such as
#      1/First/2L/Subsequent, or a named TKI);
#   2. a "First/Initial TKI" column, if present;
#   3. the imatinib exposure column - a patient recorded as "not treated with
#      imatinib" never received imatinib, so the recorded (current) TKI is
#      their FIRST TKI;
#   4. study assumption (see Fix_Data_IMatinib_Assumption.R): EVERY patient
#      started imatinib as the first TKI -> first TKI = imatinib, so the
#      recorded TKI is subsequent-line unless it is imatinib itself.
# LOT = First-line when the first TKI equals the recorded current TKI,
# otherwise Subsequent-line. With the current data this classifies: imatinib
# -> First-line; dasatinib / nilotinib / bosutinib / ponatinib ->
# Subsequent-line (prior imatinib). LOT is also the PRIMARY two-group
# survival grouping (the ponatinib group, n = 5, is too small for separate
# analysis): Section 5 (O2_02..O2_08), Table 1, Cox models.
#
# LIMITATION: a documented switch away from and back to imatinib, or a stale
# "current TKI" field, cannot be resolved from a single cross-sectional TKI
# record. Every derivation is reported per patient in
# O1_14c_LOT_Derivation_Patients.csv for verification.

# --- Priority 1: direct "Line of Therapy" column ---------------------------------------
col_lot_direct <- find_col(df,
  exact = c("Line_of_Therapy", "Line of Therapy", "Line of TKI therapy", "LOT"),
  patterns = c("^Line[_ ]?of[_ ]?(TKI )?Therapy", "^LOT$"))
lot_direct <- rep(NA_character_, nrow(df))
first_tki_direct <- rep(NA_character_, nrow(df))
if (!is.null(col_lot_direct)) {
  v <- tolower(trimws(as.character(df[[col_lot_direct]])))
  first_like <- v %in% c("1", "1l", "first", "first-line", "first line",
                         "firstline", "line 1")
  subseq_like <- v %in% c("2", "3", "4", "5", "2l", "3l", "4l", "5l", "2nd",
                          "3rd", "second", "third", "fourth", "later",
                          "later-line", "subsequent", "second-line",
                          "third-line", "subsequent-line", ">=2", ">= 2")
  lot_direct[first_like] <- "First-line TKI"
  lot_direct[subseq_like] <- "Subsequent-line TKI"
  # cells that name a TKI instead: the named drug is the first TKI
  hit <- tkis_in_text_str(v)
  named <- !is.na(hit) & hit != ""
  if (any(named)) {
    first_tki_direct[named] <- vapply(strsplit(hit[named], " \\+ "),
                                      function(p) p[1], character(1))
  }
}

# --- Priority 2: "First/Initial TKI" column ----------------------------------------------
col_first_tki <- find_col(df,
  exact = c("First_line_TKI", "First Line TKI", "First TKI", "Initial_TKI",
            "Initial TKI", "1st_line_TKI", "1L TKI", "Starting_TKI",
            "First_line_therapy", "First line therapy"),
  patterns = c("^First[_ ]?(line|Line)[_ ]?(TKI|therapy|treatment)",
               "^Initial[_ ]?(TKI|therapy|treatment)",
               "^First[_ ]?TKI", "^1L[_ ]?TKI", "^Starting[_ ]?TKI"))
first_tki_col2 <- rep(NA_character_, nrow(df))
if (!is.null(col_first_tki)) {
  hit2 <- tkis_in_text_str(tolower(str_trim(as.character(df[[col_first_tki]]))))
  named2 <- !is.na(hit2) & hit2 != ""
  if (any(named2)) {
    first_tki_col2[named2] <- vapply(strsplit(hit2[named2], " \\+ "),
                                     function(p) p[1], character(1))
  }
}

# --- Priority 3: imatinib exposure column -------------------------------------------------
im_cols <- names(df)[grepl("imatinib", names(df), ignore.case = TRUE)]
col_im_exposure <- NULL
if (length(im_cols) > 0) {
  exp_like <- im_cols[grepl("treat|use|expos|therapy|line|start|first|1l",
                            im_cols, ignore.case = TRUE)]
  col_im_exposure <- if (length(exp_like) > 0) exp_like[1] else im_cols[1]
}
# a duration/response/dose column is NOT an exposure indicator
if (!is.null(col_im_exposure) &&
    grepl("duration|months|response|mmr|cmr|relapse|dose|progression",
          col_im_exposure, ignore.case = TRUE)) {
  col_im_exposure <- NULL
}
im_exposure <- rep(NA, nrow(df))   # TRUE = imatinib documented, FALSE = none recorded
if (!is.null(col_im_exposure)) {
  ev <- tolower(trimws(as.character(df[[col_im_exposure]])))
  im_exposure[ev %in% c("yes", "y", "1", "true") |
               (grepl("^imatinib( |$)", ev) & !grepl("not|no|none", ev))] <- TRUE
  im_exposure[ev %in% c("no", "n", "0", "false", "none", "no imatinib") |
               grepl("not (treated|on|exposed)|never (treated|on|exposed)", ev)] <- FALSE
  # blank / unrecognised values stay NA (study assumption applies)
}

message("\n>>> Line of therapy (LOT) - data sources found:")
message("    Line-of-Therapy column     : ",
        if (is.null(col_lot_direct)) "not found" else col_lot_direct)
message("    First/Initial-TKI column   : ",
        if (is.null(col_first_tki)) "not found" else col_first_tki)
message("    Imatinib exposure column   : ",
        if (is.null(col_im_exposure)) "not found (study 1L-imatinib assumption used)"
        else col_im_exposure)

# --- assemble the first TKI per patient (priority 1 > 2 > 3 > assumption) ----------------
first_tki_final <- first_tki_direct
m <- is.na(first_tki_final)
first_tki_final[m] <- first_tki_col2[m]
m <- is.na(first_tki_final) & im_exposure == FALSE
first_tki_final[m] <- df$tki_clean[m]
first_tki_final[is.na(first_tki_final)] <- "imatinib"
df$first_tki <- first_tki_final

df$LOT_derivation <- case_when(
  !is.na(lot_direct) ~ "Line-of-Therapy column (line label given directly)",
  !is.na(first_tki_direct) ~ "Line-of-Therapy column (TKI named in the cell)",
  !is.na(first_tki_col2) ~ "First/Initial-TKI column",
  im_exposure == FALSE ~ "No imatinib exposure recorded -> recorded TKI is the first TKI",
  TRUE ~ "Study assumption: imatinib started as first TKI"
)

df$LOT <- factor(
  case_when(
    !is.na(lot_direct) ~ lot_direct,
    is.na(df$tki_clean) ~ NA_character_,
    df$first_tki != df$tki_clean ~ "Subsequent-line TKI",
    TRUE ~ "First-line TKI"
  ),
  levels = c("First-line TKI", "Subsequent-line TKI")
)
message(">>> Line of therapy: ",
        sum(df$LOT == "First-line TKI", na.rm = TRUE), " on FIRST-LINE TKI | ",
        sum(df$LOT == "Subsequent-line TKI", na.rm = TRUE),
        " on SUBSEQUENT-LINE TKI (any TKI given after the first TKI).")
message("    Derivation routes:\n")
print(table(df$LOT_derivation, useNA = "no"))
n_nonim_first <- sum(df$LOT == "First-line TKI" & !is.na(df$tki_clean) &
                       df$tki_clean != "imatinib", na.rm = TRUE)
if (n_nonim_first > 0) {
  message("    NOTE: ", n_nonim_first,
          " patient(s) are FIRST-LINE on a non-imatinib TKI - their recorded")
  message("    TKI is the first TKI they received for CML (check the audit file).")
}
# per-patient derivation table (built here, saved in Section 4.5b as O1_14c)
lot_audit <- df %>%
  select(any_of(c("Study_ID", "tki_clean", "first_tki", "LOT_derivation",
                  "mod_switch", "LOT"))) %>%
  rename(recorded_TKI = tki_clean, First_TKI = first_tki,
         Derivation = LOT_derivation,
         Documented_switch = mod_switch)
if (!is.null(col_im_exposure)) {
  lot_audit$Imatinib_exposure <- ifelse(im_exposure == TRUE, "Yes",
                                        ifelse(im_exposure == FALSE, "No",
                                               "Not documented"))
}

# --- group data set (AFTER all derived columns exist) -----------------------------
df_groups <- df %>% filter(!is.na(TKI_Group))
message(">>> Patients with a classifiable recorded TKI group: ", nrow(df_groups),
        " of ", nrow(df))

# --- analysis data sets (complete cases for each endpoint) -----------------------------
extra_cols <- intersect(
  c("Study_ID", "tki_clean", "TKI_Group", "LOT", "Age_num", "Age_group", "Sex",
    "Phase_raw", "Phase_group", "flag_T315I", "flag_ADR_any", "flag_ADR_sev",
    "flag_Relapse", "flag_Progression", "flag_Comorbidity", "BCR_monitored",
    "mod_switch", "mod_interruption", "any_modification"),
  names(df)
)

# NB: use select() + rename(), NOT transmute(all_of(...)) - a bare all_of()
# inside transmute() is not a column selector and fails with "must be size n
# or 1". select() is the correct context for tidyselect helpers.
os_data <- NULL
if (have_os) {
  os_data <- df %>%
    filter(!is.na(OS_Time), !is.na(OS_Status), OS_Time >= 0, OS_Status %in% 0:1) %>%
    select(all_of(c("OS_Time", "OS_Status", extra_cols))) %>%
    rename(Time = OS_Time, Event = OS_Status)
}

pfs_data <- NULL
if (have_pfs) {
  pfs_data <- df %>%
    filter(!is.na(PFS_Time), !is.na(PFS_Event), PFS_Time >= 0, PFS_Event %in% 0:1) %>%
    select(all_of(c("PFS_Time", "PFS_Event", extra_cols))) %>%
    rename(Time = PFS_Time, Event = PFS_Event)
}

message("\n>>> Data quality:")
message("    Total patients:                ", nrow(df))
message("    Classifiable recorded TKI grp: ", nrow(df_groups), " of ", nrow(df))
if (have_os)  message("    Complete OS data:              ", nrow(os_data))
if (have_pfs) message("    Complete PFS data:             ", nrow(pfs_data))

tki_dist <- df %>% count(tki_clean, name = "n", sort = TRUE)
message("\n>>> Recorded TKI distribution:\n")
print(tki_dist)

# -----------------------------------------------------------------------------
# Cohort summary (saved table)
# -----------------------------------------------------------------------------
cs_rows <- list(tibble(Variable = "Total patients (N)", Value = as.character(nrow(df))))
if (nrow(df_groups) > 0) {
  ld <- df_groups %>%
    count(TKI_Group, name = "n") %>%
    mutate(desc = group_desc[as.character(TKI_Group)],
           Variable = paste0("Recorded TKI group ", as.character(TKI_Group),
                             " (", desc, ")"),
           Value = paste0(n, " (", round(100 * n / nrow(df_groups), 1),
                          "% of classifiable)"))
  cs_rows <- c(cs_rows, list(ld %>% select(Variable, Value)))
}
if (have_os && nrow(os_data) > 0) {
  nd <- sum(os_data$Event)
  cs_rows <- c(cs_rows, list(tibble(
    Variable = "Deaths during OS follow-up",
    Value = paste0(nd, " (", round(100 * nd / nrow(os_data), 1), "%)"))))
}
if (have_pfs && nrow(pfs_data) > 0) {
  ne <- sum(pfs_data$Event)
  cs_rows <- c(cs_rows, list(tibble(
    Variable = "PFS events (progression or death)",
    Value = paste0(ne, " (", round(100 * ne / nrow(pfs_data), 1), "%)"))))
}
cohort_summary <- bind_rows(cs_rows)
save_csv(cohort_summary, "O0_Cohort_Summary.csv")

# =============================================================================
# 4. OBJECTIVE 1 - TKI TREATMENT PATTERNS
# =============================================================================
# (reordered to Objective 1 per panel comment vii)
message("\n>>> Objective 1: TKI treatment patterns...")

# --- 4.1 Recorded TKI group distribution (slide-ready) --------------------------
o1_01 <- df_groups %>%
  count(TKI_Group, name = "n") %>%
  mutate(Drugs = group_desc[as.character(TKI_Group)],
         Pct_of_cohort = round(100 * n / nrow(df), 1))
save_csv(o1_01, "O1_01_TKI_Group_Distribution.csv")
print(o1_01)

# --- 4.2 Phase distribution by recorded TKI group (panel comment ii) -------------
if (!is.null(colmap$phase)) {
  o1_02 <- df_groups %>%
    filter(!is.na(Phase_group)) %>%
    count(TKI_Group, Phase_group, name = "n") %>%
    group_by(TKI_Group) %>%
    mutate(Pct_of_group = round(100 * n / sum(n), 1)) %>%
    ungroup()
  p_phase_grp <- cat_pvalue(
    table(df_groups$TKI_Group, df_groups$Phase_group, useNA = "no"))
  o1_02_p <- tibble(
    Test = "Phase of CML vs recorded TKI group",
    P_value = p_phase_grp,
    P_report = fmt_p(p_phase_grp))
  save_csv(o1_02, "O1_02_Phase_by_TKI_Group.csv")
  save_csv(o1_02_p, "O1_02_Phase_by_TKI_Group_pvalue.csv")
  print(o1_02)
  print(o1_02_p)
} else {
  message("    Phase column not found - phase x group cross-tab skipped.")
}

# --- 4.3 Documented treatment modifications by group ------------------------------
o1_03 <- df_groups %>%
  count(TKI_Group, name = "n") %>%
  left_join(
    df_groups %>%
      group_by(TKI_Group) %>%
      summarise(n_switch = sum(mod_switch),
                n_interruption = sum(mod_interruption),
                n_dose_reduction = sum(mod_dose_reduction),
                n_any_mod = sum(any_modification),
                .groups = "drop"),
    by = "TKI_Group") %>%
  mutate(Pct_any_mod = round(100 * n_any_mod / n, 1))
save_csv(o1_03, "O1_03_Treatment_Modifications_by_Group.csv")
print(o1_03)

# --- 4.4 Documented flow within each recorded group (each group = 100%) -----------
# A blank modification field means no modification was recorded - it does NOT
# mean the patient never switched.
df_flow <- df_groups %>%
  mutate(Flow_group = case_when(
    mod_switch &  died ~ "Switched TKI - died after switch",
    mod_switch & !died ~ "Switched TKI - remained alive",
    !mod_switch & died ~ "Died without recorded switch",
    TRUE ~ "Maintained (alive, no recorded switch)"
  ),
    Flow_group = factor(Flow_group, levels = c(
      "Maintained (alive, no recorded switch)",
      "Switched TKI - remained alive",
      "Switched TKI - died after switch",
      "Died without recorded switch")))

o1_04 <- df_flow %>%
  count(TKI_Group, Flow_group, name = "n") %>%
  group_by(TKI_Group) %>%
  mutate(Percent_of_group = round(100 * n / sum(n), 1)) %>%
  ungroup()
save_csv(o1_04, "O1_04_Documented_Flow_by_Group.csv")
print(o1_04)

o1_05 <- df_flow %>%
  filter(mod_switch) %>%
  group_by(TKI_Group) %>%
  summarise(
    Switchers = n(),
    Died_after_switch = sum(died, na.rm = TRUE),
    Remained_alive = sum(!died, na.rm = TRUE),
    Pct_died_after_switch = if (n() > 0) round(100 * sum(died, na.rm = TRUE) / n(), 1) else NA_real_,
    .groups = "drop"
  )
save_csv(o1_05, "O1_05_Outcome_Among_Switchers.csv")
print(o1_05)

o1_06 <- df_flow %>%
  filter(mod_switch, died) %>%
  select(any_of(c("Study_ID", "tki_clean", "TKI_Group", "Phase_raw",
                  "reason_txt", "OS_Time"))) %>%
  rename(TKI = tki_clean, Modification_reason = reason_txt, OS_Time_Months = OS_Time)
save_csv(o1_06, "O1_06_Switched_Then_Died_Patients.csv")
print(o1_06)

if (!is.null(colmap$reason)) {
  o1_07 <- df_flow %>%
    filter(!is.na(reason_txt), reason_txt != "") %>%
    count(TKI_Group, Reason = reason_txt, name = "n") %>%
    group_by(TKI_Group) %>%
    mutate(Percent_of_group = round(100 * n / sum(n), 1)) %>%
    ungroup()
  save_csv(o1_07, "O1_07_Modification_Reasons_by_Group.csv")
  print(o1_07)
}

# --- 4.5 Sequential TKI use - indirect evidence (panel comment viii) ---------------
message("    Panel comment (viii): sequential TKI use...")

seq_detail <- df %>%
  filter(!is.na(seq_evidence)) %>%
  select(any_of(c("Study_ID", "Age_num", "Phase_raw", "tki_clean", "TKI_Group",
                  "blast_tx_raw", "blast_other_tki", "seq_evidence", "seq_inferred",
                  "OS_Status", "Cause_of_Death_clean", "OS_Time"))) %>%
  mutate(Survival = ifelse(OS_Status == 1, "Deceased",
                           ifelse(OS_Status == 0, "Alive", NA_character_))) %>%
  rename(Age = Age_num, TKI_recorded = tki_clean,
         Blast_Phase_Treatment = blast_tx_raw,
         Other_TKI_in_blast_tx = blast_other_tki,
         Cause_of_Death = Cause_of_Death_clean,
         Evidence = seq_evidence, Inferred_sequence = seq_inferred,
         OS_Time_Months = OS_Time)
save_csv(seq_detail, "O1_08_Sequential_TKI_Use_Detail.csv")
print(seq_detail)

seq_summary <- tibble(
  Evidence_source = c(
    "Different TKI in blast-phase treatment field",
    "Recorded ponatinib (prior TKI exposure implied)",
    "Total with any evidence of sequential/different TKI use"),
  n = c(
    sum(df$seq_evidence == "Different TKI documented in blast-phase treatment", na.rm = TRUE),
    sum(df$seq_evidence == "Implied (ponatinib is almost never first-line)", na.rm = TRUE),
    sum(!is.na(df$seq_evidence), na.rm = TRUE)),
  Deceased = c(
    sum(df$died & df$seq_evidence == "Different TKI documented in blast-phase treatment", na.rm = TRUE),
    sum(df$died & df$seq_evidence == "Implied (ponatinib is almost never first-line)", na.rm = TRUE),
    sum(df$died & !is.na(df$seq_evidence), na.rm = TRUE))
)
save_csv(seq_summary, "O1_09_Sequential_TKI_Use_Summary.csv")
print(seq_summary)

# --- 4.5b Line of therapy: first-line vs subsequent-line (panel comment x) -------------
message("    Panel comment (x): line of therapy (first-line vs subsequent-line)...")

lot_tab <- df_groups %>%
  count(LOT, name = "n", sort = TRUE) %>%
  mutate(Percent_of_classifiable = round(100 * n / nrow(df_groups), 1))

lot_breakdown <- df_groups %>%
  filter(LOT == "Subsequent-line TKI") %>%
  count(TKI_Group, tki_clean, name = "n", sort = TRUE) %>%
  mutate(Percent_of_subsequent_line = round(100 * n / sum(n), 1)) %>%
  rename(recorded_TKI = tki_clean)

save_csv(lot_tab, "O1_14a_LOT_Categorisation.csv")
save_csv(lot_breakdown, "O1_14b_LOT_Subsequent_Line_Breakdown.csv")
save_csv(lot_audit, "O1_14c_LOT_Derivation_Patients.csv")
print(lot_tab)
print(lot_breakdown)
message("    Per-patient LOT derivation (recorded TKI, first TKI, imatinib exposure,")
message("    documented switch, derivation route): O1_14c_LOT_Derivation_Patients.csv")
message("    NOTE: LOT = the first TKI the patient received for CML (First-line) vs any")
message("    TKI given after the first TKI (Subsequent-line). Later-line outcomes must")
message("    NOT be read as inherent drug inferiority (confounding by indication).")

# --- 4.6 Ponatinib (3G) patient detail (panel correction block) ---------------------
pon_detail <- df %>%
  filter(TKI_Group == "3G") %>%
  select(any_of(c("Study_ID", "Age_num", "Phase_raw", "tki_clean",
                  "mod_switch", "mod_interruption", "OS_Status",
                  "Cause_of_Death_clean", "OS_Time"))) %>%
  mutate(Survival = ifelse(OS_Status == 1, "Deceased",
                           ifelse(OS_Status == 0, "Alive", NA_character_))) %>%
  rename(Age = Age_num, TKI = tki_clean,
         Cause_of_Death = Cause_of_Death_clean,
         OS_Time_Months = OS_Time)
save_csv(pon_detail, "O1_10_Ponatinib_Patients.csv")
print(pon_detail)

# --- 4.7 T315I mutation (panel comment i) ---------------------------------------------
message("    Panel comment (i): T315I mutation...")

if (!is.null(colmap$t315i)) {
  t3_status <- df %>%
    mutate(T315I_cat = ifelse(is.na(T315I_cat), "Not documented", as.character(T315I_cat))) %>%
    count(T315I_cat, name = "n", sort = TRUE) %>%
    mutate(Pct_of_cohort = round(100 * n / nrow(df), 1))
  n_tested <- sum(df$T315I_cat %in% c("Yes", "No"), na.rm = TRUE)
  t3_status <- bind_rows(
    t3_status,
    tibble(T315I_cat = "TESTED (Yes or No documented)", n = n_tested,
           Pct_of_cohort = round(100 * n_tested / nrow(df), 1)))
  save_csv(t3_status, "O1_11_T315I_Status.csv")
  print(t3_status)

  t3_pos <- df %>%
    filter(T315I_cat == "Yes") %>%
    select(any_of(c("Study_ID", "Age_num", "tki_clean", "TKI_Group",
                    "Phase_raw", "OS_Status", "Cause_of_Death_clean", "OS_Time"))) %>%
    mutate(Survival = ifelse(OS_Status == 1, "Deceased",
                             ifelse(OS_Status == 0, "Alive", NA_character_))) %>%
    rename(Age = Age_num, TKI = tki_clean,
           Cause_of_Death = Cause_of_Death_clean,
           OS_Time_Months = OS_Time)
  save_csv(t3_pos, "O1_12_T315I_Patients.csv")
  print(t3_pos)

  t3_exc <- df %>%
    filter(T315I_cat == "Yes", tki_clean != "ponatinib") %>%
    select(any_of(c("Study_ID", "tki_clean", "TKI_Group", "Phase_raw", "OS_Time")))
  if (nrow(t3_exc) > 0) {
    message("    ", nrow(t3_exc),
            " T315I+ patient(s) not on ponatinib (exceptions file saved).")
    save_csv(t3_exc %>% rename(TKI = tki_clean, OS_Time_Months = OS_Time),
             "O1_13_T315I_Clinical_Exceptions.csv")
  } else {
    message("    All T315I+ patients recorded ponatinib.")
  }
} else {
  message("    T315I column not found - panel comment (i) tables skipped.")
}

# =============================================================================
# 5. OBJECTIVE 2 - FIVE-YEAR OVERALL AND PROGRESSION-FREE SURVIVAL
# =============================================================================
message("\n>>> Objective 2: five-year overall and progression-free survival...")

os_summary_tbl <- NULL
pfs_summary_tbl <- NULL
km5_by_lot <- NULL    # PRIMARY 5-yr (or max FU) estimates, first-line vs subsequent-line
lr_by_lot <- NULL     # PRIMARY log-rank p, first-line vs subsequent-line

# --- 5.1 Overall survival ---------------------------------------------------------------
if (have_os && nrow(os_data) >= 5) {
  fit_os <- survfit(Surv(Time, Event) ~ 1, data = os_data)
  os_summary_tbl <- outcome_summary(fit_os, os_data, "Overall survival (OS)")
  print(os_summary_tbl)

  p_os <- ggsurvplot(
    fit_os, data = os_data,
    title = "Overall survival - CML patients on TKIs (KNH cohort)",
    xlab = "Follow-up (months)", ylab = "Probability of survival",
    conf.int = TRUE, pval = FALSE, legend.labs = "All patients",
    risk.table = FALSE, print = FALSE
  )
  save_fig(p_os, "O2_01_KM_OS_Overall.png")
} else {
  message("    Skipped OS analysis (insufficient OS data).")
}

# --- 5.1b PRIMARY: OS by line of therapy (first-line vs subsequent-line) -----------------
# Two-group design (agreed): the ponatinib group (n = 5) is too small for any
# meaningful separate survival analysis, so all 2G/3G TKIs are pooled as the
# "subsequent-line" group. The 1G/2G/3G breakdown in 5.5 is supplementary only.
lot_leg <- c("First-line TKI" = "First-line TKI",
             "Subsequent-line TKI" = "Subsequent-line TKI")
lot_sub <- "First line = first TKI received for CML; subsequent line = any TKI after the first TKI (2nd/3rd line)"
os_lot_d <- if (have_os && "LOT" %in% names(os_data)) os_data %>% filter(!is.na(LOT)) else NULL
if (!is.null(os_lot_d) && nrow(os_lot_d) >= 5 && length(unique(os_lot_d$LOT)) >= 2) {
  fit_os_lot <- survfit(Surv(Time, Event) ~ LOT, data = os_lot_d)
  p_os_lot <- ggsurvplot(
    fit_os_lot, data = os_lot_d,
    pval = TRUE, conf.int = TRUE, risk.table = TRUE, risk.table.height = 0.2,
    legend.labs = lot_leg,
    title = "Overall survival by line of therapy (PRIMARY: first-line vs subsequent-line)",
    subtitle = lot_sub,
    xlab = "Follow-up (months)", ylab = "Probability of survival",
    print = FALSE
  )
  save_fig(p_os_lot, "O2_02_KM_OS_by_LOT.png")
}

# --- 5.2 Progression-free survival ----------------------------------------------------------
if (have_pfs && nrow(pfs_data) >= 5) {
  fit_pfs <- survfit(Surv(Time, Event) ~ 1, data = pfs_data)
  pfs_summary_tbl <- outcome_summary(fit_pfs, pfs_data, "Progression-free survival (PFS)")
  print(pfs_summary_tbl)

  p_pfs <- ggsurvplot(
    fit_pfs, data = pfs_data,
    title = "Progression-free survival - CML patients on TKIs (KNH cohort)",
    xlab = "Follow-up (months)",
    ylab = "Probability of being progression-free",
    conf.int = TRUE, pval = FALSE, legend.labs = "All patients",
    risk.table = FALSE, print = FALSE
  )
  save_fig(p_pfs, "O2_03_KM_PFS_Overall.png")
} else {
  message("    Skipped PFS analysis (insufficient PFS data).")
}

# --- 5.2b PRIMARY: PFS by line of therapy (first-line vs subsequent-line) -----------------
pfs_lot_d <- if (have_pfs && "LOT" %in% names(pfs_data)) pfs_data %>% filter(!is.na(LOT)) else NULL
if (!is.null(pfs_lot_d) && nrow(pfs_lot_d) >= 5 && length(unique(pfs_lot_d$LOT)) >= 2) {
  fit_pfs_lot <- survfit(Surv(Time, Event) ~ LOT, data = pfs_lot_d)
  p_pfs_lot <- ggsurvplot(
    fit_pfs_lot, data = pfs_lot_d,
    pval = TRUE, conf.int = TRUE, risk.table = TRUE, risk.table.height = 0.2,
    legend.labs = lot_leg,
    title = "Progression-free survival by line of therapy (PRIMARY: first-line vs subsequent-line)",
    subtitle = lot_sub,
    xlab = "Follow-up (months)", ylab = "Probability of being progression-free",
    print = FALSE
  )
  save_fig(p_pfs_lot, "O2_04_KM_PFS_by_LOT.png")
}

# --- 5.3 PRIMARY: KM slide tables by line of therapy (panel comments ii & iii) ------------
km_slide_table <- function(d, endpoint_label, group_var) {
  d <- d[!is.na(d[[group_var]]), , drop = FALSE]
  if (nrow(d) < 5 || length(unique(d[[group_var]])) < 2) return(NULL)
  grps <- levels(droplevels(d[[group_var]]))
  tmax <- max(d$Time)
  event_name <- if (endpoint_label == "OS") "Deaths" else "PFS events"

  gvals <- sapply(grps, function(g) {
    dg <- d[d[[group_var]] == g, , drop = FALSE]
    list(n = nrow(dg),
         events = sum(dg$Event),
         med = km_median(survfit(Surv(Time, Event) ~ 1, data = dg)),
         y1 = km_estimate_at(survfit(Surv(Time, Event) ~ 1, data = dg), 12)$estimate,
         y3 = km_estimate_at(survfit(Surv(Time, Event) ~ 1, data = dg), 36)$estimate,
         y5 = if (tmax >= T_FIVE_YEARS)
                km_estimate_at(survfit(Surv(Time, Event) ~ 1, data = dg), 60)$estimate
              else NA_real_)
  }, simplify = FALSE)

  mk <- function(fun) sapply(grps, function(g) fun(gvals[[g]]))
  n_row <- mk(function(x) as.character(x$n))
  death_row <- mk(function(x) sprintf("%d (%.1f%%)", x$events, 100 * x$events / x$n))
  med_row <- mk(function(x) if (is.na(x$med)) "Not reached" else sprintf("%.1f", x$med))
  y1_row <- mk(function(x) if (is.na(x$y1)) "NA" else sprintf("%.1f", 100 * x$y1))
  y3_row <- mk(function(x) if (is.na(x$y3)) "NA" else sprintf("%.1f", 100 * x$y3))
  y5_row <- mk(function(x) if (is.na(x$y5))
                 sprintf("n/a (max FU %.0f mo)", tmax)
               else sprintf("%.1f", 100 * x$y5))

  p <- logrank_p(d, group_var)

  out <- data.frame(
    Parameter = c(paste0("N (", endpoint_label, ")"),
                  paste0(event_name, ", n (%)"),
                  paste0("Median ", endpoint_label, ", months"),
                  paste0("1-year ", endpoint_label, ", %"),
                  paste0("3-year ", endpoint_label, ", %"),
                  paste0("5-year ", endpoint_label, ", %")),
    check.names = FALSE)
  for (g in grps) out[[g]] <- c(n_row[g], death_row[g], med_row[g],
                                y1_row[g], y3_row[g], y5_row[g])
  out[["p_value"]] <- c("", "", "", "", "", fmt_p(p))
  out
}

if (have_os && nrow(os_data) >= 5) {
  o2_slide_os <- km_slide_table(os_data, "OS", "LOT")
  if (!is.null(o2_slide_os)) {
    save_csv(o2_slide_os, "O2_05_KM_Slide_Table_OS_by_LOT.csv")
    message("    PRIMARY KM slide table by line of therapy (OS):\n")
    print(o2_slide_os)
  }
}
if (have_pfs && nrow(pfs_data) >= 5) {
  o2_slide_pfs <- km_slide_table(pfs_data, "PFS", "LOT")
  if (!is.null(o2_slide_pfs)) {
    save_csv(o2_slide_pfs, "O2_06_KM_Slide_Table_PFS_by_LOT.csv")
    print(o2_slide_pfs)
  }
}

# --- 5.4 PRIMARY: 5-year (or max follow-up) estimates by line of therapy ------------------
km5_lot_rows <- list()
if (have_os && nrow(os_data) >= 5) {
  t <- km_by_group(os_data, "OS", "LOT")
  if (!is.null(t)) km5_lot_rows <- c(km5_lot_rows, list(t))
}
if (have_pfs && nrow(pfs_data) >= 5) {
  t <- km_by_group(pfs_data, "PFS", "LOT")
  if (!is.null(t)) km5_lot_rows <- c(km5_lot_rows, list(t))
}
if (length(km5_lot_rows) > 0) {
  km5_by_lot <- bind_rows(km5_lot_rows)
  save_csv(km5_by_lot, "O2_07_KM_5Year_by_LOT.csv")
  print(km5_by_lot)
}

# --- 5.4b PRIMARY: log-rank p by line of therapy (panel comment x) -------------------------
lot_rows <- list()
if (!is.null(os_lot_d) && nrow(os_lot_d) >= 5 && length(unique(os_lot_d$LOT)) >= 2) {
  p_lot <- logrank_p(os_lot_d, "LOT")
  lot_rows <- c(lot_rows, list(tibble(
    Outcome = "OS", Comparison = "First-line TKI vs Subsequent-line TKI",
    N = nrow(os_lot_d), P_value = p_lot, P_report = fmt_p(p_lot))))
}
if (!is.null(pfs_lot_d) && nrow(pfs_lot_d) >= 5 && length(unique(pfs_lot_d$LOT)) >= 2) {
  p_lot <- logrank_p(pfs_lot_d, "LOT")
  lot_rows <- c(lot_rows, list(tibble(
    Outcome = "PFS", Comparison = "First-line TKI vs Subsequent-line TKI",
    N = nrow(pfs_lot_d), P_value = p_lot, P_report = fmt_p(p_lot))))
}
if (length(lot_rows) > 0) {
  lr_by_lot <- bind_rows(lot_rows)
  save_csv(lr_by_lot, "O2_08_LogRank_by_LOT.csv")
  message("    PRIMARY log-rank p by line of therapy:\n")
  print(lr_by_lot)
}

# --- 5.4c Consolidated overall summary ------------------------------------------------------
o2_parts <- Filter(Negate(is.null), list(os_summary_tbl, pfs_summary_tbl))
if (length(o2_parts) > 0) save_csv(bind_rows(o2_parts), "O2_09_5Year_OS_PFS_Summary.csv")

# --- 5.5 SUPPLEMENTARY: 1G/2G/3G breakdown (descriptive ONLY) ---------------------------------
# Kept for the treatment-patterns description (Objective 1 / panel comment ii),
# NOT as a survival comparison: the ponatinib group has only n = 5 and the
# per-drug 2G numbers are also small, so the 3-group survival estimates are
# far too unstable for any meaningful interpretation.
message("    Supplementary 3-group (1G/2G/3G) breakdown - descriptive only...")
os_grp_d <- if (have_os && nrow(os_data) >= 5) os_data %>% filter(!is.na(TKI_Group)) else NULL
if (!is.null(os_grp_d) && nrow(os_grp_d) >= 5 && length(unique(os_grp_d$TKI_Group)) >= 2) {
  fit_os_grp <- survfit(Surv(Time, Event) ~ TKI_Group, data = os_grp_d)
  present <- levels(droplevels(os_grp_d$TKI_Group))
  grp_labs <- c("1G" = "1G (Imatinib)",
                "2G" = "2G (Dasatinib/Nilotinib/Bosutinib)",
                "3G" = "3G (Ponatinib)")
  p_os_grp <- ggsurvplot(
    fit_os_grp, data = os_grp_d,
    pval = TRUE, conf.int = TRUE, risk.table = TRUE, risk.table.height = 0.2,
    legend.labs = grp_labs[present],
    title = "Overall survival by recorded TKI group (SUPPLEMENTARY, descriptive)",
    subtitle = "Recorded TKI group, not a verified treatment line; ponatinib n = 5",
    xlab = "Follow-up (months)", ylab = "Probability of survival",
    print = FALSE
  )
  save_fig(p_os_grp, "O2_10_KM_OS_by_TKI_Group.png")
}
pfs_grp_d <- if (have_pfs && nrow(pfs_data) >= 5) pfs_data %>% filter(!is.na(TKI_Group)) else NULL
if (!is.null(pfs_grp_d) && nrow(pfs_grp_d) >= 5 && length(unique(pfs_grp_d$TKI_Group)) >= 2) {
  fit_pfs_grp <- survfit(Surv(Time, Event) ~ TKI_Group, data = pfs_grp_d)
  present <- levels(droplevels(pfs_grp_d$TKI_Group))
  grp_labs <- c("1G" = "1G (Imatinib)",
                "2G" = "2G (Dasatinib/Nilotinib/Bosutinib)",
                "3G" = "3G (Ponatinib)")
  p_pfs_grp <- ggsurvplot(
    fit_pfs_grp, data = pfs_grp_d,
    pval = TRUE, conf.int = TRUE, risk.table = TRUE, risk.table.height = 0.2,
    legend.labs = grp_labs[present],
    title = "Progression-free survival by recorded TKI group (SUPPLEMENTARY, descriptive)",
    subtitle = "Recorded TKI group, not a verified treatment line; ponatinib n = 5",
    xlab = "Follow-up (months)", ylab = "Probability of being progression-free",
    print = FALSE
  )
  save_fig(p_pfs_grp, "O2_11_KM_PFS_by_TKI_Group.png")
}
if (have_os && nrow(os_data) >= 5) {
  o2_slide_os_grp <- km_slide_table(os_data, "OS", "TKI_Group")
  if (!is.null(o2_slide_os_grp)) {
    save_csv(o2_slide_os_grp, "O2_12_KM_Slide_Table_OS_by_TKI_Group.csv")
    message("    Supplementary KM slide table by TKI group (OS):\n")
    print(o2_slide_os_grp)
  }
}
if (have_pfs && nrow(pfs_data) >= 5) {
  o2_slide_pfs_grp <- km_slide_table(pfs_data, "PFS", "TKI_Group")
  if (!is.null(o2_slide_pfs_grp)) {
    save_csv(o2_slide_pfs_grp, "O2_13_KM_Slide_Table_PFS_by_TKI_Group.csv")
    print(o2_slide_pfs_grp)
  }
}
o2_grp_parts <- list()
if (have_os && nrow(os_data) >= 5) {
  t <- km_by_group(os_data, "OS", "TKI_Group")
  if (!is.null(t)) o2_grp_parts <- c(o2_grp_parts, list(t))
}
if (have_pfs && nrow(pfs_data) >= 5) {
  t <- km_by_group(pfs_data, "PFS", "TKI_Group")
  if (!is.null(t)) o2_grp_parts <- c(o2_grp_parts, list(t))
}
if (length(o2_grp_parts) > 0) save_csv(bind_rows(o2_grp_parts), "O2_14_KM_5Year_by_TKI_Group.csv")

lr_rows <- list()
if (have_os && nrow(os_data) >= 5) {
  lr_rows <- c(lr_rows, list(tibble(
    Outcome = "OS", Test = "Log-rank across recorded TKI groups (supplementary)",
    P_value = logrank_p(os_data, "TKI_Group"),
    P_report = fmt_p(logrank_p(os_data, "TKI_Group")))))
}
if (have_pfs && nrow(pfs_data) >= 5) {
  lr_rows <- c(lr_rows, list(tibble(
    Outcome = "PFS", Test = "Log-rank across recorded TKI groups (supplementary)",
    P_value = logrank_p(pfs_data, "TKI_Group"),
    P_report = fmt_p(logrank_p(pfs_data, "TKI_Group")))))
}
if (length(lr_rows) > 0) save_csv(bind_rows(lr_rows), "O2_15_LogRank_by_TKI_Group.csv")

# =============================================================================
# 6. OBJECTIVE 3 - PREDICTORS OF SURVIVAL OUTCOMES
# =============================================================================
# (panel comments v & vi: HR with 95% CI + p; full variable table; EPV control)
message("\n>>> Objective 3: predictors of survival (Cox models)...")

# --- Table 1: baseline characteristics by line of therapy (PRIMARY 2-group) -----------
# Stratified by LOT (1L imatinib vs subsequent line) so it matches the primary
# survival comparison and the Cox exposure. The 3-group TKI composition is
# kept as a descriptive row (the ponatinib group has n = 5).
d_grp <- df_groups
grp_present <- if (nrow(d_grp) > 0 && "LOT" %in% names(d_grp))
  levels(droplevels(d_grp$LOT)) else character(0)

num_specs <- list()
if ("Age_num" %in% names(d_grp)) num_specs[["Age (years)"]] <- "Age_num"

cat_specs <- list()
flag_tab <- c("T315I mutation (Yes)" = "flag_T315I",
              "Any ADR" = "flag_ADR_any",
              "Grade 3+ ADR" = "flag_ADR_sev",
              "Disease relapse" = "flag_Relapse",
              "Disease progression" = "flag_Progression",
              "Any comorbidity" = "flag_Comorbidity",
              "Documented TKI dose switch" = "mod_switch")
for (nm in names(flag_tab)) {
  if (flag_tab[[nm]] %in% names(d_grp)) {
    cat_specs[[nm]] <- list(src = flag_tab[[nm]], kind = "flag")
  }
}
if ("Sex" %in% names(d_grp)) cat_specs[["Sex"]] <- list(src = "Sex", kind = "factor")
if ("Phase_group" %in% names(d_grp)) {
  cat_specs[["CML phase at diagnosis"]] <- list(src = "Phase_group", kind = "factor")
}
if ("Age_group" %in% names(d_grp)) {
  cat_specs[["Age group"]] <- list(src = "Age_group", kind = "factor")
}
if ("BCR_monitored" %in% names(d_grp)) {
  cat_specs[["BCR-ABL monitoring done"]] <- list(src = "BCR_monitored", kind = "factor")
}
# NOTE: no "line of therapy" row - the table IS stratified by LOT, so it
# would be constant within columns.

fmt_group <- function(x, kind) {
  if (kind == "flag") fmt_pct(x) else fmt_factor_dist(x)
}

t1_rows <- list()
for (nm in names(num_specs)) {
  src <- num_specs[[nm]]
  vals <- setNames(lapply(c(grp_present, "Overall"), function(g) {
    if (g == "Overall") fmt_median_iqr(d_grp[[src]]) else
      fmt_median_iqr(d_grp[[src]][droplevels(d_grp$LOT) == g])
  }), c(grp_present, "Overall"))
  p <- pval_numeric(d_grp[[src]], droplevels(d_grp$LOT))
  t1_rows[[nm]] <- as.data.frame(c(vals, list(P_value = round(p, 4), Variable = nm)),
                                 stringsAsFactors = FALSE)
}
for (nm in names(cat_specs)) {
  sp <- cat_specs[[nm]]
  vals <- setNames(lapply(c(grp_present, "Overall"), function(g) {
    if (g == "Overall") fmt_group(d_grp[[sp$src]], sp$kind) else
      fmt_group(d_grp[[sp$src]][droplevels(d_grp$LOT) == g], sp$kind)
  }), c(grp_present, "Overall"))
  p <- pval_categorical(d_grp[[sp$src]], droplevels(d_grp$LOT))
  t1_rows[[nm]] <- as.data.frame(c(vals, list(P_value = round(p, 4), Variable = nm)),
                                 stringsAsFactors = FALSE)
}
# Descriptive row: recorded TKI composition within each LOT group (no p -
# the 1L column is imatinib by definition of the grouping).
if (nrow(d_grp) > 0 && length(grp_present) >= 2 && "tki_clean" %in% names(d_grp)) {
  vals <- setNames(lapply(grp_present, function(g)
    fmt_factor_dist(d_grp$tki_clean[droplevels(d_grp$LOT) == g])), grp_present)
  tki_row <- as.data.frame(c(vals,
                             list(Overall = fmt_factor_dist(d_grp$tki_clean),
                                  P_value = NA_real_,
                                  Variable = "Recorded TKI within group (descriptive)")),
                           stringsAsFactors = FALSE)
  t1_rows <- c(list(tki_row), t1_rows)
}
if (length(t1_rows) > 0) {
  tbl1 <- bind_rows(t1_rows)
  tbl1 <- tbl1[c("Variable", grp_present, "Overall", "P_value")]
  save_csv(tbl1, "Table1_Baseline_Characteristics.csv")
  message(">>> Table 1 (baseline characteristics by line of therapy: first-line vs subsequent-line):\n")
  print(tbl1)
}

# --- Cox models (OS and PFS) ------------------------------------------------------------
# Panel comment (xi): Grade 3+ ADR status is FORCED into the multivariable
# models (kept even if its univariable p >= 0.20); line of therapy is forced
# in as well (panel comment x).
force_vars <- intersect(c("flag_ADR_sev", "LOT"), names(os_data))
cox_os <- NULL
cox_pfs <- NULL
if (have_os && nrow(os_data) >= 10) {
  cox_os <- run_cox_models(os_data, force_vars = intersect(force_vars, names(os_data)))
}
if (have_pfs && nrow(pfs_data) >= 10) {
  cox_pfs <- run_cox_models(pfs_data, force_vars = intersect(force_vars, names(pfs_data)))
}

o3_univ <- tibble()
if (!is.null(cox_os) && nrow(cox_os$univ) > 0) {
  o3_univ <- bind_rows(o3_univ, cox_os$univ %>% mutate(Outcome = "OS", .before = 1))
}
if (!is.null(cox_pfs) && nrow(cox_pfs$univ) > 0) {
  o3_univ <- bind_rows(o3_univ, cox_pfs$univ %>% mutate(Outcome = "PFS", .before = 1))
}
if (nrow(o3_univ) > 0) save_csv(o3_univ, "O3_01_Univariate_Cox.csv")

# --- Univariable effect of Grade 3+ ADRs on OS and PFS (panel comment xi) --------------
adr_univ_out <- tibble()
for (ep in c("OS", "PFS")) {
  cx <- if (ep == "OS") cox_os else cox_pfs
  if (!is.null(cx) && nrow(cx$univ) > 0) {
    r <- cx$univ %>% filter(Covariate == "flag_ADR_sev")
    if (nrow(r) > 0) {
      adr_univ_out <- bind_rows(adr_univ_out, r %>% mutate(Outcome = ep, .before = 1))
    }
  }
}
if (nrow(adr_univ_out) > 0) {
  save_csv(adr_univ_out, "O3_08_Grade3plus_ADR_OS_PFS_Univariate.csv")
  message(">>> Univariable effect of Grade 3+ ADR on OS and PFS (panel comment xi):\n")
  print(adr_univ_out)
}

if (!is.null(cox_os) && nrow(cox_os$multiv) > 0) {
  save_csv(cox_os$multiv, "O3_02_Multivariate_Cox_OS.csv")
  message(">>> Multivariate Cox model (OS) - every HR with 95% CI and p:\n")
  print(cox_os$multiv)
  message("    Events = ", cox_os$n_events, "; covariates in final model = ",
          cox_os$n_cov_final, "; events per covariate = ",
          if (cox_os$n_cov_final > 0) round(cox_os$n_events / cox_os$n_cov_final, 1)
          else NA)
}
if (!is.null(cox_pfs) && nrow(cox_pfs$multiv) > 0) {
  save_csv(cox_pfs$multiv, "O3_03_Multivariate_Cox_PFS.csv")
}
if (!is.null(cox_os) && nrow(cox_os$zph) > 0) {
  save_csv(cox_os$zph, "O3_04_PH_Assumptions_OS.csv")
  message(">>> Proportional-hazards test (Schoenfeld residuals, OS):\n")
  print(cox_os$zph)
}
if (!is.null(cox_pfs) && nrow(cox_pfs$zph) > 0) {
  save_csv(cox_pfs$zph, "O3_05_PH_Assumptions_PFS.csv")
}

if (!is.null(cox_os) && !is.null(cox_os$cox_fit)) {
  p_zph <- tryCatch(ggcoxzph(cox_os$cox_fit), error = function(e) NULL)
  if (!is.null(p_zph)) {
    save_fig(p_zph + labs(title = "Proportional-hazards check (final OS model)"),
             "O3_06_Cox_zph_OS.png", w = 7, h = 6)
  }
}

# --- Cox variable table with coding and rationale (panel comment vi) --------------------
# NOTE: the TKI exposure in the Cox models is the TWO-GROUP line of therapy
# (LOT). The 3-level recorded TKI group (1G/2G/3G) is NOT used: the ponatinib
# level has only n = 5 (its hazard ratio is not reliably estimable) and it
# would be collinear with LOT anyway (agreed 2-group design).
cox_var_info <- tibble(
  Variable = c("Phase_group", "Age_group", "Sex", "LOT", "flag_T315I",
               "flag_ADR_any", "flag_ADR_sev", "flag_Progression", "flag_Relapse",
               "flag_Comorbidity", "BCR_monitored", "mod_switch"),
  Type = c("Categorical", "Categorical", "Binary", "Binary", "Binary",
           "Binary", "Binary", "Binary", "Binary", "Binary", "Binary", "Binary"),
  Coding = c(
    "Chronic (ref) / Accelerated / Blast",
    "<40 (ref) / 40-59 / >=60 years",
    "Female (ref) / Male",
    "First-line TKI (ref) / Subsequent-line TKI (any TKI given after the first TKI)",
    "No/other (ref) / Yes",
    "No/other (ref) / Yes",
    "No/other (ref) / Yes",
    "No/other (ref) / Yes",
    "No/other (ref) / Yes",
    "No/other (ref) / Yes",
    "No/other (ref) / Yes",
    "None (ref) / documented dose switch"),
  Rationale = c(
    "Strongest known prognostic factor at diagnosis",
    "Established prognostic factor in CML (EUTOS score)",
    "Demographic confounder",
    "Primary exposure: line of therapy (first TKI for CML vs any TKI after it). Two-group because the ponatinib group (n = 5) is too small for a separate level; first TKI = imatinib by study assumption unless documented otherwise (see O1_14c); later lines follow earlier-TKI failure - confounding by indication expected; forced into model (panel x)",
    "T315I confers resistance to all TKIs except ponatinib",
    "Toxicity/intolerance may drive modification and mortality",
    "Severe toxicity may affect survival; forced into model (panel xi)",
    "Disease course marker (progression)",
    "Relapse indicates treatment failure",
    "Comorbidity may affect treatment tolerance and mortality",
    "Proxy for quality of molecular monitoring",
    "Reflects treatment failure or intolerance"
  )
)

if (!is.null(cox_os)) {
  if (nrow(cox_os$univ) > 0) {
    univ_p <- cox_os$univ %>%
      group_by(Covariate) %>%
      summarise(Univariable_p = min(p.value, na.rm = TRUE), .groups = "drop")
  } else {
    univ_p <- tibble(Covariate = character(), Univariable_p = numeric())
  }
  final_covs <- if (nrow(cox_os$multiv) > 0)
    setdiff(unique(cox_os$multiv$Covariate), "Intercept") else character(0)
  cox_var_tbl <- cox_var_info %>%
    mutate(Available = Variable %in% names(os_data),
           Univariable_p = univ_p$Univariable_p[match(Variable, univ_p$Covariate)],
           Entered_at_p_below_0.20 = !is.na(Univariable_p) & Univariable_p < VAR_ENTRY_P,
           In_final_OS_model = Variable %in% final_covs) %>%
    mutate(Univariable_p = round(Univariable_p, 3))
  save_csv(cox_var_tbl, "O3_07_Cox_Variables_Table_OS.csv")
  message(">>> Cox variable table (OS) - panel comment (vi):\n")
  print(cox_var_tbl)
}

# --- Significance verdict: Grade 3+ ADR among the survival predictors ---------------
# Univariable AND multivariable, for OS and PFS - the direct answer to
# "is Grade 3+ ADR a significant predictor?"
adr_sig <- tibble()
for (ep in c("OS", "PFS")) {
  cx <- if (ep == "OS") cox_os else cox_pfs
  if (is.null(cx)) next
  u <- if (nrow(cx$univ) > 0) cx$univ %>% filter(Covariate == "flag_ADR_sev") else tibble()
  m <- if (nrow(cx$multiv) > 0) cx$multiv %>% filter(Covariate == "flag_ADR_sev") else tibble()
  if (nrow(u) > 0) {
    adr_sig <- bind_rows(adr_sig, tibble(
      Outcome = ep, Model = "Univariable",
      HR = round(u$HR[1], 2),
      CI95 = paste0(round(u$conf.low[1], 2), " - ", round(u$conf.high[1], 2)),
      P_value = u$p.value[1],
      P_report = fmt_p(u$p.value[1]),
      Significant_at_5pct = if (u$p.value[1] < 0.05) "YES" else "NO"))
  }
  if (nrow(m) > 0) {
    adr_sig <- bind_rows(adr_sig, tibble(
      Outcome = ep, Model = "Multivariable (adjusted)",
      HR = round(m$HR[1], 2),
      CI95 = paste0(round(m$conf.low[1], 2), " - ", round(m$conf.high[1], 2)),
      P_value = m$p.value[1],
      P_report = fmt_p(m$p.value[1]),
      Significant_at_5pct = if (m$p.value[1] < 0.05) "YES" else "NO"))
  }
}
if (nrow(adr_sig) > 0) {
  save_csv(adr_sig, "O3_09_Grade3plus_ADR_Significance_Summary.csv")
  message("\n>>> GRADE 3+ ADR AS A PREDICTOR - significance summary (panel xi):\n")
  print(adr_sig)
} else {
  message("\n>>> Grade 3+ ADR significance summary not available (insufficient data).")
}

# =============================================================================
# 7. OBJECTIVE 4 - PREVALENCE OF ADVERSE DRUG REACTIONS
# =============================================================================
# ADRs in this dataset are recorded for the CURRENT TKI (cross-sectional).
message("\n>>> Objective 4: ADR prevalence...")

has_adr_any <- !is.null(colmap$adr_any)
has_adr_sev <- !is.null(colmap$adr_sev)

if (!has_adr_any) {
  message("    'Any ADRs' column not found - ADR prevalence cannot be computed.")
} else {
  n_adr_known <- sum(!is.na(df$flag_ADR_any))
  n_adr <- sum(df$flag_ADR_any, na.rm = TRUE)
  n_adr_undoc <- nrow(df) - n_adr_known
  ci_adr <- prop_ci(n_adr, n_adr_known)

  o4_01 <- tibble(
    Measure = c("Any ADR (documented)", "Grade 3+ (severe) ADR (documented)",
                "No ADR documentation"),
    N_documented = c(n_adr_known, n_adr_known, n_adr_undoc),
    n_ADR = c(n_adr,
              if (has_adr_sev) sum(df$flag_ADR_sev, na.rm = TRUE) else NA_integer_,
              NA_integer_),
    Pct_of_documented = c(
      if (n_adr_known > 0) round(100 * n_adr / n_adr_known, 1) else NA_real_,
      if (has_adr_sev && n_adr_known > 0)
        round(100 * sum(df$flag_ADR_sev, na.rm = TRUE) / n_adr_known, 1)
      else NA_real_,
      NA_real_),
    CI95_lower_pct = c(if (is.na(ci_adr[1])) NA_real_ else round(100 * ci_adr[1], 1),
                       NA_real_, NA_real_),
    CI95_upper_pct = c(if (is.na(ci_adr[2])) NA_real_ else round(100 * ci_adr[2], 1),
                       NA_real_, NA_real_),
    Pct_of_cohort = c(if (nrow(df) > 0) round(100 * n_adr / nrow(df), 1) else NA_real_,
                      if (has_adr_sev && nrow(df) > 0)
                        round(100 * sum(df$flag_ADR_sev, na.rm = TRUE) / nrow(df), 1)
                      else NA_real_,
                      round(100 * n_adr_undoc / nrow(df), 1))
  )
  save_csv(o4_01, "O4_01_ADR_Prevalence_Overall.csv")
  print(o4_01)

  # --- 7.1 ADR prevalence by individual TKI -------------------------------------
  o4_02 <- df %>%
    filter(!is.na(tki_clean)) %>%
    group_by(TKI = tki_clean) %>%
    summarise(
      n = n(),
      n_ADR_known = sum(!is.na(flag_ADR_any)),
      n_ADR = sum(flag_ADR_any, na.rm = TRUE),
      Prevalence_ADR_pct = if (n_ADR_known > 0) round(100 * n_ADR / n_ADR_known, 1) else NA_real_,
      .groups = "drop"
    )
  if (has_adr_sev) {
    o4_02_sev <- df %>%
      filter(!is.na(tki_clean)) %>%
      group_by(TKI = tki_clean) %>%
      summarise(
        n_Severe = sum(flag_ADR_sev, na.rm = TRUE),
        Prevalence_Severe_pct = if (sum(!is.na(flag_ADR_any)) > 0)
          round(100 * sum(flag_ADR_sev, na.rm = TRUE) / sum(!is.na(flag_ADR_any)), 1)
        else NA_real_,
        .groups = "drop"
      )
    o4_02 <- o4_02 %>% left_join(o4_02_sev, by = "TKI")
  } else {
    o4_02 <- o4_02 %>% mutate(n_Severe = NA_integer_, Prevalence_Severe_pct = NA_real_)
  }
  o4_02 <- o4_02 %>% arrange(desc(n))
  p_adr_by_tki <- cat_pvalue(table(df$tki_clean, df$flag_ADR_any, useNA = "no"))
  o4_02$p_value_TKI_comparison <- p_adr_by_tki
  o4_02$p_value_report <- fmt_p(p_adr_by_tki)
  save_csv(o4_02, "O4_02_ADR_Prevalence_by_TKI.csv")
  print(o4_02)

  # --- 7.1b ADR prevalence by recorded TKI group ---------------------------------
  o4_03 <- df_groups %>%
    group_by(TKI_Group) %>%
    summarise(
      n = n(),
      n_ADR_known = sum(!is.na(flag_ADR_any)),
      n_ADR = sum(flag_ADR_any, na.rm = TRUE),
      Prevalence_ADR_pct = if (n_ADR_known > 0) round(100 * n_ADR / n_ADR_known, 1) else NA_real_,
      n_Severe = if (has_adr_sev) sum(flag_ADR_sev, na.rm = TRUE) else NA_integer_,
      Prevalence_Severe_pct = if (has_adr_sev && n_ADR_known > 0)
        round(100 * sum(flag_ADR_sev, na.rm = TRUE) / n_ADR_known, 1)
      else NA_real_,
      .groups = "drop"
    )
  save_csv(o4_03, "O4_03_ADR_Prevalence_by_Group.csv")
  print(o4_03)

  # --- 7.2 Grade 3+ ADR types (panel comment iv) ----------------------------------
  if (has_adr_sev && !is.null(colmap$adr_sev_sp)) {
    g3_rows <- df %>%
      filter(flag_ADR_sev == TRUE) %>%
      transmute(TKI = tki_clean, g3txt = as.character(.[[colmap$adr_sev_sp]]))
    g3_long <- g3_rows %>%
      filter(!is.na(g3txt)) %>%
      separate_rows(g3txt, sep = "[,;]") %>%
      mutate(g3txt = trimws(g3txt), g3txt = ifelse(g3txt == "", NA_character_, g3txt)) %>%
      filter(!is.na(g3txt))
    if (nrow(g3_long) > 0) {
      n_g3_patients <- nrow(g3_rows)
      o4_04 <- g3_long %>%
        count(g3txt, sort = TRUE, name = "n") %>%
        mutate(Pct_of_grade3_cases = round(100 * n / nrow(g3_long), 1),
               Pct_of_grade3_patients = round(100 * n / max(1, n_g3_patients), 1))
      o4_04 <- bind_rows(
        tibble(g3txt = "TOTAL specified entries (patients: ",
               n = nrow(g3_long),
               Pct_of_grade3_cases = NA_real_,
               Pct_of_grade3_patients = NA_real_),
        o4_04)
      save_csv(o4_04, "O4_04_Grade3_Specified_ADRs.csv")
      message("    Panel comment (iv): Grade 3+ ADR types (",
              n_g3_patients, " patients with Grade 3+ ADR = Yes):\n")
      print(o4_04)
    } else {
      message("    No Grade 3 ADR specification text found.")
    }
  } else {
    message("    Grade 3 ADR 'specified' column not found - panel comment (iv) type table skipped.")
  }

  # --- 7.3 Drug-specific ADR profiles (panel comment ix) ----------------------------
  adr_source_cols <- intersect(
    c(colmap$adr_name, colmap$toxicity_cols), names(df))
  if (length(adr_source_cols) > 0) {
    adr_long <- df %>%
      select(any_of(c("Study_ID", "tki_clean", adr_source_cols))) %>%
      filter(!is.na(tki_clean)) %>%
      pivot_longer(all_of(adr_source_cols), names_to = "source", values_to = "adr_raw") %>%
      filter(!is.na(adr_raw)) %>%
      mutate(adr_raw = trimws(as.character(adr_raw))) %>%
      filter(nzchar(adr_raw),
             !tolower(adr_raw) %in%
               c("not applicable", "not documented", "none", "0", "na", "-", "no")) %>%
      separate_rows(adr_raw, sep = "[,;]") %>%
      mutate(adr_raw = trimws(adr_raw), adr_type = norm_adr(adr_raw)) %>%
      filter(nzchar(adr_raw), !is.na(adr_type)) %>%
      distinct(Study_ID, tki_clean, adr_type)

    adr_denom <- df %>%
      filter(!is.na(tki_clean)) %>%
      count(tki_clean, name = "n_tki")

    o4_05 <- adr_long %>%
      group_by(tki_clean, adr_type) %>%
      summarise(n = n(), .groups = "drop") %>%
      left_join(adr_denom, by = "tki_clean") %>%
      mutate(pct = round(100 * n / n_tki, 1),
             cell = paste0(n, " (", pct, "%)")) %>%
      pivot_wider(names_from = tki_clean, values_from = cell, values_fill = "0")

    # add class-signal labels and the "any ADR" row
    o4_05$Class_signal <- vapply(o4_05$adr_type, class_signal, character(1))
    o4_05 <- o4_05 %>% select(adr_type, Class_signal, everything())

    any_row <- df %>%
      filter(!is.na(tki_clean)) %>%
      group_by(tki_clean) %>%
      summarise(any_n = sum(flag_ADR_any, na.rm = TRUE),
                any_known = sum(!is.na(flag_ADR_any)),
                .groups = "drop") %>%
      mutate(any_pct = if (any_known > 0) round(100 * any_n / any_known, 1) else NA_real_,
             cell = paste0(any_n, " (", any_pct, "%)")) %>%
      pivot_wider(names_from = tki_clean, values_from = cell, values_fill = "0") %>%
      mutate(adr_type = "ANY ADR REPORTED",
             Class_signal = "(any ADR documented)") %>%
      select(-any_n, -any_known, -any_pct)
    o4_05 <- bind_rows(any_row, o4_05)
    o4_05 <- o4_05 %>%
      arrange(adr_type != "ANY ADR REPORTED") %>%
      select(adr_type, Class_signal, everything())
    save_csv(o4_05, "O4_05_ADR_by_Individual_Drug.csv")
    message("    Panel comment (ix): drug-specific ADR profile (n and % within each TKI):\n")
    print(o4_05)
  } else {
    message("    No individual ADR text columns found (Name_of_ADR / Toxicity_1..4);")
    message("    panel comment (ix) table skipped.")
  }

  # --- 7.4 ADR-related TKI modifications --------------------------------------------
  if (!is.null(colmap$reason)) {
    sw <- df_groups %>% filter(mod_switch)
    o4_06 <- sw %>%
      group_by(TKI_Group) %>%
      summarise(
        n_switch = n(),
        n_adr_reason = sum(reason_adr_related, na.rm = TRUE),
        Pct_adr_reason = if (n() > 0)
          round(100 * sum(reason_adr_related, na.rm = TRUE) / n(), 1)
        else NA_real_,
        .groups = "drop"
      )
    all_row <- sw %>%
      summarise(
        TKI_Group = "All",
        n_switch = n(),
        n_adr_reason = sum(reason_adr_related, na.rm = TRUE),
        Pct_adr_reason = if (n() > 0)
          round(100 * sum(reason_adr_related, na.rm = TRUE) / n(), 1)
        else NA_real_
      )
    o4_06 <- bind_rows(o4_06 %>% mutate(TKI_Group = as.character(TKI_Group)), all_row)
    save_csv(o4_06, "O4_06_ADR_related_TKI_Modifications.csv")
    print(o4_06)
  }

  # --- 7.5 ADR figure -----------------------------------------------------------------
  p_adr <- o4_02 %>%
    mutate(TKI = reorder(factor(TKI), Prevalence_ADR_pct)) %>%
    ggplot(aes(x = TKI, y = Prevalence_ADR_pct)) +
    geom_col(fill = "steelblue", width = 0.65) +
    geom_text(aes(label = paste0(n_ADR, " / ", n)), vjust = -0.7, size = 3.5) +
    labs(title = "Prevalence of any ADR by current TKI",
         subtitle = sprintf("p across TKIs = %s", fmt_p(p_adr_by_tki)),
         x = "Current TKI", y = "% patients with any ADR") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
  save_fig(p_adr, "O4_07_ADR_Prevalence_by_TKI.png", w = 9, h = 5.5)
}

# =============================================================================
# 8. DRAFT RESULTS NARRATIVE (PANEL-RESPONSE VERSION)
# =============================================================================
cat("\n")
cat("=====================================================================\n")
cat(" DRAFT RESULTS NARRATIVE - PANEL RESPONSE VERSION\n")
cat("=====================================================================\n\n")

# --- Objective 1 (treatment patterns) --------------------------------------------
cat("OBJECTIVE 1 - TKI TREATMENT PATTERNS\n\n")

if (nrow(o1_01) > 0) {
  cat("Recorded TKI groups: ")
  parts <- sprintf("%s (%s) n = %d (%.1f%%)",
                   o1_01$TKI_Group, o1_01$Drugs, o1_01$n, o1_01$Pct_of_cohort)
  cat(paste(parts, collapse = "; "), ".\n\n", sep = "")
}

n_seq <- sum(!is.na(df$seq_evidence), na.rm = TRUE)
n_seq_dead <- sum(df$died & !is.na(df$seq_evidence), na.rm = TRUE)
cat("Sequential TKI use (panel comment viii):\n")
cat("Complete dated sequential treatment histories were not captured. Indirect\n")
cat("evidence of sequential/different TKI exposure was identified in ",
    n_seq, " patients (different TKI documented in the\n", sep = "")
cat("blast-phase treatment field; ponatinib recipients by clinical convention).\n")
cat("Of these, ", n_seq_dead, " were deceased. This is a marker of treatment-resistant,\n", sep = "")
cat("advanced-phase disease, not evidence that switching itself worsened outcomes.\n\n")

if (!is.null(colmap$t315i)) {
  cat("T315I mutation (panel comment i):\n")
  cat("Mutation testing was documented in ",
      sum(df$T315I_cat %in% c("Yes", "No"), na.rm = TRUE), " of ", nrow(df),
      " patients (", round(100 * sum(df$T315I_cat %in% c('Yes','No'), na.rm = TRUE) / nrow(df), 1),
      "%). Among those tested, T315I was detected in ",
      sum(df$T315I_cat == 'Yes', na.rm = TRUE), " patients. See O1_12_T315I_Patients.csv for the\n", sep = "")
  cat("per-patient detail (TKI used, phase, cause of death). Limitation: T315I was\n")
  cat("the only mutation captured; other BCR-ABL1 mutations could not be assessed.\n\n")
}

if ("LOT" %in% names(df)) {
  n_first <- sum(df$LOT == "First-line TKI", na.rm = TRUE)
  n_subseq <- sum(df$LOT == "Subsequent-line TKI", na.rm = TRUE)
  n_lot <- max(1, n_first + n_subseq)
  cat("Line of therapy (panel comment x):\n")
  cat(paste0("Patients were categorised by line of therapy: ", n_first,
             " (", round(100 * n_first / n_lot, 1),
             "%) received the FIRST-LINE TKI (the first TKI they received for CML) and ",
             n_subseq, " (", round(100 * n_subseq / n_lot, 1),
             "%) were on a SUBSEQUENT-LINE TKI (any TKI given after the first TKI,\n",
             "e.g. 2nd-line or 3rd-line). Categorisation: O1_14a_LOT_Categorisation.csv;\n"))
  cat("subsequent-line breakdown: O1_14b_LOT_Subsequent_Line_Breakdown.csv;\n")
  cat("per-patient derivation: O1_14c_LOT_Derivation_Patients.csv; log-rank p-values:\n")
  cat("O2_08_LogRank_by_LOT.csv.\n")
  cat("With this dataset the first TKI is imatinib by study assumption (unless a\n")
  cat("patient is recorded as never treated with imatinib), so imatinib = first line\n")
  cat("and any 2G/3G TKI = subsequent line with prior imatinib. LOT is derived from\n")
  cat("the first TKI received - it is NOT a verified dated treatment sequence, and\n")
  cat("later-line outcomes are confounded by indication.\n")
  cat("LOT (first-line vs subsequent-line) is the PRIMARY two-group survival\n")
  cat("comparison: the ponatinib group (n = 5) is too small for any meaningful\n")
  cat("separate survival analysis, so all 2G/3G TKIs are pooled as 'subsequent line'.\n\n")
}

# --- Objective 2 (survival) -----------------------------------------------------------
cat("OBJECTIVE 2 - FIVE-YEAR OVERALL AND PROGRESSION-FREE SURVIVAL\n\n")

if (!is.null(os_summary_tbl)) {
  r <- os_summary_tbl
  cat("Among ", r$N[1], " patients with documented follow-up, ", r$Events[1],
      " (", r$Event_rate_pct[1], "%) died (median time to death/censoring: ",
      r$Median_time_months[1], " months). ",
      sep = "")
  if (r$Five_year_estimable[1]) {
    cat("The 5-year (60-month) overall survival was ", r$Estimate_pct[1],
        "% (95% CI ", r$CI95_lower_pct[1], " - ", r$CI95_upper_pct[1],
        "%); the restricted mean survival time over 60 months was ",
        r$RMST_months[1], " months.\n\n", sep = "")
  } else {
    cat("Observed follow-up did not reach 60 months (maximum ",
        r$Max_followup_months[1],
        " months), so a 5-year OS estimate is not available; the OS estimate at the ",
        "longest follow-up (", r$Estimate_at_months[1], " months) was ",
        r$Estimate_pct[1], "% (95% CI ", r$CI95_lower_pct[1], " - ",
        r$CI95_upper_pct[1], "%).\n\n", sep = "")
  }
} else {
  cat("OS analysis could not be performed (time-to-event data not found).\n\n")
}

if (!is.null(pfs_summary_tbl)) {
  r <- pfs_summary_tbl
  cat("PFS: ", r$Events[1], " of ", r$N[1], " patients (", r$Event_rate_pct[1],
      "%) experienced a PFS event (median time ", r$Median_time_months[1],
      " months). ", sep = "")
  if (r$Five_year_estimable[1]) {
    cat("The 5-year (60-month) progression-free survival was ", r$Estimate_pct[1],
        "% (95% CI ", r$CI95_lower_pct[1], " - ", r$CI95_upper_pct[1], "%).\n\n",
        sep = "")
  } else {
    cat("Observed follow-up did not reach 60 months (maximum ",
        r$Max_followup_months[1],
        " months); the PFS estimate at the longest follow-up (",
        r$Estimate_at_months[1], " months) was ", r$Estimate_pct[1],
        "% (95% CI ", r$CI95_lower_pct[1], " - ", r$CI95_upper_pct[1], "%).\n\n",
        sep = "")
  }
}

cat("PRIMARY COMPARISON - first-line (imatinib) vs subsequent line (2G/3G TKI):\n")
if (!is.null(km5_by_lot) && nrow(km5_by_lot) > 0) {
  for (ep in c("OS", "PFS")) {
    r1 <- km5_by_lot[km5_by_lot$Outcome == ep & km5_by_lot$Group == "First-line TKI", , drop = FALSE]
    r2 <- km5_by_lot[km5_by_lot$Outcome == ep & km5_by_lot$Group == "Subsequent-line TKI", , drop = FALSE]
    if (nrow(r1) != 1 || nrow(r2) != 1) next
    p_rep <- "NA"
    if (!is.null(lr_by_lot)) {
      sel_p <- lr_by_lot$P_report[lr_by_lot$Outcome == ep]
      if (length(sel_p) >= 1) p_rep <- sel_p[1]
    }
    cat("  ", ep, " at ", r1$Estimate_at_months[1],
        " months (60 where follow-up allowed):\n", sep = "")
    cat("     First-line TKI: n = ", r1$N[1], " (", r1$Events[1],
        " events); ", r1$Estimate_pct[1], "% (95% CI ",
        r1$CI95_lower_pct[1], " - ", r1$CI95_upper_pct[1], ")\n", sep = "")
    cat("     Subsequent-line TKI: n = ", r2$N[1], " (", r2$Events[1],
        " events); ", r2$Estimate_pct[1], "% (95% CI ",
        r2$CI95_lower_pct[1], " - ", r2$CI95_upper_pct[1], ")\n", sep = "")
    cat("     Log-rank p = ", p_rep, "\n\n", sep = "")
  }
}
cat("The subsequent-line group is selected for prior treatment failure, relapse\n")
cat("or progression (confounding by indication) - any survival difference\n")
cat("reflects selection, not intrinsic drug inferiority.\n\n")
cat("Supplementary 3-group breakdown (descriptive only - the ponatinib group has\n")
cat("n = 5 and per-drug 2G numbers are also small, so it is too unstable for\n")
cat("meaningful separate survival analysis): O2_10 to O2_15. The\n")
cat("generation-outcome discrepancy (panel comment ii) is explained by\n")
cat("confounding by indication plus advanced-phase disease at presentation\n")
cat("(see O1_02_Phase_by_TKI_Group.csv) and shorter follow-up in later groups.\n\n")

# --- Objective 3 (predictors) -----------------------------------------------------------
cat("OBJECTIVE 3 - PREDICTORS OF SURVIVAL (panel comments v & vi)\n\n")
for (ep in c("OS", "PFS")) {
  cx <- if (ep == "OS") cox_os else cox_pfs
  if (is.null(cx) || is.null(cx$multiv) || nrow(cx$multiv) == 0) next
  m <- cx$multiv %>% filter(Covariate != "Intercept")
  if (nrow(m) == 0) next
  cat("Multivariate Cox model for ", ep, " (", cx$n_used,
      " patients with complete covariate data; entry at univariate p < ",
      VAR_ENTRY_P, "; ", cx$n_events, " events / ", cx$n_cov_final,
      " covariates = ",
      if (cx$n_cov_final > 0) round(cx$n_events / cx$n_cov_final, 1) else NA,
      " events per covariate): \n", sep = "")
  for (i in seq_len(nrow(m))) {
    lab <- if (m$Level[i] == "") m$Label[i] else paste0(m$Label[i], " (", m$Level[i], ")")
    cat("   - ", lab, ": ", m$Report[i], "\n", sep = "")
  }
  cat("\n")
}
if (is.null(cox_os) && is.null(cox_pfs)) {
  cat("Cox models could not be performed (time-to-event data not found).\n\n")
}

if (nrow(adr_univ_out) > 0) {
  cat("Univariable effect of Grade 3+ ADR on survival (panel comment xi):\n")
  for (i in seq_len(nrow(adr_univ_out))) {
    cat("   - ", adr_univ_out$Outcome[i], ": ", adr_univ_out$Report[i], "\n", sep = "")
  }
  cat("Grade 3+ ADR status (yes/no) was also FORCED into the multivariable OS and PFS\n")
  cat("models above (retained even if its univariable p >= 0.20); line of therapy was\n")
  cat("forced in as well (panel comment x).\n\n")
}

if (nrow(adr_sig) > 0) {
  cat("IS Grade 3+ ADR A SIGNIFICANT PREDICTOR? (univariable vs multivariable)\n")
  for (i in seq_len(nrow(adr_sig))) {
    cat("   - ", adr_sig$Outcome[i], " / ", adr_sig$Model[i], ": HR ",
        adr_sig$HR[i], " (95% CI ", adr_sig$CI95[i], "), p = ",
        adr_sig$P_report[i], " -> ", adr_sig$Significant_at_5pct[i],
        " significant at the 5% level\n", sep = "")
  }
  cat("Full table: O3_09_Grade3plus_ADR_Significance_Summary.csv\n\n")
}

# --- Objective 4 (ADRs) -------------------------------------------------------------------
cat("OBJECTIVE 4 - ADR PREVALENCE (panel comments iv & ix)\n\n")
if (!has_adr_any) {
  cat("ADR analysis could not be performed (no 'Any ADRs' column found).\n\n")
} else {
  cat("Any ADR: ", o4_01$n_ADR[1], " of ", o4_01$N_documented[1],
      " documented patients (", o4_01$Pct_of_documented[1],
      "%; 95% CI ", o4_01$CI95_lower_pct[1], " - ", o4_01$CI95_upper_pct[1],
      "%). Grade 3+ ADR: ",
      o4_01$n_ADR[2], " (", o4_01$Pct_of_documented[2],
      "%). No ADR documentation in ",
      o4_01$N_documented[3], " patients (", o4_01$Pct_of_cohort[3],
      "%) - severe-ADR prevalence is likely underestimated.\n", sep = "")
  if (!is.na(p_adr_by_tki)) {
    if (p_adr_by_tki < 0.05) {
      cat("ADR prevalence differed significantly across TKIs (p = ",
          fmt_p(p_adr_by_tki), ").\n", sep = "")
    } else {
      cat("No significant difference in ADR prevalence across TKIs (p = ",
          fmt_p(p_adr_by_tki), ").\n", sep = "")
    }
  }
  cat("Drug-specific profile (n, % within each TKI): O4_05_ADR_by_Individual_Drug.csv.\n\n")
}

# --- Caveats ---------------------------------------------------------------------------------
cat("---------------------------------------------------------------------\n")
cat("CAVEATS FOR THE DISCUSSION\n")
cat("---------------------------------------------------------------------\n")
cat("- The 1G/2G/3G labels are RECORDED TKI GROUPS, not verified treatment\n")
cat("  lines; they must not be presented as a confirmed treatment flow.\n")
cat("- Line of therapy (LOT) was derived from the first TKI the patient\n")
cat("  received for CML: First-line TKI = the first TKI; Subsequent-line TKI =\n")
cat("  any TKI given after the first TKI (2nd/3rd line). With this dataset the\n")
cat("  first TKI is imatinib by study assumption (unless a patient is recorded\n")
cat("  as never treated with imatinib), so imatinib = first line and any 2G/3G\n")
cat("  TKI = subsequent line with prior imatinib presumed. Per-patient\n")
cat("  derivation: O1_14c_LOT_Derivation_Patients.csv. Worse outcomes on later\n")
cat("  lines reflect confounding by indication, not drug inferiority.\n")
cat("- The PRIMARY survival comparison is the two-group line of therapy\n")
cat("  (first-line imatinib vs subsequent line, all 2G/3G TKIs pooled); 2G\n")
cat("  TKIs are not analysed separately either (numbers too small). The\n")
cat("  ponatinib group (n = 5) is too small for any meaningful separate\n")
cat("  survival analysis; the 1G/2G/3G KM breakdown (O2_10..O2_15) is\n")
cat("  supplementary/descriptive only and must not be used for inference.\n")
cat("- A blank modification field means no modification was recorded, not\n")
cat("  that the patient never switched; verify the 'Dose switch' code against\n")
cat("  the study codebook before interpreting it as a TKI change.\n")
cat("- Sequential TKI use is inferred indirectly; no ordered, dated\n")
cat("  three-drug sequence can be established. Mortality alone does not\n")
cat("  establish whether switching was or was not a rescue.\n")
cat("- Only T315I was captured; other BCR-ABL1 mutations were not assessed.\n")
cat("- Patients on 2G/3G TKIs likely represent treatment resistance,\n")
cat("  intolerance or advanced disease (confounding by indication); poorer\n")
cat("  outcomes on later groups should not be read as inferior drug efficacy.\n")
cat("- ADRs are recorded for the current TKI (cross-sectional) and grading is\n")
cat("  undocumented in a substantial fraction of the cohort.\n")
est5 <- c()
if (!is.null(os_summary_tbl))  est5 <- c(est5, os_summary_tbl$Five_year_estimable[1])
if (!is.null(pfs_summary_tbl)) est5 <- c(est5, pfs_summary_tbl$Five_year_estimable[1])
if (length(est5) > 0 && !all(est5)) {
  cat("- Some 5-year estimates are limited by follow-up shorter than 60 months.\n")
}

# =============================================================================
# 9. OUTPUT INDEX
# =============================================================================
cat("\n---------------------------------------------------------------------\n")
cat("OUTPUT FILES (saved in: ", out_dir, ")\n", sep = "")
cat("---------------------------------------------------------------------\n")
for (f in saved_files) cat("  ", basename(f), "\n")
cat("\n>>> Analysis complete.\n")
