# =============================================================================
# CML THESIS - OBJECTIVE 2 (survival) ANALYSIS  -  LOT + GRADE >=3 ADR
# "To investigate predictors of survival outcomes among chronic myeloid
#  leukaemia patients managed with TKIs at Kenyatta National Hospital."
#
# Standalone R script. Reads the Excel workbook (trial 5.xlsx), then produces
# EVERY table and graph for the Objective-2 analysis:
#
#   1. Line of Therapy (LOT):  First-line = Imatinib
#                              Subsequent-line = Dasatinib/Nilotinib/
#                                                Bosutinib/Ponatinib
#   2. Grade >=3 ADR as a Yes/No variable
#   3. OS event-coding check
#   4. PFS event-coding check (death MUST be a PFS event)
#   5. LOT descriptive tables
#   6. Grade >=3 ADR descriptive tables
#   7. KM curves: OS by LOT; PFS by LOT; OS by Grade>=3 ADR; PFS by Grade>=3 ADR
#   8. Univariable Cox: LOT->OS, LOT->PFS, Grade>=3 ADR->OS, Grade>=3 ADR->PFS
#   9. Multivariable Cox (OS and PFS) with the FIXED adjusted set:
#         CML phase + gender + age + LOT + Grade>=3 ADR   (no stepwise)
#  10. Proportional-hazards (Schoenfeld) tests for all of the above models
#
# HOW TO RUN
#   1. RStudio -> open this file -> Source.
#   2. Pick the Excel file in the window (or type its path).
#   3. Everything is written to a results folder next to the Excel file:
#        * KM_*.png graphs
#        * Objective2_Analysis_Tables.xlsx (all tables in one workbook)
#        * the same tables as individual .csv files
#        * model + PH summaries as .txt
#
# CODING DECISIONS (must match the rest of the thesis analysis):
#   * Grade >=3 ADR: "Yes" -> Yes; "No", "Not documented", "Not applicable"
#     -> No.  (Assumption, stated in the thesis: 31 of the 47 "Not documented"
#     patients had no ADR recorded at all; grading was simply absent.)
#   * Inferential analyses (KM comparisons, Cox) use the LOT-classifiable
#     cohort (n = 162 in trial 5; 1 patient with a blank TKI field excluded).
#   * PFS event: death is an event. The script VERIFYS this and stops with a
#     message if any death is coded PFS_event = 0.
#
# CROSS-CHECK VALUES (what trial 5.xlsx must produce - verified independently
# with the Python companion analysis, analysis_objective3.py):
#   LOT: first-line 97 (59.9%), subsequent-line 65 (40.1%)
#   Grade>=3 ADR: Yes 37 (22.8%), No 125
#   OS events 70; PFS events 92 (70 deaths + 22 non-death); 0 deaths missed
#   Univariable: LOT->OS  HR 2.17 (1.35-3.50) p=0.001 ; LOT->PFS HR 1.54
#                (1.02-2.33) p=0.040 ; G3ADR->OS HR 0.91 (0.50-1.64) p=0.75 ;
#                G3ADR->PFS HR 1.09 (0.67-1.79) p=0.72
#   Multivariable OS: Accelerated HR ~8.0, Blast HR ~12.4 (p<0.001),
#                LOT HR ~1.56 (p~0.09), G3ADR HR ~1.11 (ns)
#   Log-rank: OS by LOT p ~0.001 ; PFS by LOT p ~0.04 ; by G3ADR p ~0.75/0.72
# =============================================================================

required <- c("readxl", "dplyr", "tidyr", "stringr", "survival", "survminer",
              "ggplot2", "broom", "openxlsx", "gridExtra")
missing <- required[!sapply(required, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) install.packages(missing, dependencies = TRUE)
invisible(lapply(required, library, character.only = TRUE))

options(width = 110)

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

# ---------------------------------------------------------------- 1. FILE ----
excel_file <- tryCatch(file.choose(), error = function(e) NA_character_)
if (is.na(excel_file) || !nzchar(excel_file) || !file.exists(excel_file)) {
  typed <- tryCatch({
    cat("Type the FULL path to the .xlsx file and press Enter:\n")
    gsub('"', "", trimws(readLines(con = stdin(), n = 1)))
  }, error = function(e) "")
  if (nzchar(typed) && file.exists(typed)) excel_file <- typed
}
if (is.na(excel_file) && file.exists("trial 5.xlsx")) excel_file <- "trial 5.xlsx"
if (is.na(excel_file) || !file.exists(excel_file))
  stop("No Excel file selected/found.", call. = FALSE)

sheets <- readxl::excel_sheets(excel_file)
sheet_to_use <- if ("Sheet2" %in% sheets) "Sheet2" else sheets[1]
dat <- readxl::read_excel(excel_file, sheet = sheet_to_use)
cat("\n>>> File:", basename(excel_file), "| sheet:", sheet_to_use,
    "|", nrow(dat), "rows x", ncol(dat), "cols\n")

# ------------------------------------------- 2. STANDARDISE COLUMN NAMES ----
clean_names <- function(x) {
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  stringr::str_replace_all(x, "^_|_$", "")
}
names(dat) <- clean_names(names(dat))

find_col <- function(candidates, required = TRUE) {
  n <- names(dat)
  hit <- n[tolower(n) %in% tolower(candidates)]
  if (length(hit)) return(hit[1])
  for (z in candidates) {
    hit <- n[stringr::str_detect(tolower(n), stringr::fixed(tolower(z)))]
    if (length(hit)) return(hit[1])
  }
  if (required)
    stop("Required column not found. Tried: ", paste(candidates, collapse = ", "))
  NA_character_
}

# -------------------------------------------- 3. IDENTIFY VARIABLES ----
tki_col    <- find_col(c("TKI_Used", "TKI"))
phase_col  <- find_col(c("Phase_of_CML", "Phase", "CML_Phase"))
gender_col <- find_col(c("Gender", "Sex"))
age_col    <- find_col(c("Age"))
os_time_col  <- find_col(c("OS_Time_Months", "OS_time", "OS_Time"), FALSE)
os_event_col <- find_col(c("OS_Event", "OS_event"), FALSE)
pfs_time_col  <- find_col(c("PFS_Time_Months", "PFS_time", "PFS_Time"), FALSE)
pfs_event_col <- find_col(c("PFS_Event", "PFS_event"), FALSE)
adr_col    <- find_col(c("Grade_3_ADR", "Grade_3plus_ADR"), FALSE)

# -------------------------------------- 4. TKI CLEANING AND LOT (spec) ----
dat <- dat %>%
  mutate(
    TKI_clean = stringr::str_to_lower(stringr::str_trim(as.character(.data[[tki_col]]))),
    Phase_clean = stringr::str_to_lower(stringr::str_trim(as.character(.data[[phase_col]]))),
    Gender_clean = as.character(.data[[gender_col]]),
    Age_num = suppressWarnings(as.numeric(.data[[age_col]]))
  ) %>%
  mutate(
    TKI_clean = case_when(
      str_detect(TKI_clean, "imatin") ~ "Imatinib",
      str_detect(TKI_clean, "dasatin") ~ "Dasatinib",
      str_detect(TKI_clean, "nilotin") ~ "Nilotinib",
      str_detect(TKI_clean, "bosutin") ~ "Bosutinib",
      str_detect(TKI_clean, "ponatin") ~ "Ponatinib",
      TRUE ~ NA_character_)
  )

# Supervisor-requested LOT grouping
dat <- dat %>%
  mutate(
    LOT = factor(
      case_when(
        TKI_clean == "Imatinib" ~ "First-line TKI",
        TKI_clean %in% c("Dasatinib", "Nilotinib", "Bosutinib", "Ponatinib") ~
          "Subsequent-line TKI",
        TRUE ~ NA_character_),
      levels = c("First-line TKI", "Subsequent-line TKI"))
  )

# ------------------------------------------- 5. GRADE >=3 ADR (Yes/No) ----
if (!is.na(adr_col)) {
  dat <- dat %>%
    mutate(
      ADR_raw = str_to_lower(str_trim(as.character(.data[[adr_col]]))),
      # Yes/No per spec: only "yes" is Yes; everything else recorded is No.
      Grade3_ADR = factor(
        ifelse(ADR_raw %in% c("yes", "y", "1", "true", "3", "4", "5"),
               "Yes", "No"),
        levels = c("No", "Yes"))
    )
} else {
  dat$ADR_raw <- NA_character_
  dat$Grade3_ADR <- NA
  warning("No Grade >=3 ADR column identified - ADR analyses will be skipped.")
}
has_adr <- !is.na(adr_col)

# --------------------------------------------- 6. SURVIVAL VARIABLES ----
if (is.na(os_time_col) || is.na(os_event_col))
  stop("OS time/event columns not found (need OS_Time_Months / OS_Event).")
if (is.na(pfs_time_col) || is.na(pfs_event_col))
  stop("PFS time/event columns not found (need PFS_Time_Months / PFS_Event).")

dat <- dat %>%
  mutate(
    OS_time  = suppressWarnings(as.numeric(.data[[os_time_col]])),
    OS_event = suppressWarnings(as.numeric(.data[[os_event_col]])),
    PFS_time = suppressWarnings(as.numeric(.data[[pfs_time_col]])),
    PFS_event = suppressWarnings(as.numeric(.data[[pfs_event_col]]))
  )

# ---------------------------------------- 7. EVENT-CODING CHECKS (spec) ----
n_death        <- sum(dat$OS_event == 1, na.rm = TRUE)
n_os_event     <- sum(dat$OS_event == 1, na.rm = TRUE)
n_pfs_event    <- sum(dat$PFS_event == 1, na.rm = TRUE)
n_pfs_death    <- sum(dat$OS_event == 1 & dat$PFS_event == 1, na.rm = TRUE)
n_death_notpfs <- sum(dat$OS_event == 1 & dat$PFS_event != 1, na.rm = TRUE)
n_pfs_nondeath <- sum(dat$OS_event != 1 & dat$PFS_event == 1, na.rm = TRUE)

event_checks <- data.frame(
  Check = c("OS events (= deaths)",
            "PFS events (total)",
            "PFS events that are deaths",
            "Deaths NOT coded as PFS event (MUST be 0)",
            "Non-death PFS events (progression/relapse while alive)"),
  n = c(n_os_event, n_pfs_event, n_pfs_death, n_death_notpfs, n_pfs_nondeath))
cat("\n================ EVENT-CODING CHECK ================\n")
print(event_checks)
if (n_death_notpfs > 0)
  stop("PFS coding error: ", n_death_notpfs,
       " death(s) are NOT coded as a PFS event. Fix the PFS event column.",
       call. = FALSE)
cat("OK: every death is coded as a PFS event.\n")

# ------------------------- 8. ANALYSIS COHORT (LOT-classifiable patients) ----
an <- dat %>% filter(!is.na(LOT))
cat("\n>>> Analysis cohort:", nrow(an), "patients (",
    sum(is.na(dat$LOT)), "with blank TKI excluded )\n")

# ------------------------------------------------------ 9. OUTPUT FOLDER ----
base_dir <- dirname(normalizePath(excel_file))
out_dir <- file.path(base_dir,
                     paste0("CML_Objective2_Results_",
                            format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

save_tbl <- function(x, name) {
  openxlsx::write.xlsx(x, file.path(out_dir, paste0(name, ".xlsx")))
  write.csv(x, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  invisible(x)
}
save_tbl(event_checks, "O3_00_Event_Coding_Checks")

# -------------------------------------------- 10. LOT DESCRIPTIVE TABLES ----
lot_summary <- an %>%
  count(LOT, name = "Frequency") %>%
  mutate(Percentage = round(100 * Frequency / sum(Frequency), 1))

lot_breakdown <- an %>%
  filter(LOT == "Subsequent-line TKI") %>%
  count(recorded_TKI = TKI_clean, name = "n", sort = TRUE) %>%
  mutate(Pct_of_subsequent = round(100 * n / sum(n), 1))

lot_by_factors <- bind_rows(
  an %>% count(LOT, Phase_clean, name = "n") %>%
    mutate(Variable = "CML phase", Level = Phase_clean),
  an %>% count(LOT, Gender_clean, name = "n") %>%
    mutate(Variable = "Gender", Level = Gender_clean)
) %>% select(Variable, Level, LOT, n)

save_tbl(lot_summary, "O3_01a_LOT_Distribution")
save_tbl(lot_breakdown, "O3_01b_LOT_Subsequent_Line_Breakdown")
save_tbl(lot_by_factors, "O3_01c_LOT_by_Baseline_Factors")

# -------------------------------------- 11. GRADE>=3 ADR DESCRIPTIVE TABLES ----
adr_raw_tab <- an %>%
  count(Recorded = coalesce(ADR_raw, "NA"), name = "n") %>%
  mutate(Pct = round(100 * n / sum(n), 1))
adr_summary <- an %>%
  count(Grade3_ADR, name = "Frequency") %>%
  mutate(Percentage = round(100 * Frequency / sum(Frequency), 1))
adr_by_lot <- an %>%
  count(LOT, Grade3_ADR, name = "n") %>%
  group_by(LOT) %>%
  mutate(Pct_of_line = round(100 * n / sum(n), 1)) %>%
  ungroup()

save_tbl(adr_raw_tab, "O3_02a_Grade3ADR_Recorded_Categories")
save_tbl(adr_summary, "O3_02b_Grade3ADR_Binary")
save_tbl(adr_by_lot, "O3_02c_Grade3ADR_by_LOT")

cat("\n>>> LOT distribution:\n");        print(lot_summary)
cat("\n>>> Grade >=3 ADR (Yes/No):\n"); print(adr_summary)

# --------------------------------------------- 12. KAPLAN-MEIER CURVES ----
km_plot <- function(data, time, event, group, fname, title, ylab) {
  d <- data %>% filter(!is.na(.data[[group]]))
  fit <- survfit(as.formula(paste0("Surv(", time, ", ", event, ") ~ ", group)),
                 data = d)
  p <- ggsurvplot(
    fit, data = d, risk.table = TRUE, pval = TRUE, conf.int = TRUE,
    xlab = "Follow-up (months)", ylab = ylab, title = title,
    risk.table.height = 0.22, legend.labs = levels(droplevels(factor(d[[group]]))))
  ggsave(file.path(out_dir, fname),
         plot = gridExtra::arrangeGrob(p$plot, p$table, ncol = 1,
                                       heights = c(3, 1)),
         width = 9, height = 8, dpi = 300)
  lr <- survdiff(as.formula(paste0("Surv(", time, ", ", event, ") ~ ", group)),
                 data = d)
  1 - pchisq(lr$chisq, df = length(lr$n) - 1)
}

# helper: KM summary (n, events, median, 1/3/5-yr %) for one group
one_km <- function(dd, time, event) {
  fit <- survfit(as.formula(paste0("Surv(", time, ", ", event, ") ~ 1")), data = dd)
  med <- tryCatch(quantile(fit, probs = 0.5, conf.int = FALSE)[1],
                  error = function(e) NA_real_)
  st <- function(t) tryCatch(100 * summary(fit, times = t, extend = TRUE)$surv[1],
                             error = function(e) NA_real_)
  data.frame(n = nrow(dd), events = sum(dd[[event]], na.rm = TRUE),
             median_months = round(med, 1),
             S_12mo = round(st(12), 1), S_36mo = round(st(36), 1),
             S_60mo = round(st(60), 1))
}
km_summary_table <- function(data, time, event, group, outcome) {
  d <- data %>% filter(!is.na(.data[[group]]))
  rows <- lapply(c("Overall", levels(droplevels(factor(d[[group]])))), function(g) {
    dd <- if (g == "Overall") d else d %>% filter(.data[[group]] == g)
    r <- one_km(dd, time, event)
    r$Outcome <- outcome; r$Group <- g
    r[, c("Outcome", "Group", "n", "events", "median_months",
          "S_12mo", "S_36mo", "S_60mo")]
  })
  bind_rows(rows)
}

p_os_lot  <- km_plot(an, "OS_time", "OS_event", "LOT",
                     "KM_OS_by_LOT.png",
                     "Overall survival by line of therapy",
                     "Survival probability")
p_pfs_lot <- km_plot(an, "PFS_time", "PFS_event", "LOT",
                     "KM_PFS_by_LOT.png",
                     "Progression-free survival by line of therapy",
                     "Progression-free probability")
p_os_adr <- p_pfs_adr <- NA_real_
if (has_adr) {
  p_os_adr  <- km_plot(an, "OS_time", "OS_event", "Grade3_ADR",
                       "KM_OS_by_Grade3plus_ADR.png",
                       "Overall survival by Grade >=3 ADR",
                       "Survival probability")
  p_pfs_adr <- km_plot(an, "PFS_time", "PFS_event", "Grade3_ADR",
                       "KM_PFS_by_Grade3plus_ADR.png",
                       "Progression-free survival by Grade >=3 ADR",
                       "Progression-free probability")
}

km_tables <- bind_rows(
  km_summary_table(an, "OS_time", "OS_event", "LOT", "OS"),
  km_summary_table(an, "PFS_time", "PFS_event", "LOT", "PFS"))
if (has_adr) {
  km_tables <- bind_rows(
    km_tables,
    km_summary_table(an, "OS_time", "OS_event", "Grade3_ADR", "OS"),
    km_summary_table(an, "PFS_time", "PFS_event", "Grade3_ADR", "PFS"))
}
ps <- c(p_os_lot, p_pfs_lot)
if (has_adr) ps <- c(ps, p_os_adr, p_pfs_adr)
km_tables$logrank_p <- fmt_p(rep(ps, each = 3))
save_tbl(km_tables, "O3_07_KM_Summary_Table")

# --------------------------------------------- 13. UNIVARIABLE COX ----
uni_cox <- function(data, time, event, variable) {
  f <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ", variable))
  d <- data %>% filter(!is.na(.data[[variable]]))
  fit <- coxph(f, data = d)
  tt <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  z <- tryCatch(cox.zph(fit), error = function(e) NULL)
  list(fit = fit,
       table = tt %>%
         mutate(Predictor = variable, N = fit$n, Events = fit$nevent,
                HR = round(estimate, 2),
                CI95 = paste0(round(conf.low, 2), " - ", round(conf.high, 2)),
                P = ifelse(p.value < 0.001, "<0.001",
                           format.pval(p.value, digits = 3, eps = 0.001))) %>%
         select(Predictor, term, HR, CI95, P, N, Events),
       zph = if (!is.null(z)) data.frame(Term = rownames(z$table),
                                         Chi2 = round(z$table[, "chisq"], 2),
                                         P_value = round(z$table[, "p"], 4))
             else NULL)
}

u_lot_os  <- uni_cox(an, "OS_time", "OS_event", "LOT")
u_lot_pfs <- uni_cox(an, "PFS_time", "PFS_event", "LOT")
u_adr_os  <- if (has_adr) uni_cox(an, "OS_time", "OS_event", "Grade3_ADR")  else NULL
u_adr_pfs <- if (has_adr) uni_cox(an, "PFS_time", "PFS_event", "Grade3_ADR") else NULL

univ_all <- bind_rows(u_lot_os$table, u_lot_pfs$table,
                      if (!is.null(u_adr_os))  u_adr_os$table  else NULL,
                      if (!is.null(u_adr_pfs)) u_adr_pfs$table else NULL)
save_tbl(univ_all, "O3_08_Univariate_Cox_LOT_Grade3ADR")

# --------------------------------------- 14. MULTIVARIABLE COX (fixed set) ----
model_data <- an %>%
  mutate(
    Phase_model = factor(
      case_when(
        str_detect(Phase_clean, "chron") ~ "Chronic",
        str_detect(Phase_clean, "accel") ~ "Accelerated",
        str_detect(Phase_clean, "blast") ~ "Blast",
        TRUE ~ NA_character_),
      levels = c("Chronic", "Accelerated", "Blast")),
    Gender_model = factor(Gender_clean)
  )

fit_multi <- function(data, time, event) {
  rhs <- c("Phase_model", "Gender_model", "Age_num", "LOT")
  if (has_adr) rhs <- c(rhs, "Grade3_ADR")
  f <- as.formula(paste0("Surv(", time, ", ", event, ") ~ ",
                         paste(rhs, collapse = " + ")))
  fit <- coxph(f, data = data, na.action = na.omit)
  z <- cox.zph(fit)
  list(fit = fit,
       table = broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
         mutate(HR = round(estimate, 2),
                CI95 = paste0(round(conf.low, 2), " - ", round(conf.high, 2)),
                P = ifelse(p.value < 0.001, "<0.001",
                           format.pval(p.value, digits = 3, eps = 0.001))) %>%
         select(term, HR, CI95, P),
       zph = data.frame(Term = rownames(z$table),
                        Chi2 = round(z$table[, "chisq"], 2),
                        P_value = round(z$table[, "p"], 4)),
       n = fit$n, events = fit$nevent)
}

m_os  <- fit_multi(model_data, "OS_time", "OS_event")
m_pfs <- fit_multi(model_data, "PFS_time", "PFS_event")

save_tbl(m_os$table,  "O3_09a_Multivariate_Cox_OS")
save_tbl(m_pfs$table, "O3_09b_Multivariate_Cox_PFS")

# ------------------------------------- 15. PH TESTS (all models) ----
ph_parts <- list(
  data.frame(Model = "Univariable LOT -> OS",  u_lot_os$zph),
  data.frame(Model = "Univariable LOT -> PFS", u_lot_pfs$zph),
  data.frame(Model = "Multivariable OS",  m_os$zph),
  data.frame(Model = "Multivariable PFS", m_pfs$zph))
if (has_adr) {
  ph_parts <- c(list(
    data.frame(Model = "Univariable G3ADR -> OS",  u_adr_os$zph),
    data.frame(Model = "Univariable G3ADR -> PFS", u_adr_pfs$zph)), ph_parts)
}
ph_all <- bind_rows(ph_parts) %>%
  mutate(PH_OK_at_5pct = ifelse(P_value > 0.05, "yes", "NO"))
save_tbl(ph_all, "O3_10_Proportional_Hazards_Tests")

# --------------------------------------- 16. EXCEL WORKBOOK + TXT SUMMARIES ----
openxlsx::write.xlsx(
  list(Event_Coding_Checks = event_checks,
       LOT_Summary = lot_summary,
       LOT_Subsequent_Breakdown = lot_breakdown,
       LOT_by_Baseline_Factors = lot_by_factors,
       Grade3ADR_Recorded = adr_raw_tab,
       Grade3ADR_Binary = adr_summary,
       Grade3ADR_by_LOT = adr_by_lot,
       KM_Summary = km_tables,
       Univariable_Cox = univ_all,
       Multivariable_Cox_OS = m_os$table,
       Multivariable_Cox_PFS = m_pfs$table,
       PH_Assumption_Tests = ph_all),
  file = file.path(out_dir, "Objective2_Analysis_Tables.xlsx"), overwrite = TRUE)

capture.output(summary(m_os$fit),  file = file.path(out_dir, "Multivariable_OS_Model.txt"))
capture.output(summary(m_pfs$fit), file = file.path(out_dir, "Multivariable_PFS_Model.txt"))
openxlsx::write.xlsx(an, file.path(out_dir, "Derived_CML_Objective2_Dataset.xlsx"))

# ------------------------------------------------- 17. CONSOLE SUMMARY ----
cat("\n============================================================\n")
cat(" OBJECTIVE 2 ANALYSIS COMPLETED\n")
cat("============================================================\n")
cat("Output folder:", out_dir, "\n\n")
cat("Cohort: n =", nrow(an), "| OS events =", sum(an$OS_event == 1),
    "| PFS events =", sum(an$PFS_event == 1), "\n\n")
cat("Univariable Cox (HR, 95% CI, p):\n"); print(univ_all, row.names = FALSE)
cat("\nMultivariable OS (n =", m_os$n, ", events =", m_os$events, "):\n")
print(m_os$table, row.names = FALSE)
cat("\nMultivariable PFS (n =", m_pfs$n, ", events =", m_pfs$events, "):\n")
print(m_pfs$table, row.names = FALSE)
cat("\nProportional-hazards tests:\n"); print(ph_all, row.names = FALSE)
cat("\nGraphs saved: KM_OS_by_LOT.png, KM_PFS_by_LOT.png,")
cat(" KM_OS_by_Grade3plus_ADR.png, KM_PFS_by_Grade3plus_ADR.png\n")
cat("Tables saved: Objective2_Analysis_Tables.xlsx + individual .csv/.xlsx\n")
