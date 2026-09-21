# ================================================================
# CML THESIS - OBJECTIVE 2 SURVIVAL ANALYSIS
# Predictors of survival outcomes among CML patients managed with TKIs
# ================================================================
#
# HOW TO USE
# 1. Install R and RStudio.
# 2. Open this file in RStudio.
# 3. Click "Source".
# 4. A file-upload window will open.
# 5. Select the Excel 5 file (.xlsx).
# 6. The script automatically:
#    - reads the selected Excel workbook
#    - creates First-line vs Subsequent-line LOT
#    - creates Grade >=3 ADR Yes/No
#    - checks OS and PFS variables
#    - runs Kaplan-Meier analyses
#    - runs univariable Cox models
#    - runs multivariable Cox models
#    - produces publication-ready tables
#    - saves graphs and results in an output folder
#
# IMPORTANT:
# The script does NOT invent missing treatment dates or ADR dates.
# If chronological treatment history is available in the workbook,
# the LOT section can be replaced by the chronological algorithm
# provided below.
# ================================================================

required <- c(
  "readxl","dplyr","tidyr","stringr","survival","survminer",
  "ggplot2","broom","openxlsx","scales"
)

missing <- required[!sapply(required, requireNamespace, quietly=TRUE)]
if(length(missing) > 0){
  install.packages(missing, dependencies=TRUE)
}
invisible(lapply(required, library, character.only=TRUE))

# -------------------- 1. SELECT EXCEL FILE -----------------------
excel_file <- file.choose()

# -------------------- 2. READ WORKBOOK ----------------------------
sheets <- readxl::excel_sheets(excel_file)

# Prefer Sheet2, otherwise use the first sheet
sheet_to_use <- if("Sheet2" %in% sheets) "Sheet2" else sheets[1]

dat <- readxl::read_excel(excel_file, sheet=sheet_to_use)

cat("\nSelected file:", basename(excel_file), "\n")
cat("Sheet:", sheet_to_use, "\n")
cat("Rows:", nrow(dat), " Columns:", ncol(dat), "\n\n")

# -------------------- 3. STANDARDISE COLUMN NAMES ----------------
clean_names <- function(x){
  x <- stringr::str_trim(x)
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  x
}
names(dat) <- clean_names(names(dat))

# helper to find columns despite minor naming differences
find_col <- function(candidates, required=TRUE){
  n <- names(dat)
  hit <- n[tolower(n) %in% tolower(candidates)]
  if(length(hit)) return(hit[1])

  # fuzzy matching
  for(z in candidates){
    hit <- n[stringr::str_detect(
      tolower(n), stringr::fixed(tolower(z))
    )]
    if(length(hit)) return(hit[1])
  }

  if(required)
    stop("Required column not found. Tried: ", paste(candidates, collapse=", "))
  return(NA_character_)
}

# -------------------- 4. IDENTIFY IMPORTANT VARIABLES ------------
tki_col <- find_col(c("TKI_Used","TKI_Used_1","TKI"))
phase_col <- find_col(c("Phase","CML_Phase","Disease_Phase"))
gender_col <- find_col(c("Gender","Sex"))
age_col <- find_col(c("Age_clean","Age","Age_years"))

# OS
os_time_col <- find_col(
  c("OS_time","OS_Time","Overall_Survival_Time","OS_months","OS_Months"),
  required=FALSE
)
os_event_col <- find_col(
  c("OS_event","OS_Event","Death_event","Death_Event","Vital_Status"),
  required=FALSE
)

# PFS
pfs_time_col <- find_col(
  c("PFS_time","PFS_Time","Progression_Free_Survival_Time","PFS_months","PFS_Months"),
  required=FALSE
)
pfs_event_col <- find_col(
  c("PFS_event","PFS_Event","Progression_or_Death","PFS_Event_Status"),
  required=FALSE
)

# ADR
adr_col <- find_col(
  c("Grade_3_ADR","Grade_3plus_ADR","Grade_3_or_greater_ADR",
    "Grade_3_5_ADR","ADR_Grade"),
  required=FALSE
)

# -------------------- 5. TKI CLEANING AND LOT --------------------
dat <- dat %>%
  mutate(
    TKI_clean = stringr::str_to_lower(stringr::str_trim(as.character(.data[[tki_col]]))),
    Phase_clean = stringr::str_to_lower(stringr::str_trim(as.character(.data[[phase_col]]))),
    Gender_clean = as.character(.data[[gender_col]]),
    Age_clean_num = suppressWarnings(as.numeric(.data[[age_col]]))
  )

# Canonical TKI names
dat <- dat %>%
  mutate(
    TKI_clean = case_when(
      str_detect(TKI_clean, "imatin") ~ "Imatinib",
      str_detect(TKI_clean, "dasatin") ~ "Dasatinib",
      str_detect(TKI_clean, "nilotin") ~ "Nilotinib",
      str_detect(TKI_clean, "bosutin") ~ "Bosutinib",
      str_detect(TKI_clean, "ponatin") ~ "Ponatinib",
      TRUE ~ stringr::str_to_title(TKI_clean)
    )
  )

# Supervisor-requested LOT grouping
dat <- dat %>%
  mutate(
    LOT = case_when(
      TKI_clean == "Imatinib" ~ "First-line TKI",
      TKI_clean %in% c("Dasatinib","Nilotinib","Bosutinib","Ponatinib") ~
        "Subsequent-line TKI",
      TRUE ~ NA_character_
    ),
    LOT = factor(
      LOT,
      levels=c("First-line TKI","Subsequent-line TKI")
    )
  )

# -------------------- 6. GRADE >=3 ADR ----------------------------
if(!is.na(adr_col)){
  dat <- dat %>%
    mutate(
      ADR_raw = str_to_lower(str_trim(as.character(.data[[adr_col]]))),
      Grade3plus_ADR = case_when(
        ADR_raw %in% c("yes","y","1","true","grade 3","grade 4","grade 5",
                       "grade >=3","grade 3 or greater","3","4","5") ~ "Yes",
        ADR_raw %in% c("no","n","0","false","none","grade 1","grade 2") ~ "No",
        TRUE ~ NA_character_
      ),
      Grade3plus_ADR = factor(Grade3plus_ADR, levels=c("No","Yes"))
    )
} else {
  dat$Grade3plus_ADR <- NA
  warning("No Grade >=3 ADR column was identified. ADR analyses will be skipped.")
}

# -------------------- 7. SURVIVAL VARIABLES ----------------------
# Use existing time/event columns if present.
if(!is.na(os_time_col)){
  dat$OS_time <- suppressWarnings(as.numeric(dat[[os_time_col]]))
} else {
  stop("OS time variable was not found. Add OS_time (months) to the workbook.")
}

if(!is.na(os_event_col)){
  raw <- dat[[os_event_col]]
  dat$OS_event <- case_when(
    as.character(raw) %in% c("1","Yes","YES","yes","Death","Dead","Died") ~ 1,
    as.character(raw) %in% c("0","No","NO","no","Alive","Censored") ~ 0,
    TRUE ~ suppressWarnings(as.numeric(raw))
  )
} else {
  stop("OS event variable was not found. Add OS_event to the workbook.")
}

if(!is.na(pfs_time_col)){
  dat$PFS_time <- suppressWarnings(as.numeric(dat[[pfs_time_col]]))
} else {
  stop("PFS time variable was not found. Add PFS_time (months) to the workbook.")
}

if(!is.na(pfs_event_col)){
  raw <- dat[[pfs_event_col]]
  # CRITICAL: death MUST be an event for PFS.
  dat$PFS_event <- case_when(
    as.character(raw) %in% c("1","Yes","YES","yes","Progression","Death",
                              "Dead","Died","Progression/Death") ~ 1,
    as.character(raw) %in% c("0","No","NO","no","Alive","Censored") ~ 0,
    TRUE ~ suppressWarnings(as.numeric(raw))
  )
} else {
  stop("PFS event variable was not found. Add PFS_event to the workbook.")
}

# Safety check: event values must be 0/1
dat <- dat %>%
  mutate(
    OS_event = ifelse(OS_event %in% c(0,1), OS_event, NA),
    PFS_event = ifelse(PFS_event %in% c(0,1), PFS_event, NA)
  )

# -------------------- 8. DATA QUALITY CHECKS ---------------------
cat("\n================ DATA QUALITY CHECK ================\n")
cat("Total patients:", nrow(dat), "\n")
cat("OS events:", sum(dat$OS_event==1, na.rm=TRUE), "\n")
cat("PFS events:", sum(dat$PFS_event==1, na.rm=TRUE), "\n")
cat("\nLOT:\n")
print(table(dat$LOT, useNA="ifany"))
cat("\nGrade >=3 ADR:\n")
print(table(dat$Grade3plus_ADR, useNA="ifany"))

if(any(dat$PFS_event==0 & dat$OS_event==1, na.rm=TRUE)){
  cat("\nWARNING: Some patients have OS death but PFS event=0.\n")
  cat("Review these records. Death should be coded as a PFS event.\n")
}

# -------------------- 9. OUTPUT FOLDER ----------------------------
base_dir <- dirname(excel_file)
out_dir <- file.path(
  base_dir,
  paste0("CML_Objective2_Results_", format(Sys.time(), "%Y%m%d_%H%M%S"))
)
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# -------------------- 10. LOT SUMMARY -----------------------------
lot_summary <- dat %>%
  count(LOT, name="Frequency") %>%
  mutate(
    Percentage=round(100*Frequency/sum(Frequency),1)
  )

# -------------------- 11. ADR SUMMARY -----------------------------
adr_summary <- dat %>%
  count(Grade3plus_ADR, name="Frequency") %>%
  mutate(
    Percentage=round(100*Frequency/sum(Frequency),1)
  )

# -------------------- 12. KAPLAN-MEIER: LOT -----------------------
km_lot_os <- survfit(Surv(OS_time, OS_event) ~ LOT, data=dat)
km_lot_pfs <- survfit(Surv(PFS_time, PFS_event) ~ LOT, data=dat)

km_lot_os_plot <- ggsurvplot(
  km_lot_os, data=dat, risk.table=TRUE, pval=TRUE,
  conf.int=TRUE, xlab="Time (months)", ylab="Overall survival probability",
  title="Overall Survival by Line of Therapy",
  legend.title="Line of Therapy"
)
ggsave(file.path(out_dir,"KM_OS_by_LOT.png"),
       km_lot_os_plot$plot, width=9, height=7, dpi=300)

km_lot_pfs_plot <- ggsurvplot(
  km_lot_pfs, data=dat, risk.table=TRUE, pval=TRUE,
  conf.int=TRUE, xlab="Time (months)", ylab="Progression-free survival probability",
  title="Progression-Free Survival by Line of Therapy",
  legend.title="Line of Therapy"
)
ggsave(file.path(out_dir,"KM_PFS_by_LOT.png"),
       km_lot_pfs_plot$plot, width=9, height=7, dpi=300)

# -------------------- 13. KAPLAN-MEIER: ADR -----------------------
adr_data_os <- dat %>% filter(!is.na(Grade3plus_ADR))
adr_data_pfs <- dat %>% filter(!is.na(Grade3plus_ADR))

if(nrow(adr_data_os)>0 && length(unique(adr_data_os$Grade3plus_ADR))==2){
  km_adr_os <- survfit(Surv(OS_time, OS_event) ~ Grade3plus_ADR, data=adr_data_os)
  p <- ggsurvplot(
    km_adr_os, data=adr_data_os, risk.table=TRUE, pval=TRUE,
    conf.int=TRUE, xlab="Time (months)",
    ylab="Overall survival probability",
    title="Overall Survival by Grade >=3 ADR Status",
    legend.title="Grade >=3 ADR"
  )
  ggsave(file.path(out_dir,"KM_OS_by_Grade3plus_ADR.png"),
         p$plot, width=9, height=7, dpi=300)
}

if(nrow(adr_data_pfs)>0 && length(unique(adr_data_pfs$Grade3plus_ADR))==2){
  km_adr_pfs <- survfit(Surv(PFS_time, PFS_event) ~ Grade3plus_ADR, data=adr_data_pfs)
  p <- ggsurvplot(
    km_adr_pfs, data=adr_data_pfs, risk.table=TRUE, pval=TRUE,
    conf.int=TRUE, xlab="Time (months)",
    ylab="Progression-free survival probability",
    title="Progression-Free Survival by Grade >=3 ADR Status",
    legend.title="Grade >=3 ADR"
  )
  ggsave(file.path(out_dir,"KM_PFS_by_Grade3plus_ADR.png"),
         p$plot, width=9, height=7, dpi=300)
}

# -------------------- 14. UNIVARIABLE COX -------------------------
uni_cox <- function(data, time, event, variable){
  f <- as.formula(paste0("Surv(",time,",",event,") ~ ",variable))
  fit <- coxph(f, data=data, na.action=na.omit)
  broom::tidy(fit, exponentiate=TRUE, conf.int=TRUE) %>%
    mutate(
      Predictor=variable,
      N=nobs(fit),
      Events=sum(model.response(model.frame(fit))[,2])
    ) %>%
    select(Predictor, term, estimate, conf.low, conf.high, p.value, N, Events)
}

uni_lot_os <- uni_cox(dat,"OS_time","OS_event","LOT")
uni_lot_pfs <- uni_cox(dat,"PFS_time","PFS_event","LOT")

uni_adr_os <- if(sum(!is.na(dat$Grade3plus_ADR))>0)
  uni_cox(dat,"OS_time","OS_event","Grade3plus_ADR") else data.frame()

uni_adr_pfs <- if(sum(!is.na(dat$Grade3plus_ADR))>0)
  uni_cox(dat,"PFS_time","PFS_event","Grade3plus_ADR") else data.frame()

# -------------------- 15. MULTIVARIABLE COX -----------------------
# CML phase is retained as a clinical predictor.
# Gender and age are also adjusted for.
model_data <- dat %>%
  mutate(
    Phase_model = case_when(
      str_detect(Phase_clean,"chronic") ~ "Chronic",
      str_detect(Phase_clean,"accelerated") ~ "Accelerated",
      str_detect(Phase_clean,"blast") ~ "Blast",
      TRUE ~ NA_character_
    ),
    Phase_model=factor(Phase_model, levels=c("Chronic","Accelerated","Blast")),
    Gender_model=factor(Gender_clean)
  )

# OS
cox_os <- coxph(
  Surv(OS_time,OS_event) ~ Phase_model + Gender_model +
    Age_clean_num + LOT + Grade3plus_ADR,
  data=model_data, na.action=na.omit
)

# PFS
cox_pfs <- coxph(
  Surv(PFS_time,PFS_event) ~ Phase_model + Gender_model +
    Age_clean_num + LOT + Grade3plus_ADR,
  data=model_data, na.action=na.omit
)

multi_os <- broom::tidy(cox_os, exponentiate=TRUE, conf.int=TRUE)
multi_pfs <- broom::tidy(cox_pfs, exponentiate=TRUE, conf.int=TRUE)

# -------------------- 16. PROPORTIONAL HAZARDS CHECK --------------
ph_os <- cox.zph(cox_os)
ph_pfs <- cox.zph(cox_pfs)

# -------------------- 17. TABLES ------------------------------
format_cox <- function(x){
  x %>%
    mutate(
      HR=round(estimate,3),
      CI_95=paste0(round(conf.low,3),"–",round(conf.high,3)),
      P_value=ifelse(p.value<0.001,"<0.001",
                     format.pval(p.value,digits=3,eps=0.001))
    ) %>%
    select(term,HR,CI_95,P_value)
}

tables <- list(
  LOT_Summary=lot_summary,
  Grade3plus_ADR_Summary=adr_summary,
  Univariable_Cox_LOT_OS=format_cox(uni_lot_os),
  Univariable_Cox_LOT_PFS=format_cox(uni_lot_pfs),
  Univariable_Cox_ADR_OS=format_cox(uni_adr_os),
  Univariable_Cox_ADR_PFS=format_cox(uni_adr_pfs),
  Multivariable_Cox_OS=format_cox(multi_os),
  Multivariable_Cox_PFS=format_cox(multi_pfs),
  PH_Assumption_OS=as.data.frame(ph_os$table),
  PH_Assumption_PFS=as.data.frame(ph_pfs$table)
)

# -------------------- 18. EXPORT EXCEL RESULTS -------------------
openxlsx::write.xlsx(
  tables,
  file=file.path(out_dir,"Objective2_Analysis_Tables.xlsx"),
  overwrite=TRUE
)

# Save the derived dataset
openxlsx::write.xlsx(
  dat,
  file=file.path(out_dir,"Derived_CML_Objective2_Dataset.xlsx"),
  overwrite=TRUE
)

# Save model summaries
capture.output(summary(cox_os),
               file=file.path(out_dir,"Multivariable_OS_Model.txt"))
capture.output(summary(cox_pfs),
               file=file.path(out_dir,"Multivariable_PFS_Model.txt"))
capture.output(ph_os,
               file=file.path(out_dir,"PH_Assumption_OS.txt"))
capture.output(ph_pfs,
               file=file.path(out_dir,"PH_Assumption_PFS.txt"))

# -------------------- 19. CONSOLE SUMMARY -------------------------
cat("\n============================================================\n")
cat("OBJECTIVE 2 ANALYSIS COMPLETED\n")
cat("============================================================\n")
cat("Output folder:\n",out_dir,"\n\n")

cat("LOT distribution:\n")
print(lot_summary)

cat("\nGrade >=3 ADR distribution:\n")
print(adr_summary)

cat("\nUnivariable LOT -> OS:\n")
print(format_cox(uni_lot_os))

cat("\nUnivariable LOT -> PFS:\n")
print(format_cox(uni_lot_pfs))

cat("\nUnivariable Grade >=3 ADR -> OS:\n")
print(format_cox(uni_adr_os))

cat("\nUnivariable Grade >=3 ADR -> PFS:\n")
print(format_cox(uni_adr_pfs))

cat("\nMultivariable OS:\n")
print(format_cox(multi_os))

cat("\nMultivariable PFS:\n")
print(format_cox(multi_pfs))

cat("\nPH test - OS:\n")
print(ph_os)

cat("\nPH test - PFS:\n")
print(ph_pfs)

cat("\nAll tables, graphs and derived data have been saved in the output folder.\n")
