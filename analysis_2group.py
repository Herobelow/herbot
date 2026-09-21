#!/usr/bin/env python3
# =============================================================================
#  CML TKI STUDY (KNH) - REGROUPED ANALYSIS, trial 5.xlsx
#
#  Regrouping requested by the supervisor (Sep 2026):
#    * 2nd-generation TKIs are analysed TOGETHER (never separately -
#      numbers are too small).
#    * PRIMARY comparison = TWO GROUPS ONLY:
#          FIRST LINE      = recorded imatinib          (n = 97)
#          SUBSEQUENT LINE = recorded 2G TKI (dasatinib / nilotinib /
#                            bosutinib) or ponatinib    (n = 65)
#      Rationale: only 5 patients were recorded on ponatinib, so a
#      separate "3rd-line / 3G" group is too small for any meaningful
#      analysis. Ponatinib patients are pooled into "subsequent line".
#
#  Line of therapy is INFERRED from the recorded (current) TKI: every
#  patient on a 2G/3G TKI is presumed to have received imatinib first
#  (see Fix_Data_IMatinib_Assumption.R). Later-line outcomes are
#  confounded by indication and must not be read as drug inferiority.
#
#  Outputs -> ./output/  (CSV tables + PNG figures + summary MD)
# =============================================================================

import os
import contextlib
import warnings

import numpy as np
import pandas as pd

from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test

warnings.filterwarnings("ignore")

XLSX = "trial 5.xlsx"
OUT = "output"
os.makedirs(OUT, exist_ok=True)
T_FIVE_YEARS = 60.0
VAR_ENTRY_P = 0.20
EPV_MIN = 10

saved = []


def save_csv(df, name):
    p = os.path.join(OUT, name)
    df.to_csv(p, index=False)
    saved.append(name)
    print(f">>> saved table:  {name}")


def save_fig(fig, name):
    p = os.path.join(OUT, name)
    fig.savefig(p, dpi=300, bbox_inches="tight")
    plt.close(fig)
    saved.append(name)
    print(f">>> saved figure: {name}")


def fmt_p(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "NA"
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


def median_iqr(x):
    x = pd.Series(x).dropna()
    if len(x) == 0:
        return "NA"
    q1, q3 = np.percentile(x, [25, 75])
    return f"{np.median(x):.1f} ({q1:.1f} - {q3:.1f})"


def n_pct(x):
    x = pd.Series(x).dropna()
    n_known = len(x)
    if n_known == 0:
        return "NA"
    n_yes = int(x.astype(bool).sum()) if x.dtype != object else int((x == x.mode()[0]).sum())
    return f"{n_yes} ({100 * n_yes / n_known:.1f}%)"


def flag_n_pct(x):
    x = pd.Series(x).dropna()
    n_known = len(x)
    if n_known == 0:
        return "NA"
    n_yes = int(x.sum())
    return f"{n_yes} ({100 * n_yes / n_known:.1f}%)"


def factor_dist(x):
    x = pd.Series(x).dropna()
    if len(x) == 0:
        return "NA"
    vc = x.value_counts()
    return " | ".join(f"{k}: {v} ({100 * v / len(x):.1f}%)" for k, v in vc.items())


def cat_pvalue(tab):
    """Fisher for 2x2 or small cells, else chi-square."""
    tab = np.asarray(tab)
    tab = tab[tab.sum(axis=1) > 0][:, tab.sum(axis=0) > 0]
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return np.nan
    if tab.shape == (2, 2):
        return stats.fisher_exact(tab)[1]
    if (tab < 5).any():
        try:
            res = stats.fisher_exact(tab, alternative="two-sided") if tab.shape == (2, 2) else None
        except Exception:
            res = None
        if res is not None:
            return res[1]
        try:
            return stats.chi2_contingency(tab, correction=False)[1]
        except ValueError:
            return np.nan
    try:
        return stats.chi2_contingency(tab, correction=False)[1]
    except ValueError:
        return np.nan


def prop_ci(x, n):
    if n < 1:
        return np.nan, np.nan
    if x == 0:
        return 0.0, min(1.0, 3 / n)
    if x >= n:
        return max(0.0, 1 - 3 / n), 1.0
    p = x / n
    z = 1.959964
    se = np.sqrt(p * (1 - p) / n)
    return max(0.0, p - z * se), min(1.0, p + z * se)


# =============================================================================
# 1. LOAD + DERIVE VARIABLES
# =============================================================================
print(">>> Loading", XLSX)
raw = pd.read_excel(XLSX, sheet_name="Sheet2")
raw.columns = [str(c).strip() for c in raw.columns]
df = raw.copy()
print(f">>> Loaded {df.shape[0]} rows x {df.shape[1]} columns")

# raw 'Event' duplicates OS_Event (identical values) - drop to avoid a
# duplicated column name after renaming below
df = df.drop(columns=["Event"], errors="ignore")

# --- recorded TKI + generation (descriptive only) --------------------------
df["tki_clean"] = df["TKI Used"].astype(str).str.strip().str.lower()
df.loc[df["tki_clean"].isin(["", "nan", "none"]), "tki_clean"] = np.nan

GEN = {"imatinib": "1G", "dasatinib": "2G", "nilotinib": "2G",
       "bosutinib": "2G", "ponatinib": "3G"}
df["TKI_Gen"] = df["tki_clean"].map(GEN)

# --- PRIMARY GROUPING: first line vs subsequent line -----------------------
line_map = {"1G": "First line (imatinib)"}
df["LINE"] = df["TKI_Gen"].map(line_map)
df.loc[df["TKI_Gen"].isin(["2G", "3G"]), "LINE"] = "Subsequent line (2G/3G TKI)"
LINE_LEVELS = ["First line (imatinib)", "Subsequent line (2G/3G TKI)"]

# --- time-to-event ----------------------------------------------------------
df["OS_Time"] = pd.to_numeric(df["OS_Time_Months"], errors="coerce")
df["OS_Event"] = pd.to_numeric(df["OS_Event"], errors="coerce").astype("Int64")
df["PFS_Time"] = pd.to_numeric(df["PFS_Time_Months"], errors="coerce")
df["PFS_Event"] = pd.to_numeric(df["PFS_Event"], errors="coerce").astype("Int64")

# --- demographics -------------------------------------------------------------
df["Age_num"] = pd.to_numeric(df["Age"], errors="coerce")
df["Age_group"] = pd.cut(df["Age_num"], [-np.inf, 39, 59, np.inf],
                         labels=["<40", "40-59", ">=60"])
df["Sex"] = df["Gender"].astype(str).str.strip()
df.loc[~df["Sex"].isin(["Male", "Female"]), "Sex"] = np.nan

ph = df["Phase of CML"].astype(str).str.strip().str.lower()
df["Phase_group"] = pd.Series(np.nan, index=df.index, dtype=object)
df.loc[ph.str.contains("chron", na=False), "Phase_group"] = "Chronic"
df.loc[ph.str.contains("accel", na=False), "Phase_group"] = "Accelerated"
df.loc[ph.str.contains("blast", na=False), "Phase_group"] = "Blast"
# explicit category order so Chronic is the Cox reference level
df["Phase_group"] = pd.Categorical(
    df["Phase_group"], categories=["Chronic", "Accelerated", "Blast"])

# --- clinical flags -----------------------------------------------------------
def to_flag(series, positives=("yes", "y", "1", "true", "present")):
    # Same semantics as the R script's to_flag(): blank/NA stays NA;
    # any other text (incl. 'not documented' / 'not applicable') counts
    # as FALSE unless it is a positive word.
    v = series.astype(str).str.strip().str.lower()
    v = v.replace({"nan": np.nan, "": np.nan})
    out = pd.Series(False, index=series.index)
    out[v.isin(positives)] = True
    return out.astype("boolean").where(v.notna(), pd.NA)

df["flag_T315I"] = to_flag(df["T315I_Mutation"])
df["flag_ADR_any"] = to_flag(df["Any ADRs"])
df["flag_ADR_sev"] = to_flag(df["Grade_3_ADR"])
df["flag_Relapse"] = to_flag(df["Disease_Relapse"])
df["flag_Comorbidity"] = to_flag(df["Any Comorbidity"])

# NOTE: the 'Disease Progression' column in this dataset actually holds
# Alive/Dead status, not a progression flag, so it is NOT used as a
# covariate (it duplicates Survival_Status).

# BCR-ABL monitoring: any recorded value (incl. numeric transcript levels)
# = monitored; 'not done' / 'not documented' / blank = not monitored.
b = df["BCR_ABL_Monitoring"].astype(str).str.strip().str.lower()
b = b.replace({"nan": np.nan, "": np.nan})
df["BCR_monitored"] = pd.Series(np.nan, index=df.index, dtype=object)
df.loc[b.isin(["not done", "not documented", "no", "0", "none"]), "BCR_monitored"] = "No"
df.loc[b.notna() & ~b.isin(["not done", "not documented", "no", "0", "none"]),
       "BCR_monitored"] = "Yes"

modif = df["TKI_Treatment_Modification"].astype(str).str.strip().str.lower()
modif = modif.replace({"nan": ""})
df["mod_switch"] = modif.str.contains("switch")
df["any_modification"] = modif.str.contains(
    "switch|interruption|reduction|adjust|discontin")
df["died"] = df["OS_Event"] == 1

# --- analysis data sets -------------------------------------------------------
base_covs = ["Study ID", "tki_clean", "TKI_Gen", "LINE", "Age_num", "Age_group",
             "Sex", "Phase_group", "flag_T315I", "flag_ADR_any", "flag_ADR_sev",
             "flag_Relapse", "flag_Comorbidity", "BCR_monitored", "mod_switch"]

os_data = df.dropna(subset=["OS_Time", "OS_Event"]).copy()
os_data = os_data.rename(columns={"OS_Time": "Time", "OS_Event": "Event"})
pfs_data = df.dropna(subset=["PFS_Time", "PFS_Event"]).copy()
pfs_data = pfs_data.rename(columns={"PFS_Time": "Time", "PFS_Event": "Event"})

n_class = df["LINE"].notna().sum()
print(f">>> Cohort: {len(df)} patients | classifiable by line: {n_class}")
print(f"    OS complete: {len(os_data)} | PFS complete: {len(pfs_data)}")
print(">>> Recorded TKI distribution:")
print(df["tki_clean"].value_counts(dropna=False).to_string())
print(">>> PRIMARY grouping:")
print(df["LINE"].value_counts(dropna=False).to_string())

# =============================================================================
# 2. GROUPING TABLES (Objective 1)
# =============================================================================
o1_01 = (df.groupby("tki_clean", dropna=False).size().reset_index(name="n")
         .rename(columns={"tki_clean": "TKI"}))
o1_01["Generation"] = o1_01["TKI"].map(GEN)
o1_01["LINE_group"] = np.where(o1_01["Generation"] == "1G", LINE_LEVELS[0],
                               np.where(o1_01["Generation"].isin(["2G", "3G"]),
                                        LINE_LEVELS[1], "Unclassified (blank TKI)"))
o1_01["Pct_of_cohort"] = (100 * o1_01["n"] / len(df)).round(1)
o1_01 = o1_01.sort_values("n", ascending=False)
save_csv(o1_01, "O1_01_TKI_Distribution_and_Grouping.csv")

lot_tab = (df.groupby("LINE").size().reset_index(name="n"))
lot_tab["Percent_of_classifiable"] = (100 * lot_tab["n"] / n_class).round(1)
lot_tab["LINE"] = pd.Categorical(lot_tab["LINE"], LINE_LEVELS)
lot_tab = lot_tab.sort_values("LINE")
save_csv(lot_tab, "O1_14a_Line_Group_Categorisation.csv")

lot_break = (df[df["LINE"] == LINE_LEVELS[1]]
             .groupby(["TKI_Gen", "tki_clean"]).size().reset_index(name="n")
             .rename(columns={"tki_clean": "recorded_TKI"}))
lot_break["Percent_of_subsequent_line"] = (
    100 * lot_break["n"] / lot_break["n"].sum()).round(1)
lot_break = lot_break.sort_values("n", ascending=False)
save_csv(lot_break, "O1_14b_Subsequent_Line_Breakdown.csv")

# ponatinib detail (n=5 -> why pooling is necessary)
pon = df[df["tki_clean"] == "ponatinib"][
    ["Study ID", "Age_num", "Phase_group", "tki_clean", "LINE",
     "OS_Event", "OS_Time"]].rename(
    columns={"Age_num": "Age", "tki_clean": "TKI", "OS_Time": "OS_Time_Months"})
pon["Survival"] = np.where(pon["OS_Event"] == 1, "Deceased", "Alive")
save_csv(pon, "O1_10_Ponatinib_Patients_n5.csv")

# =============================================================================
# 3. TABLE 1 - baseline characteristics by line group (2 groups)
# =============================================================================
d_grp = df[df["LINE"].notna()].copy()

def t1_row(variable, kind, col, levels=None):
    row = {"Variable": variable}
    for g in LINE_LEVELS + ["Overall"]:
        sub = d_grp[col] if g == "Overall" else d_grp.loc[d_grp["LINE"] == g, col]
        if kind == "num":
            row[g] = median_iqr(sub)
        elif kind == "flag":
            row[g] = flag_n_pct(sub)
        else:
            row[g] = factor_dist(sub)
    if kind == "num":
        a = [d_grp.loc[d_grp["LINE"] == g, col].dropna().values for g in LINE_LEVELS]
        a = [x for x in a if len(x) >= 1]
        row["P_value"] = stats.kruskal(*a).pvalue if len(a) == 2 else np.nan
    else:
        sub2 = d_grp[[col, "LINE"]].dropna()
        ct = pd.crosstab(sub2["LINE"], sub2[col])
        row["P_value"] = cat_pvalue(ct.values)
    return row

rows = []
rows.append(t1_row("Age (years)", "num", "Age_num"))
rows.append(t1_row("Age group", "factor", "Age_group"))
rows.append(t1_row("Sex", "factor", "Sex"))
rows.append(t1_row("CML phase at diagnosis", "factor", "Phase_group"))
rows.append(t1_row("T315I mutation (Yes)", "flag", "flag_T315I"))
rows.append(t1_row("Any ADR", "flag", "flag_ADR_any"))
rows.append(t1_row("Grade 3+ ADR", "flag", "flag_ADR_sev"))
rows.append(t1_row("Disease relapse", "flag", "flag_Relapse"))
rows.append(t1_row("Any comorbidity", "flag", "flag_Comorbidity"))
rows.append(t1_row("BCR-ABL monitoring done", "factor", "BCR_monitored"))
rows.append(t1_row("Documented TKI dose switch", "flag", "mod_switch"))

tbl1 = pd.DataFrame(rows)
tbl1["P_value"] = tbl1["P_value"].round(4)
tbl1 = tbl1[["Variable", LINE_LEVELS[0], LINE_LEVELS[1], "Overall", "P_value"]]
tbl1.columns = ["Variable", "First line (imatinib)", "Subsequent line (2G/3G TKI)",
                "Overall", "P_value"]
save_csv(tbl1, "Table1_Baseline_by_Line_Group.csv")
print("\n>>> Table 1 (baseline by line group):")
print(tbl1.to_string(index=False))

# phase x line cross-tab (panel comment ii, now on the 2-group split)
o1_02 = (d_grp.dropna(subset=["Phase_group"])
         .groupby(["LINE", "Phase_group"]).size().reset_index(name="n"))
o1_02["Pct_of_group"] = o1_02.groupby("LINE")["n"].transform(
    lambda s: (100 * s / s.sum()).round(1))
save_csv(o1_02, "O1_02_Phase_by_Line_Group.csv")
o1_02_p = pd.DataFrame([{
    "Test": "Phase of CML vs line group (first line vs subsequent)",
    "P_value": cat_pvalue(pd.crosstab(d_grp["LINE"], d_grp["Phase_group"]).values),
}])
o1_02_p["P_report"] = o1_02_p["P_value"].map(fmt_p)
save_csv(o1_02_p, "O1_02_Phase_by_Line_Group_pvalue.csv")

# treatment modifications by line group
o1_03 = d_grp.groupby("LINE").agg(
    n=("Study ID", "size"),
    n_switch=("mod_switch", "sum"),
    n_any_mod=("any_modification", "sum")).reset_index()
o1_03["Pct_any_mod"] = (100 * o1_03["n_any_mod"] / o1_03["n"]).round(1)
save_csv(o1_03, "O1_03_Treatment_Modifications_by_Line_Group.csv")

# =============================================================================
# 4. OBJECTIVE 2 - SURVIVAL BY LINE GROUP (primary analysis)
# =============================================================================
def km_at(kmf, t):
    """KM estimate at t with log-transformed 95% CI (Greenwood)."""
    sf = kmf.survival_function_.iloc[:, 0]
    if t < sf.index.min():
        return 1.0, 1.0, 1.0
    s = float(sf.loc[:t].iloc[-1])
    if s <= 0 or np.isnan(s):
        return np.nan, np.nan, np.nan
    ci = kmf.confidence_interval_survival_function_
    lo = float(ci.iloc[:, 0].loc[:t].iloc[-1])
    hi = float(ci.iloc[:, 1].loc[:t].iloc[-1])
    return s, lo, hi


def km_rmst(kmf, tau):
    try:
        from lifelines.utils import restricted_mean_survival_time
        return float(restricted_mean_survival_time(kmf, tau))
    except Exception:
        sf = kmf.survival_function_.iloc[:, 0]
        times = np.concatenate([[0.0], sf.index.values, [tau]])
        surv = np.concatenate([[1.0], sf.values, [sf.values[-1]]])
        keep = times <= tau
        return float(np.trapezoid(surv[keep], times[keep]))


def group_km_table(d, outcome_label):
    rows = []
    overall_kmf = KaplanMeierFitter().fit(d["Time"], d["Event"])
    specs = [("Overall", overall_kmf, len(d), int(d["Event"].sum()))]
    for g in LINE_LEVELS:
        dg = d[d["LINE"] == g]
        kmf = KaplanMeierFitter().fit(dg["Time"], dg["Event"], label=g)
        specs.append((g, kmf, len(dg), int(dg["Event"].sum())))
    tmax = float(d["Time"].max())
    for name, kmf, n, ev in specs:
        med = kmf.median_survival_time_
        row = {
            "Group": name, "Outcome": outcome_label, "N": n,
            "Events_n": ev,
            "Events_pct": round(100 * ev / n, 1) if n else np.nan,
            "Median_months": np.nan if np.isinf(med) else round(float(med), 1),
        }
        for lab, t in [("S_12mo", 12), ("S_36mo", 36), ("S_60mo", 60)]:
            s, lo, hi = km_at(kmf, t)
            row[lab] = np.nan if np.isnan(s) else round(100 * s, 1)
            row[lab + "_CI95"] = ("NA" if np.isnan(s)
                                  else f"{100 * lo:.1f} - {100 * hi:.1f}")
        row["RMST_60mo_months"] = round(km_rmst(kmf, min(T_FIVE_YEARS, tmax)), 1)
        row["Max_followup_months"] = round(tmax, 1)
        rows.append(row)
    return pd.DataFrame(rows)


def km_plot(d, outcome_label, ylabel, fname, title):
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(9, 6.6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.06)
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax)
    colors = {LINE_LEVELS[0]: "#1f77b4", LINE_LEVELS[1]: "#d62728"}
    for g in LINE_LEVELS:
        dg = d[d["LINE"] == g]
        kmf = KaplanMeierFitter().fit(dg["Time"], dg["Event"], label=g)
        kmf.plot_survival_function(ax=ax, ci_show=True, color=colors[g],
                                   linewidth=2)
    lr = multivariate_logrank_test(d["Time"], d["LINE"], d["Event"])
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="lower left", frameon=False)
    ax.text(0.98, 0.97, f"Log-rank p {fmt_p(lr.p_value)}", transform=ax.transAxes,
            fontsize=11, ha="right", va="top",
            bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9))
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.25)
    plt.setp(ax.get_xticklabels(), visible=False)
    # --- number-at-risk table ---
    times = np.arange(0, min(181, d["Time"].max()) + 1, 24)
    for i, g in enumerate(LINE_LEVELS):
        dg = d[d["LINE"] == g]
        risk = [(dg["Time"] >= t).sum() for t in times]
        y = 1 - 0.45 * i
        ax2.text(-0.02, y, "First line" if i == 0 else "Subsequent",
                 transform=ax2.get_yaxis_transform(), ha="right", va="center",
                 fontsize=9, color=colors[g])
        for t, r in zip(times, risk):
            ax2.text(t, y, str(r), ha="center", va="center", fontsize=8,
                     color=colors[g])
    ax2.set_yticks([])
    ax2.set_ylim(0, 1)
    ax2.set_xlim(0, min(181, d["Time"].max()))
    ax2.set_xticks(times)
    ax2.set_xlabel("Follow-up (months)")
    ax2.set_title("Number at risk", fontsize=9, loc="left")
    for s in ["top", "right", "left"]:
        ax2.spines[s].set_visible(False)
    save_fig(fig, fname)
    return lr.p_value


print("\n>>> Objective 2: survival by line group...")
os_grp = os_data[os_data["LINE"].notna()]
pfs_grp = pfs_data[pfs_data["LINE"].notna()]

km_os = group_km_table(os_grp, "OS")
km_pfs = group_km_table(pfs_grp, "PFS")
km_tbl = pd.concat([km_os, km_pfs], ignore_index=True)
save_csv(km_tbl, "O2_01_KM_by_Line_Group.csv")
print(km_tbl.to_string(index=False))

p_os_lr = km_plot(os_grp, "OS", "Probability of survival",
                  "O2_02_KM_OS_by_Line_Group.png",
                  "Overall survival by line of therapy\n(first line vs subsequent line, trial 5)")
p_pfs_lr = km_plot(pfs_grp, "PFS", "Probability of being progression-free",
                   "O2_03_KM_PFS_by_Line_Group.png",
                   "Progression-free survival by line of therapy\n(first line vs subsequent line, trial 5)")

lr_tbl = pd.DataFrame([
    {"Outcome": "OS", "Comparison": "First line vs subsequent line",
     "N": len(os_grp), "P_value": p_os_lr, "P_report": fmt_p(p_os_lr)},
    {"Outcome": "PFS", "Comparison": "First line vs subsequent line",
     "N": len(pfs_grp), "P_value": p_pfs_lr, "P_report": fmt_p(p_pfs_lr)},
])
save_csv(lr_tbl, "O2_04_LogRank_by_Line_Group.csv")

# --- PH-free effect measure: difference in restricted mean survival (60 mo) ---
# NOTE: lifelines' analytic RMST variance is numerically unstable on this
# censored tail (IntegrationWarning, absurd CIs); use a seeded bootstrap.
from lifelines.utils import restricted_mean_survival_time


def rmst60(times, events):
    kmf = KaplanMeierFitter().fit(times, events)
    return float(restricted_mean_survival_time(kmf, T_FIVE_YEARS))


rmst_rows = []
boot_rng = np.random.default_rng(20260921)
B = 1000
for d, out in [(os_grp, "OS"), (pfs_grp, "PFS")]:
    vals = {}
    for g in LINE_LEVELS:
        dg = d[d["LINE"] == g]
        vals[g] = rmst60(dg["Time"], dg["Event"])
    diff = vals[LINE_LEVELS[0]] - vals[LINE_LEVELS[1]]
    boot = []
    for _ in range(B):
        b0 = d[d["LINE"] == LINE_LEVELS[0]].sample(frac=1, replace=True, random_state=boot_rng)
        b1 = d[d["LINE"] == LINE_LEVELS[1]].sample(frac=1, replace=True, random_state=boot_rng)
        boot.append(rmst60(b0["Time"], b0["Event"]) - rmst60(b1["Time"], b1["Event"]))
    boot = np.asarray(boot)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    p_boot = 2 * min((boot <= 0).mean(), (boot >= 0).mean())
    rmst_rows.append({
        "Outcome": out,
        "RMST60_first_line_mo": round(vals[LINE_LEVELS[0]], 1),
        "RMST60_subsequent_line_mo": round(vals[LINE_LEVELS[1]], 1),
        "Difference_mo": round(diff, 1),
        "CI95_boot": f"{lo:.1f} - {hi:.1f}",
        "P_boot": p_boot,
    })
rmst_df = pd.DataFrame(rmst_rows)
rmst_df["P_report"] = rmst_df["P_boot"].map(fmt_p)
save_csv(rmst_df, "O2_07_RMST_Difference_by_Line.csv")

# overall KM (all patients) for reference
for d, lab, fname in [(os_data, "Overall survival", "O2_05_KM_OS_Overall.png"),
                      (pfs_data, "Progression-free survival", "O2_06_KM_PFS_Overall.png")]:
    fig, ax = plt.subplots(figsize=(8, 5.5), constrained_layout=True)
    kmf = KaplanMeierFitter().fit(d["Time"], d["Event"])
    kmf.plot_survival_function(ax=ax, ci_show=True, color="#2ca02c", linewidth=2)
    ax.set_xlabel("Follow-up (months)")
    ax.set_ylabel("Probability")
    ax.set_title(f"{lab} - all CML patients on TKIs (KNH cohort, n={len(d)})")
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.25)
    save_fig(fig, fname)

# =============================================================================
# 5. OBJECTIVE 3 - COX REGRESSION (2-group design)
# =============================================================================
# Candidates: TKI generation is NOT included (collinear with line group and
# the 3G cell has only 5 patients). LINE (first vs subsequent) is the
# primary exposure.
COX_CANDIDATES = ["Phase_group", "Age_group", "Sex", "LINE", "flag_T315I",
                  "flag_ADR_any", "flag_ADR_sev", "flag_Relapse",
                  "flag_Comorbidity", "BCR_monitored", "mod_switch"]
FORCE_VARS = ["flag_ADR_sev", "LINE"]

VAR_LABELS = {
    "Phase_group": "CML phase (ref Chronic)",
    "Age_group": "Age group (ref <40)",
    "Sex": "Sex (Male vs Female)",
    "LINE": "Line of therapy (subsequent vs first)",
    "flag_T315I": "T315I mutation",
    "flag_ADR_any": "Any ADR",
    "flag_ADR_sev": "Grade 3+ ADR",
    "flag_Relapse": "Disease relapse",
    "flag_Comorbidity": "Any comorbidity",
    "BCR_monitored": "BCR-ABL monitoring done",
    "mod_switch": "Documented TKI dose switch",
}


def cox_ready(d):
    x = d.copy()
    x["LINE"] = (x["LINE"] == LINE_LEVELS[1]).astype(float)
    x["Sex"] = (x["Sex"] == "Male").astype(float).where(x["Sex"].notna(), np.nan)
    for c in ["flag_T315I", "flag_ADR_any", "flag_ADR_sev", "flag_Relapse",
              "flag_Comorbidity", "mod_switch"]:
        x[c] = x[c].astype("boolean").astype(float)
    x["BCR_monitored"] = ((x["BCR_monitored"] == "Yes").astype(float)
                          .where(x["BCR_monitored"].notna(), np.nan))
    # reference levels: Chronic (phase) and <40 (age) are dropped
    x["Age_group"] = pd.Categorical(
        x["Age_group"], categories=["<40", "40-59", ">=60"])
    x = pd.get_dummies(x, columns=["Phase_group", "Age_group"],
                       drop_first=True, dtype=float)
    return x


def usable(d, var):
    if var in ("Phase_group", "Age_group"):
        s = d[var].dropna()
        return len(s) >= 10 and s.nunique() >= 2
    s = d[var].dropna()
    return len(s) >= 10 and s.nunique() >= 2


def univariable_cox(d, outcome):
    rows = []
    dx = cox_ready(d)
    dummy_map = {}
    for v in COX_CANDIDATES:
        if v in ("Phase_group", "Age_group"):
            dummy_map[v] = [c for c in dx.columns if c.startswith(v + "_")]
        else:
            dummy_map[v] = [v]
    for v in COX_CANDIDATES:
        if not usable(d, v):
            continue
        terms = dummy_map[v]
        sub = dx[["Time", "Event"] + terms].dropna()
        if len(sub) < 10 or sub["Event"].sum() < 5:
            continue
        # drop zero-variance dummy levels
        terms2 = [t for t in terms if sub[t].nunique() > 1]
        if not terms2:
            continue
        sub = sub.drop(columns=[t for t in terms if t not in terms2])
        try:
            cph = CoxPHFitter(penalizer=1e-4)
            cph.fit(sub, "Time", "Event")
        except Exception as e:
            print(f"    cox failed for {v}: {e}")
            continue
        s = cph.summary
        for term in terms2:
            r = s.loc[term]
            rows.append({
                "Outcome": outcome, "Covariate": v, "Label": VAR_LABELS[v],
                "Term": term,
                "HR": round(float(np.exp(r["coef"])), 2),
                "CI95_lower": round(float(np.exp(r["coef lower 95%"])), 2),
                "CI95_upper": round(float(np.exp(r["coef upper 95%"])), 2),
                "P_value": float(r["p"]),
                "P_report": fmt_p(float(r["p"])),
                "Report": (f"HR {np.exp(r['coef']):.2f} "
                           f"(95% CI {np.exp(r['coef lower 95%']):.2f} - "
                           f"{np.exp(r['coef upper 95%']):.2f}), "
                           f"p = {fmt_p(float(r['p']))}"),
            })
    univ = pd.DataFrame(rows)
    # minimum p per covariate (for selection)
    if len(univ):
        univ_min = univ.groupby("Covariate")["P_value"].min()
    else:
        univ_min = pd.Series(dtype=float)
    return univ, univ_min


def fit_block(d, terms):
    try:
        cph = CoxPHFitter(penalizer=1e-4)
        cph.fit(d[["Time", "Event"] + terms].dropna(), "Time", "Event")
        return cph
    except Exception:
        return None


def multivariable_cox(d, univ_min, outcome):
    dx = cox_ready(d)
    dummy_map = {}
    for v in COX_CANDIDATES:
        if v in ("Phase_group", "Age_group"):
            dummy_map[v] = [c for c in dx.columns if c.startswith(v + "_")]
        else:
            dummy_map[v] = [v]

    sel = [v for v in univ_min.index if univ_min[v] < VAR_ENTRY_P and usable(d, v)]
    for fv in FORCE_VARS:
        if fv not in sel and usable(d, fv):
            sel.append(fv)
    if not sel:
        return None, None, []

    n_events = int(d["Event"].sum())
    max_cov = max(2, n_events // EPV_MIN)
    if len(sel) > max_cov:
        forced = [v for v in FORCE_VARS if v in sel]
        rest = [v for v in sel if v not in forced]
        rest.sort(key=lambda v: univ_min[v])
        sel = forced + rest[: max(0, max_cov - len(forced))]
        print(f"    [{outcome}] EPV guard: limited to {max_cov} covariates "
              f"(events={n_events}); kept: {sel}")

    terms = [t for v in sel for t in dummy_map[v]]
    sub_all = dx[["Time", "Event"] + terms].dropna()
    # drop zero-variance levels within the selected subset
    terms = [t for t in terms if sub_all[t].nunique() > 1]
    sub_all = sub_all[["Time", "Event"] + terms]

    # backward elimination by AIC (block-wise)
    cur_vars = list(sel)
    cur_terms = list(terms)
    fit = fit_block(sub_all, cur_terms)
    if fit is None:
        return None, None, sel
    best_aic = fit.AIC_partial_
    improved = True
    while improved and len(cur_vars) > 1:
        improved = False
        best_drop = None
        for v in cur_vars:
            if v in FORCE_VARS and len(cur_vars) > len(FORCE_VARS):
                continue  # keep forced vars unless nothing else left
            cand_vars = [x for x in cur_vars if x != v]
            cand_terms = [t for x in cand_vars for t in dummy_map[x]]
            f2 = fit_block(sub_all, cand_terms)
            if f2 is not None and f2.AIC_partial_ < best_aic - 1e-9:
                best_aic = f2.AIC_partial_
                best_drop = v
        if best_drop is not None:
            cur_vars = [x for x in cur_vars if x != best_drop]
            cur_terms = [t for x in cur_vars for t in dummy_map[x]]
            fit = fit_block(sub_all, cur_terms)
            improved = True
    # final refit without penalizer for reporting (separation-safe fallback)
    try:
        final = CoxPHFitter()
        final.fit(sub_all, "Time", "Event")
    except Exception:
        final = fit

    s = final.summary
    rows = []
    for term in cur_terms:
        r = s.loc[term]
        v = next(x for x in cur_vars if term in dummy_map[x])
        rows.append({
            "Outcome": outcome, "Covariate": v, "Label": VAR_LABELS[v],
            "Term": term,
            "HR": round(float(np.exp(r["coef"])), 2),
            "CI95_lower": round(float(np.exp(r["coef lower 95%"])), 2),
            "CI95_upper": round(float(np.exp(r["coef upper 95%"])), 2),
            "P_value": float(r["p"]),
            "P_report": fmt_p(float(r["p"])),
            "Report": (f"HR {np.exp(r['coef']):.2f} "
                       f"(95% CI {np.exp(r['coef lower 95%']):.2f} - "
                       f"{np.exp(r['coef upper 95%']):.2f}), "
                       f"p = {fmt_p(float(r['p']))}"),
        })
    multiv = pd.DataFrame(rows)

    # PH assumption check (Schoenfeld-type test, lifelines)
    zph = pd.DataFrame()
    try:
        from lifelines.statistics import proportional_hazard_test
        res = proportional_hazard_test(final, sub_all, time_transform="km")
        zph = res.summary.reset_index()
        zph = zph.rename(columns={"index": "Term", "p": "P_value"})
        zph["Outcome"] = outcome
        zph["Chi2"] = zph["test_statistic"].round(2)
        zph["Passes_PH"] = zph["P_value"] > 0.05
        zph = zph[["Outcome", "Term", "Chi2", "P_value", "Passes_PH"]]
        viol = zph[~zph["Passes_PH"]]
        if len(viol):
            print(f"    [{outcome}] PH violations (p<0.05): " +
                  ", ".join(f"{r.Term} (p={fmt_p(r.P_value)})" for r in viol.itertuples()))
    except Exception as e:
        print(f"    PH check failed for {outcome}: {e}")
    print(f"    [{outcome}] final multivariable model: n={sub_all.shape[0]}, "
          f"events={int(sub_all['Event'].sum())}, covariates={cur_vars}, "
          f"EPV={sub_all['Event'].sum() / max(1, len(cur_vars)):.1f}")
    return multiv, zph, sel


print("\n>>> Objective 3: Cox models (2-group design)...")
cox_frames, zph_frames = [], []
for d, outcome in [(os_grp, "OS"), (pfs_grp, "PFS")]:
    univ, univ_min = univariable_cox(d, outcome)
    if len(univ):
        cox_frames.append(univ)
    multiv, zph, sel = multivariable_cox(d, univ_min, outcome)
    if multiv is not None and len(multiv):
        cox_frames.append(multiv.assign(Model="multivariable"))
        univ2 = univ.assign(Model="univariable") if len(univ) else univ
        print(f"\n    [{outcome}] multivariable covariates: {sel}")
        print(multiv[["Covariate", "Term", "HR", "CI95_lower", "CI95_upper",
                      "P_report"]].to_string(index=False))
    if zph is not None and len(zph):
        zph_frames.append(zph)

def sep_note(row):
    if row["CI95_upper"] > 1e5 or (row["HR"] == 0.0):
        return ("quasi-complete separation - HR unstable, interpret with caution")
    return ""

for i in range(len(cox_frames)):
    if len(cox_frames[i]):
        cox_frames[i] = cox_frames[i].assign(
            Note=cox_frames[i].apply(sep_note, axis=1))

univ_all = pd.concat([f for f in cox_frames if "Model" not in f.columns],
                     ignore_index=True) if cox_frames else pd.DataFrame()
mult_all = pd.concat([f for f in cox_frames if "Model" in f.columns],
                     ignore_index=True) if cox_frames else pd.DataFrame()
if len(univ_all):
    save_csv(univ_all, "O3_01_Univariate_Cox.csv")
if len(mult_all):
    save_csv(mult_all, "O3_02_Multivariate_Cox.csv")
if zph_frames:
    save_csv(pd.concat(zph_frames, ignore_index=True), "O3_03_PH_Assumption_Checks.csv")

# =============================================================================
# 6. OBJECTIVE 4 - ADRs
# =============================================================================
print("\n>>> Objective 4: ADR prevalence...")
n_adr_known = int(df["flag_ADR_any"].notna().sum())
n_adr = int(df["flag_ADR_any"].fillna(False).sum())
n_sev_known = int(df["flag_ADR_sev"].notna().sum())
n_sev = int(df["flag_ADR_sev"].fillna(False).sum())
ci_lo, ci_hi = prop_ci(n_adr, n_adr_known)

o4_01 = pd.DataFrame([
    {"Measure": "Any ADR (documented)", "N_documented": n_adr_known, "n_ADR": n_adr,
     "Pct_of_documented": round(100 * n_adr / n_adr_known, 1),
     "CI95_pct": f"{100 * ci_lo:.1f} - {100 * ci_hi:.1f}"},
    {"Measure": "Grade 3+ (severe) ADR (documented)", "N_documented": n_sev_known,
     "n_ADR": n_sev, "Pct_of_documented": round(100 * n_sev / n_sev_known, 1),
     "CI95_pct": " - ".join(f"{100 * x:.1f}" for x in prop_ci(n_sev, n_sev_known))},
    {"Measure": "No ADR documentation", "N_documented": len(df) - n_adr_known,
     "n_ADR": np.nan, "Pct_of_documented": np.nan, "CI95_pct": ""},
])
save_csv(o4_01, "O4_01_ADR_Prevalence_Overall.csv")

# ADR by line group (2 groups)
o4_03 = d_grp.groupby("LINE").apply(
    lambda g: pd.Series({
        "n": len(g),
        "n_ADR_known": int(g["flag_ADR_any"].notna().sum()),
        "n_ADR": int(g["flag_ADR_any"].fillna(False).sum()),
        "Prevalence_ADR_pct": round(100 * g["flag_ADR_any"].fillna(False).sum() /
                                    max(1, g["flag_ADR_any"].notna().sum()), 1),
        "n_Severe": int(g["flag_ADR_sev"].fillna(False).sum()),
        "Prevalence_Severe_pct": round(100 * g["flag_ADR_sev"].fillna(False).sum() /
                                       max(1, g["flag_ADR_sev"].notna().sum()), 1),
    }), include_groups=False).reset_index()
tab_adr = pd.crosstab(d_grp["LINE"], d_grp["flag_ADR_any"])
o4_03["p_value_line_comparison"] = cat_pvalue(tab_adr.values)
o4_03["p_value_report"] = o4_03["p_value_line_comparison"].map(fmt_p)
save_csv(o4_03, "O4_03_ADR_Prevalence_by_Line_Group.csv")

# ADR by individual TKI (descriptive only - panel comment ix)
o4_02 = df.dropna(subset=["tki_clean"]).groupby("tki_clean").apply(
    lambda g: pd.Series({
        "n": len(g),
        "n_ADR_known": int(g["flag_ADR_any"].notna().sum()),
        "n_ADR": int(g["flag_ADR_any"].fillna(False).sum()),
        "Prevalence_ADR_pct": round(100 * g["flag_ADR_any"].fillna(False).sum() /
                                    max(1, g["flag_ADR_any"].notna().sum()), 1),
        "n_Severe": int(g["flag_ADR_sev"].fillna(False).sum()),
    }), include_groups=False).reset_index().sort_values("n", ascending=False)
save_csv(o4_02, "O4_02_ADR_Prevalence_by_TKI_descriptive.csv")

# grade 3+ ADR types
g3txt = df.loc[df["flag_ADR_sev"] == True, "Grade_3_ADR_Specified"].dropna()  # noqa E712
g3_long = (g3txt.astype(str).str.split(r"[,;]").explode().str.strip()
           .replace("", np.nan).dropna())
if len(g3_long):
    o4_04 = (g3_long.value_counts().rename_axis("ADR_type").reset_index(name="n"))
    n_g3_pat = int((df["flag_ADR_sev"] == True).sum())  # noqa E712
    o4_04["Pct_of_grade3_patients"] = (100 * o4_04["n"] / n_g3_pat).round(1)
    save_csv(o4_04, "O4_04_Grade3_Specified_ADRs.csv")

# =============================================================================
# 7. SUMMARY REPORT
# =============================================================================
def get_km(tbl, outcome, group, field):
    r = tbl[(tbl["Outcome"] == outcome) & (tbl["Group"] == group)]
    if len(r) == 0:
        return "NA"
    v = r.iloc[0][field]
    return v if not pd.isna(v) else "NA"


lr_os = lr_tbl.loc[lr_tbl["Outcome"] == "OS", "P_value"].iloc[0]
lr_pfs = lr_tbl.loc[lr_tbl["Outcome"] == "PFS", "P_value"].iloc[0]

md = []
md.append("# CML TKI cohort (KNH) - regrouped analysis, trial 5\n")
md.append("**New grouping (per supervisor feedback):** 2nd-generation TKIs are analysed "
          "together, and the primary comparison is **two groups only** - "
          "**first line vs subsequent line** - because only 5 patients were recorded "
          "on ponatinib, which is too few for any meaningful separate analysis.\n")
md.append("| Group | Definition (recorded TKI) | n |")
md.append("|---|---|---|")
md.append(f"| First line | imatinib | {int((df['LINE'] == LINE_LEVELS[0]).sum())} |")
md.append(f"| Subsequent line | dasatinib / nilotinib / bosutinib / ponatinib "
          f"(2G TKIs grouped + ponatinib) | {int((df['LINE'] == LINE_LEVELS[1]).sum())} |")
md.append(f"\nClassifiable patients: {n_class} of {len(df)} "
          f"(1 patient with blank TKI excluded from group comparisons; retained in overall estimates).\n")

md.append("## Survival by line group\n")
md.append("| Outcome | Group | n | Events | Median (mo) | 1-yr % | 3-yr % | 5-yr % (95% CI) | "
          "RMST to 60 mo | Log-rank p |")
md.append("|---|---|---|---|---|---|---|---|---|")
for outcome, tbl, lrp in [("OS", km_os, lr_os), ("PFS", km_pfs, lr_pfs)]:
    for g in ["Overall"] + LINE_LEVELS:
        r = tbl[tbl["Group"] == g].iloc[0]
        med = r["Median_months"] if not pd.isna(r["Median_months"]) else "not reached"
        md.append(f"| {outcome} | {g} | {int(r['N'])} | {int(r['Events_n'])} "
                  f"({r['Events_pct']}%) | {med} | {r['S_12mo']} | {r['S_36mo']} "
                  f"| {r['S_60mo']} ({r['S_60mo_CI95']}) | {r['RMST_60mo_months']} "
                  f"| {fmt_p(lrp)} |")
md.append("")
md.append("PH-free between-group difference (restricted mean survival time to 60 "
          "months, bootstrap 95% CI, B=1000):")
for _, r in rmst_df.iterrows():
    md.append(f"- {r['Outcome']}: first line minus subsequent line = "
              f"{r['Difference_mo']} months (95% CI {r['CI95_boot']}), "
              f"p = {r['P_report']}")
md.append("")

if len(mult_all):
    md.append("## Multivariable Cox regression (adjusted)\n")
    md.append("| Outcome | Covariate | Term | HR | 95% CI | p |")
    md.append("|---|---|---|---|---|---|")
    for _, r in mult_all.iterrows():
        note = r.get("Note", "")
        md.append(f"| {r['Outcome']} | {r['Label']} | {r['Term']} | {r['HR']} "
                  f"| {r['CI95_lower']} - {r['CI95_upper']} | {r['P_report']}"
                  f"{'; ' + note if note else ''} |")
    md.append("")

if zph_frames:
    md.append("### Proportional-hazards assumption (Schoenfeld-type test)\n")
    md.append("| Outcome | Term | p | PH OK? |")
    md.append("|---|---|---|---|")
    zph_all = pd.concat(zph_frames, ignore_index=True)
    for _, r in zph_all.iterrows():
        md.append(f"| {r['Outcome']} | {r['Term']} | {fmt_p(r['P_value'])} "
                  f"| {'yes' if r['Passes_PH'] else 'NO - see note'} |")
    md.append("")
    md.append("Where the line-group term violates proportionality, the single hazard "
              "ratio is an *average* effect; the log-rank test and the RMST difference "
              "above are the PH-free primary comparisons.")
    md.append("")

md.append("## Key points\n")
md.append(f"- First line (imatinib) n={int((df['LINE']==LINE_LEVELS[0]).sum())} vs "
          f"subsequent line n={int((df['LINE']==LINE_LEVELS[1]).sum())}: both groups are "
          "now large enough for meaningful comparison.")
md.append("- Ponatinib patients (n=5) are pooled into the subsequent-line group; "
          "their individual details remain in O1_10_Ponatinib_Patients_n5.csv.")
md.append("- Line of therapy is INFERRED from the recorded current TKI (2G/3G patients "
          "are presumed imatinib-exposed). Worse subsequent-line outcomes reflect "
          "confounding by indication, not drug inferiority.")
md.append("- 2nd-generation TKIs (dasatinib/nilotinib/bosutinib) are never split in the "
          "comparative analyses; per-drug tables are kept as descriptive only.")

with open(os.path.join(OUT, "ANALYSIS_SUMMARY_2GROUP.md"), "w") as f:
    f.write("\n".join(md))
saved.append("ANALYSIS_SUMMARY_2GROUP.md")

print("\n>>> ALL SAVED FILES:")
for s in saved:
    print("   ", s)
print("\nDONE.")
