#!/usr/bin/env python3
# =============================================================================
#  OBJECTIVE 3 - PREDICTORS OF SURVIVAL (OS and PFS), trial 5.xlsx
#  CML patients managed with TKIs at Kenyatta National Hospital
#
#  EXACT SPECIFICATION (supervisor, Sep 2026):
#   * LOT: First line = imatinib;
#          Subsequent line = dasatinib/nilotinib/bosutinib/ponatinib.
#   * Grade >=3 ADR: Yes/No binary.
#   * Check OS event coding; check PFS event coding (death MUST be a PFS event).
#   * LOT descriptive tables; Grade >=3 ADR descriptive tables.
#   * KM: OS by LOT; PFS by LOT; OS by Grade>=3 ADR; PFS by Grade>=3 ADR.
#   * Univariable Cox: LOT->OS, LOT->PFS, Grade>=3 ADR->OS, Grade>=3 ADR->PFS.
#   * Multivariable Cox (OS and PFS) with the FIXED covariate set:
#         CML phase + gender + age + LOT + Grade>=3 ADR  (no stepwise).
#   * Proportional-hazards assumption tests for the models.
#
#  Analysis cohort: the 162 patients with a classifiable LOT (1 patient with a
#  blank TKI field is excluded). Grade>=3 ADR coding: 'Yes' -> Yes;
#  'No'/'Not documented'/'Not applicable' -> No (assumption stated in report;
#  31 of the 47 'Not documented' patients had no ADR at all).
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
import matplotlib.gridspec as gridspec

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import multivariate_logrank_test, proportional_hazard_test

warnings.filterwarnings("ignore")

OUT = "output"
os.makedirs(OUT, exist_ok=True)
saved = []

FIRST = "First line (imatinib)"
SUBSQ = "Subsequent line (2G/3G TKI)"


def save_csv(df, name):
    df.to_csv(os.path.join(OUT, name), index=False)
    saved.append(name)
    print(">>> saved table:  ", name)


def save_fig(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=300, bbox_inches="tight")
    plt.close(fig)
    saved.append(name)
    print(">>> saved figure: ", name)


def fmt_p(p):
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "NA"
    return "<0.001" if p < 0.001 else f"{p:.3f}"


def median_iqr(x):
    x = pd.Series(x).dropna()
    if len(x) == 0:
        return "NA"
    q1, q3 = np.percentile(x, [25, 75])
    return f"{np.median(x):.1f} ({q1:.1f} - {q3:.1f})"


def cat_pvalue(tab):
    tab = np.asarray(tab)
    tab = tab[tab.sum(axis=1) > 0][:, tab.sum(axis=0) > 0]
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return np.nan
    if tab.shape == (2, 2):
        return stats.fisher_exact(tab)[1]
    try:
        return stats.chi2_contingency(tab, correction=False)[1]
    except ValueError:
        return np.nan


# =============================================================================
# 1. LOAD + VARIABLES
# =============================================================================
print(">>> Loading trial 5.xlsx")
raw = pd.read_excel("trial 5.xlsx", sheet_name="Sheet2")
raw.columns = [str(c).strip() for c in raw.columns]
df = raw.drop(columns=["Event"], errors="ignore")  # raw 'Event' duplicates OS_Event

# --- LOT (spec) ---
tki = df["TKI Used"].astype(str).str.strip().str.lower()
tki = tki.replace({"nan": np.nan, "none": np.nan, "": np.nan})
df["LOT"] = pd.Series(np.nan, index=df.index, dtype=object)
df.loc[tki == "imatinib", "LOT"] = FIRST
df.loc[tki.isin(["dasatinib", "nilotinib", "bosutinib", "ponatinib"]), "LOT"] = SUBSQ

# --- Grade >=3 ADR: Yes/No (spec) ---
g3raw = df["Grade_3_ADR"].astype(str).str.strip()
df["Grade3_ADR_raw"] = g3raw
df["G3ADR"] = (g3raw.str.lower() == "yes")

# --- endpoints ---
df["OS_Time"] = pd.to_numeric(df["OS_Time_Months"], errors="coerce")
df["OS_Event"] = pd.to_numeric(df["OS_Event"], errors="coerce")
df["PFS_Time"] = pd.to_numeric(df["PFS_Time_Months"], errors="coerce")
df["PFS_Event"] = pd.to_numeric(df["PFS_Event"], errors="coerce")

# --- covariates for the adjusted models (spec) ---
ph = df["Phase of CML"].astype(str).str.strip().str.lower()
df["Phase_group"] = pd.Series(np.nan, index=df.index, dtype=object)
df.loc[ph.str.contains("chron", na=False), "Phase_group"] = "Chronic"
df.loc[ph.str.contains("accel", na=False), "Phase_group"] = "Accelerated"
df.loc[ph.str.contains("blast", na=False), "Phase_group"] = "Blast"
df["Phase_group"] = pd.Categorical(df["Phase_group"],
                                   categories=["Chronic", "Accelerated", "Blast"])
df["Sex"] = df["Gender"].astype(str).str.strip()
df.loc[~df["Sex"].isin(["Male", "Female"]), "Sex"] = np.nan
df["Age_num"] = pd.to_numeric(df["Age"], errors="coerce")

# =============================================================================
# 2. EVENT-CODING CHECKS
# =============================================================================
print("\n>>> Event-coding checks")
ss = df["Survival_Status"].astype(str).str.strip()
death_txt = ss == "Deceased"
alive_txt = ss.isin(["Alive", "Alive/Censored"])
dod = df["Date_of_Death"].notna()
rel_yes = df["Disease_Relapse"].astype(str).str.strip().str.lower() == "yes"

chk = pd.DataFrame([
    ["OS", "Patients", len(df), ""],
    ["OS", "OS_Event=1 (death)", int((df.OS_Event == 1).sum()), ""],
    ["OS", "Survival_Status='Deceased'", int(death_txt.sum()), ""],
    ["OS", "Deaths with a recorded death date", int((death_txt & dod).sum()), ""],
    ["OS", "Mismatch OS_Event vs Survival_Status",
     int(((df.OS_Event == 1) != death_txt).sum()), "must be 0"],
    ["OS", "OS_Event=0 (censored) matching Alive/Alive-Censored",
     int(((df.OS_Event == 0) & alive_txt).sum()), ""],
    ["OS", "OS time <=0 or missing", int((df.OS_Time.isna() | (df.OS_Time <= 0)).sum()),
     "must be 0"],
    ["PFS", "Patients", len(df), ""],
    ["PFS", "PFS_Event=1 (event)", int((df.PFS_Event == 1).sum()), ""],
    ["PFS", "  of which deaths (OS_Event=1)",
     int(((df.PFS_Event == 1) & (df.OS_Event == 1)).sum()), ""],
    ["PFS", "  of which non-death events (alive, progression/relapse)",
     int(((df.PFS_Event == 1) & (df.OS_Event == 0)).sum()), ""],
    ["PFS", "DEATHS NOT CODED AS PFS EVENT",
     int(((df.OS_Event == 1) & (df.PFS_Event != 1)).sum()), "must be 0"],
    ["PFS", "Non-death PFS events with documented relapse",
     int(((df.PFS_Event == 1) & (df.OS_Event == 0) & rel_yes).sum()), ""],
    ["PFS", "PFS time <=0 or missing",
     int((df.PFS_Time.isna() | (df.PFS_Time <= 0)).sum()), "must be 0"],
], columns=["Endpoint", "Check", "n", "Requirement"])
print(chk.to_string(index=False))
save_csv(chk, "O3_00_Event_Coding_Checks.csv")

# =============================================================================
# 3. ANALYSIS COHORT (classifiable LOT)
# =============================================================================
d = df[df["LOT"].notna()].copy()
print(f"\n>>> Analysis cohort: {len(d)} patients "
      f"(excluded {len(df) - len(d)} with blank TKI)")

# =============================================================================
# 4. DESCRIPTIVE TABLES: LOT
# =============================================================================
print("\n>>> LOT descriptive tables")
lot = (d["LOT"].value_counts().rename_axis("LOT").reset_index(name="n"))
lot["Pct"] = (100 * lot["n"] / len(d)).round(1)
lot = lot.set_index("LOT").loc[[FIRST, SUBSQ]].reset_index()
lot_break = (d[d["LOT"] == SUBSQ]["TKI Used"].astype(str).str.strip()
             .value_counts().rename_axis("recorded_TKI").reset_index(name="n"))
lot_break["Pct_of_subsequent"] = (100 * lot_break["n"] / lot_break["n"].sum()).round(1)

rows = [{"LOT table": "Overall", "First line": lot.n[0], "Subsequent line": lot.n[1]}]
desc_frames = [lot.assign(table="LOT_distribution")]
desc_frames.append(lot_break.assign(table="Subsequent_line_breakdown"))

lot_x = []
for var, lab in [("Phase_group", "Phase"), ("Sex", "Sex")]:
    ct = pd.crosstab(d["LOT"], d[var])
    p = cat_pvalue(ct.values)
    for level in ct.columns:
        lot_x.append({"Variable": lab, "Level": level,
                      FIRST: int(ct.loc[FIRST, level]),
                      SUBSQ: int(ct.loc[SUBSQ, level]),
                      "P_value": p})
lot_x.append({"Variable": "Age (years), median (IQR)", "Level": "",
              FIRST: median_iqr(d.loc[d.LOT == FIRST, "Age_num"]),
              SUBSQ: median_iqr(d.loc[d.LOT == SUBSQ, "Age_num"]),
              "P_value": stats.mannwhitneyu(
                  d.loc[d.LOT == FIRST, "Age_num"].dropna(),
                  d.loc[d.LOT == SUBSQ, "Age_num"].dropna()).pvalue})
lot_desc = pd.concat(desc_frames, ignore_index=True)
lot_cross = pd.DataFrame(lot_x)
save_csv(lot, "O3_01a_LOT_Distribution.csv")
save_csv(lot_break, "O3_01b_LOT_Subsequent_Line_Breakdown.csv")
save_csv(lot_cross, "O3_01c_LOT_by_Baseline_Factors.csv")
print(lot.to_string(index=False)); print(lot_break.to_string(index=False))

# =============================================================================
# 5. DESCRIPTIVE TABLES: GRADE >=3 ADR
# =============================================================================
print("\n>>> Grade>=3 ADR descriptive tables")
g3_raw_tab = (d["Grade3_ADR_raw"].value_counts().rename_axis("Grade_3_ADR_recorded")
              .reset_index(name="n"))
g3_raw_tab["Pct"] = (100 * g3_raw_tab["n"] / len(d)).round(1)
g3_bin = (d["G3ADR"].map({True: "Yes", False: "No"})
          .value_counts().rename_axis("Grade3_ADR_binary").reset_index(name="n"))
g3_bin["Pct"] = (100 * g3_bin["n"] / len(d)).round(1)
g3_by_lot = pd.crosstab(d["LOT"], d["G3ADR"].map({True: "Yes", False: "No"}))
g3_by_lot["Pct_Yes"] = (100 * g3_by_lot["Yes"] / g3_by_lot.sum(axis=1)).round(1)
g3_by_lot["P_value_fisher"] = cat_pvalue(g3_by_lot[["No", "Yes"]].values)
g3_types = (d.loc[d["G3ADR"], "Grade_3_ADR_Specified"].dropna().astype(str)
            .str.split(r"[,;]").explode().str.strip()
            .replace("", np.nan).dropna().value_counts()
            .rename_axis("ADR_type").reset_index(name="n"))
g3_types["Pct_of_G3_patients"] = (100 * g3_types["n"] / d["G3ADR"].sum()).round(1)

save_csv(g3_raw_tab, "O3_02a_Grade3ADR_Recorded_Categories.csv")
save_csv(g3_bin, "O3_02b_Grade3ADR_Binary.csv")
save_csv(g3_by_lot.reset_index(), "O3_02c_Grade3ADR_by_LOT.csv")
save_csv(g3_types, "O3_02d_Grade3ADR_Types.csv")
print(g3_raw_tab.to_string(index=False)); print(g3_bin.to_string(index=False))
print(g3_by_lot.to_string())

# =============================================================================
# 6. KAPLAN-MEIER CURVES (by LOT and by Grade>=3 ADR)
# =============================================================================
print("\n>>> Kaplan-Meier curves")


def km_figure(d, var, levels, colors, ylabel, fname, title):
    fig = plt.figure(figsize=(9, 6.6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.06)
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax)
    sub = d[d[var].notna()]
    lab_of = {FIRST: "First line", SUBSQ: "Subsequent",
              True: "Grade>=3 ADR", False: "No grade>=3 ADR"}
    for g, c in zip(levels, colors):
        dg = sub[sub[var] == g]
        kmf = KaplanMeierFitter().fit(dg["Time"], dg["Event"], label=lab_of[g])
        kmf.plot_survival_function(ax=ax, ci_show=True, color=c, linewidth=2)
    lr = multivariate_logrank_test(sub["Time"], sub[var], sub["Event"])
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="lower left", frameon=False)
    ax.text(0.98, 0.97, f"Log-rank p {fmt_p(lr.p_value)}", transform=ax.transAxes,
            fontsize=11, ha="right", va="top",
            bbox=dict(boxstyle="round", fc="white", ec="0.6", alpha=0.9))
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.25)
    plt.setp(ax.get_xticklabels(), visible=False)
    tmax = min(181, sub["Time"].max())
    times = np.arange(0, tmax + 1, 24)
    for i, (g, c) in enumerate(zip(levels, colors)):
        dg = sub[sub[var] == g]
        risk = [(dg["Time"] >= t).sum() for t in times]
        y = 1 - 0.45 * i
        ax2.text(-0.02, y, lab_of[g], transform=ax2.get_yaxis_transform(),
                 ha="right", va="center", fontsize=9, color=c)
        for t, r in zip(times, risk):
            ax2.text(t, y, str(r), ha="center", va="center", fontsize=8, color=c)
    ax2.set_yticks([])
    ax2.set_ylim(0, 1)
    ax2.set_xlim(0, tmax)
    ax2.set_xticks(times)
    ax2.set_xlabel("Follow-up (months)")
    ax2.set_title("Number at risk", fontsize=9, loc="left")
    for s in ["top", "right", "left"]:
        ax2.spines[s].set_visible(False)
    save_fig(fig, fname)
    return lr.p_value


def endpoint(d, time_col, event_col):
    x = d.rename(columns={time_col: "Time", event_col: "Event"})
    return x


os_d = endpoint(d, "OS_Time", "OS_Event")
pfs_d = endpoint(d, "PFS_Time", "PFS_Event")
g3map = {True: True, False: False}

p_os_lot = km_figure(os_d, "LOT", [FIRST, SUBSQ], ["#1f77b4", "#d62728"],
                     "Probability of survival", "O3_03_KM_OS_by_LOT.png",
                     "Overall survival by line of therapy\n(first line vs subsequent line)")
p_pfs_lot = km_figure(pfs_d, "LOT", [FIRST, SUBSQ], ["#1f77b4", "#d62728"],
                      "Probability of being progression-free",
                      "O3_04_KM_PFS_by_LOT.png",
                      "Progression-free survival by line of therapy\n(first line vs subsequent line)")
p_os_g3 = km_figure(os_d, "G3ADR", [False, True], ["#2ca02c", "#9467bd"],
                    "Probability of survival", "O3_05_KM_OS_by_Grade3ADR.png",
                    "Overall survival by Grade >=3 ADR")
p_pfs_g3 = km_figure(pfs_d, "G3ADR", [False, True], ["#2ca02c", "#9467bd"],
                     "Probability of being progression-free",
                     "O3_06_KM_PFS_by_Grade3ADR.png",
                     "Progression-free survival by Grade >=3 ADR")

km_summary = pd.DataFrame([
    ["OS", "LOT", fmt_p(p_os_lot)], ["PFS", "LOT", fmt_p(p_pfs_lot)],
    ["OS", "Grade>=3 ADR", fmt_p(p_os_g3)], ["PFS", "Grade>=3 ADR", fmt_p(p_pfs_g3)]],
    columns=["Outcome", "Grouping", "Log_rank_p"])
save_csv(km_summary, "O3_07_KM_LogRank_Summary.csv")

# =============================================================================
# 7. COX REGRESSION
# =============================================================================
print("\n>>> Cox regression")


def cox_design(d):
    """Numeric design matrix for the specified covariates."""
    x = d.copy()
    x["LOT_subsequent"] = (x["LOT"] == SUBSQ).astype(float)
    x["G3ADR_yes"] = x["G3ADR"].astype(float)
    x["Male"] = (x["Sex"] == "Male").astype(float).where(x["Sex"].notna(), np.nan)
    x["Phase_Accelerated"] = (x["Phase_group"] == "Accelerated").astype(float)
    x["Phase_Blast"] = (x["Phase_group"] == "Blast").astype(float)
    return x


def ph_table(fit, sub, outcome, model_name):
    res = proportional_hazard_test(fit, sub, time_transform="km")
    z = res.summary.reset_index().rename(columns={"index": "Term", "p": "P_value"})
    z.insert(0, "Outcome", outcome)
    z.insert(1, "Model", model_name)
    z["Chi2"] = z["test_statistic"].round(2)
    z["PH_ok"] = z["P_value"] > 0.05
    return z[["Outcome", "Model", "Term", "Chi2", "P_value", "PH_ok"]]


def tidy(fit, outcome, model_name):
    s = fit.summary
    rows = []
    for term in s.index:
        r = s.loc[term]
        hr = np.exp(r["coef"])
        lo = np.exp(r["coef lower 95%"])
        hi = np.exp(r["coef upper 95%"])
        rows.append({"Outcome": outcome, "Model": model_name, "Term": term,
                     "HR": round(hr, 2), "CI95_lower": round(lo, 2),
                     "CI95_upper": round(hi, 2), "P_value": round(float(r["p"]), 4),
                     "P_report": fmt_p(float(r["p"])),
                     "Report": f"HR {hr:.2f} (95% CI {lo:.2f} - {hi:.2f}), "
                               f"p = {fmt_p(float(r['p']))}"})
    return pd.DataFrame(rows)


univ_frames, ph_frames = [], []
LABEL = {"LOT_subsequent": "LOT (subsequent vs first line)",
         "G3ADR_yes": "Grade>=3 ADR (Yes vs No)"}
for d_ep, oc, time_col, ev_col in [(d, "OS", "OS_Time", "OS_Event"),
                                   (d, "PFS", "PFS_Time", "PFS_Event")]:
    x = cox_design(d_ep)
    for var in ["LOT_subsequent", "G3ADR_yes"]:
        sub = x[[time_col, ev_col, var]].dropna()
        sub = sub.rename(columns={time_col: "Time", ev_col: "Event"})
        fit = CoxPHFitter(penalizer=1e-4).fit(sub, "Time", "Event")
        t = tidy(fit, oc, f"Univariable: {LABEL[var]}")
        t["Covariate"] = LABEL[var]
        univ_frames.append(t)
        ph_frames.append(ph_table(fit, sub, oc, f"Univariable {LABEL[var]}"))

univ = pd.concat(univ_frames, ignore_index=True)
save_csv(univ, "O3_08_Univariate_Cox_LOT_Grade3ADR.csv")
print(univ[["Outcome", "Covariate", "HR", "CI95_lower", "CI95_upper",
            "P_report"]].to_string(index=False))

# --- multivariable: FIXED set (phase + gender + age + LOT + Grade>=3 ADR) ---
mult_frames = []
for d_ep, oc, time_col, ev_col in [(d, "OS", "OS_Time", "OS_Event"),
                                   (d, "PFS", "PFS_Time", "PFS_Event")]:
    x = cox_design(d_ep)
    covs = ["Phase_Accelerated", "Phase_Blast", "Male", "Age_num",
            "LOT_subsequent", "G3ADR_yes"]
    sub = x[[time_col, ev_col] + covs].dropna()
    sub = sub.rename(columns={time_col: "Time", ev_col: "Event"})
    n_events = int(sub["Event"].sum())
    try:
        fit = CoxPHFitter().fit(sub, "Time", "Event")
    except Exception:
        fit = CoxPHFitter(penalizer=1e-3).fit(sub, "Time", "Event")
    t = tidy(fit, oc, "Multivariable (phase+sex+age+LOT+G3ADR)")
    # add per-10-years HR for age
    age = t[t["Term"] == "Age_num"]
    if len(age):
        hr10 = age.iloc[0]["HR"] ** 10
        t.loc[t["Term"] == "Age_num", "Note"] = (
            f"per 10 yrs: HR {hr10:.2f}")
    t["events"] = n_events
    t["n"] = len(sub)
    t["EPV"] = round(n_events / 6, 1)
    mult_frames.append(t)
    ph_frames.append(ph_table(fit, sub, oc,
                              "Multivariable phase+sex+age+LOT+G3ADR"))
    print(f"\n    [{oc}] multivariable: n={len(sub)}, events={n_events}, "
          f"EPV={n_events / 6:.1f}")
    print(t[["Term", "HR", "CI95_lower", "CI95_upper", "P_report"]].to_string(
        index=False))

mult = pd.concat(mult_frames, ignore_index=True)
save_csv(mult, "O3_09_Multivariate_Cox_OS_PFS.csv")

ph = pd.concat(ph_frames, ignore_index=True)
save_csv(ph, "O3_10_Proportional_Hazards_Tests.csv")
print("\n>>> Proportional-hazards tests (violations, p<0.05):")
viol = ph[~ph["PH_ok"]]
print(viol.to_string(index=False) if len(viol) else "    none")

# =============================================================================
# 8. SUMMARY MD
# =============================================================================
md = ["# Objective 3 - predictors of survival (OS & PFS), trial 5\n",
      "Analysis cohort: n=%d (1 patient with blank TKI excluded). " % len(d) +
      "Grade>=3 ADR coded Yes/No ('Not documented'/'Not applicable' -> No; "
      "31 of those 47 had no ADR recorded at all).\n",
      "## Event-coding checks\n",
      "- OS: %d deaths = %d OS_Event=1; 0 mismatches vs Survival_Status; "
      % (int((d.OS_Event == 1).sum()), int((d.OS_Event == 1).sum())) +
      "every death has a recorded death date; 0 non-positive OS times.",
      "- PFS: %d events = %d deaths + %d non-death progression/relapse events; "
      % (int((d.PFS_Event == 1).sum()),
         int(((d.PFS_Event == 1) & (d.OS_Event == 1)).sum()),
         int(((d.PFS_Event == 1) & (d.OS_Event == 0)).sum())) +
      "**0 deaths missed as PFS events** (death is treated as a PFS event).",
      "- Full check table: O3_00_Event_Coding_Checks.csv.\n",
      "## Descriptive\n",
      f"- LOT: first line n={int(lot.n[0])} ({lot.Pct[0]}%), subsequent line "
      f"n={int(lot.n[1])} ({lot.Pct[1]}%).",
      f"- Grade>=3 ADR: Yes n={int(d.G3ADR.sum())} "
      f"({100 * d.G3ADR.mean():.1f}%), No n={int((~d.G3ADR).sum())}.\n",
      "## Univariable Cox\n",
      "| Outcome | Covariate | HR | 95% CI | p |", "|---|---|---|---|---|"]
for _, r in univ.iterrows():
    md.append(f"| {r.Outcome} | {r.Covariate} | {r.HR} | {r.CI95_lower} - "
              f"{r.CI95_upper} | {r.P_report} |")
md += ["\n## Multivariable Cox (fixed set: CML phase + gender + age + LOT + "
       "Grade>=3 ADR)\n",
       "| Outcome | Term | HR | 95% CI | p | Note |", "|---|---|---|---|---|---|"]
for _, r in mult.iterrows():
    note = r.get("Note", "")
    if not isinstance(note, str) or pd.isna(note):
        note = ""
    md.append(f"| {r.Outcome} | {r.Term} | {r.HR} | {r.CI95_lower} - "
              f"{r.CI95_upper} | {r.P_report} | {note} |")
md += ["\n## Proportional hazards\n"]
if len(viol):
    md += [f"- Violation (p<0.05): {r.Model} / {r.Term}, p={fmt_p(r.P_value)}"
           for r in viol.itertuples()]
    md.append("- Where LOT or another term violates PH, the HR is an average "
              "effect; the log-rank KM comparison is the PH-free counterpart.")
else:
    md.append("- No PH violations at the 5% level.")
md.append("")

with open(os.path.join(OUT, "OBJECTIVE3_SUMMARY.md"), "w") as f:
    f.write("\n".join(md))
saved.append("OBJECTIVE3_SUMMARY.md")

print("\n>>> ALL SAVED:")
for s in saved:
    print("   ", s)
print("DONE.")
