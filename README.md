# Volatility Spillovers Between Housing and the Stock Market in Iran
**Replication package (EViews 13+)**

This repository reproduces the paper’s baseline:
- **ARMA–GARCH(1,1)** margins for monthly returns of **Tehran housing** (CBI average price per sqm) and **TSE all-share index**.
- **Two-step DCC(1,1)** (Engle, 2002) on **standardized residuals**.
- **Cheung–Ng** variance-causality on squared standardized residuals with **m = 12** (monthly one-year horizon).

All scripts are written for **EViews 13 or later**.

---

## 1. Data
- Input file: `mgarch_gregorian_unmodified.xlsx` (place **next to** the PRG).
- Sheet: `gregorian`
- Columns:
  - `date_gregorian` (YYYY-MM, range: 2016-04 to 2023-10; month-end)
  - `stock_index_close` (TSE all-share, end-of-month)
  - `tehran_price_per_m2_rial` (CBI average transaction price/sqm, Tehran, end-of-month)
- The PRG computes monthly log returns as `100 * Δln(P_t)` and aligns the sample to the non-missing overlap.

**Data note.** The housing return series exhibits extreme tails (min ~ −229%, max ~ +240%), consistent with abrupt market adjustments and retained outliers. Results are robust to these features.

---

## 2. Software and versions
- **EViews** 13 (build with MGARCH DCC support)
- OS: any
- No external add-ins are required.

---

## 3. How to run
1. Put `MGARCH_IR_main.prg` (this PRG) and `mgarch_gregorian_unmodified.xlsx` in the **same folder**.
2. Open the PRG in EViews and run.  
   The script will:
   - create a monthly workfile (2016M04–2023M10),
   - import data,
   - compute returns and align samples,
   - select **ARMA(p,q)** by AIC over `{0,1,2}`,
   - fit **GARCH(1,1)** to ARMA residuals (Gaussian QML),
   - standardize residuals and estimate **DCC(1,1)**,
   - compute **Cheung–Ng** portmanteau statistics (m = 12, plus robustness m = 6,18),
   - produce tables/figures and export them.

Outputs:
- Workfile: `MGARCH_IR_EViews.wf1`
- Tables (Excel): `results_tables_from_eviews.xlsx`
- Figures (PNG): `Figure1_conditional_vols.png`, `Figure2_dcc_correlation.png`
- Spool: `SP_ALL` (contains all tables, diagnostics, and figures)

---

## 4. Robustness
We assess sensitivity along three dimensions (as in the paper):
1. **CCC vs DCC.** Likelihood ratio test of `H0: a = b = 0` (CCC) vs DCC; see `tbl_ccc_dcc` and `tbl_ccc_dcc_summary`.
2. **Causality horizon.** Cheung–Ng with `m ∈ {6, 12, 18}`; see `tbl_cn_12`, `tbl_cn_multi`, and `tbl_cn_summary`.
3. **Heavy tails.** Re-estimate GARCH margins under **Student-t QML** (objects `g_stock_t`, `g_house_t`), rebuild standardized residuals, re-estimate DCC (`dcc_eq_t`), and recompute causality; see `tbl_dcc_t` and `tbl_cn_summary_t`.  
   *Note:* If your EViews build uses a different keyword for t-distribution, change `(tdist)` to `(t)` or `dist=t`.

These checks **do not alter** the paper’s identification or conclusions based on the **ARMA–GARCH(1,1)+DCC(1,1)** baseline.

---

## 5. Mapping to the paper
- **Table 1–3 (Descriptives & ADF):** `sp_desc`, `sp_adf`, plus S5 for PP (`sp_pp`).
- **Table 4–5 (ARMA):** `eq_rcap_mean.output`, `eq_rest_mean.output`, and `tbl_arma`.
- **Table 6 (GARCH):** `g_stock.output`, `g_house.output`, and persistence panel `tbl_garch_pers`.
- **Table 7 (DCC):** `dcc_out` and summary `tbl_dcc` (includes `a`, `b`, `a+b`, mean/min/max/sd of ρ_t).
- **Table 8 (Cheung–Ng):** `tbl_cn_12`; lag-wise correlations in `tbl_cn_lag_12`.
- **Robustness:** `tbl_ccc_dcc`, `tbl_ccc_dcc_summary`, `tbl_cn_multi`, `dcc_out_t`, `tbl_dcc_t`, `tbl_cn_summary`, `tbl_cn_summary_t`.

---

## 6. Methodological references
- Akaike (1974) — AIC  
- Bollerslev (1986) — GARCH(1,1)  
- Bollerslev & Wooldridge (1992) — QML & robust SE (Econometric Reviews, 11(2), 143–172)  
- Cheung & Ng (1996) — Variance causality  
- Dickey & Fuller (1979); Phillips & Perron (1988) — Unit roots  
- Engle (1982) — ARCH; Engle (2002) — DCC(1,1)  
- Ljung & Box (1978) — Portmanteau  
- Schwarz (1978) — BIC  
- White (1982) — QML consistency under misspecification

---

## 7. Reproducibility notes
- The PRG uses `@runpath` to resolve file locations relative to the script.
- Effective sample in Cheung–Ng uses `T_eff = T − m`.
- Conditional correlation stability requires `a + b < 1`; the PRG reports this in `tbl_dcc` and `tbl_dcc_t`.
- If the input file is updated, re-running the PRG regenerates all outputs deterministically.
