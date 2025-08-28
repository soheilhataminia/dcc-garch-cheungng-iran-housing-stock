'====================[ 10) ROBUSTNESS BLOCK ]====================

'---------------------------------------------------------------
' 10.A) CCC vs DCC (summary panel) — uses objects from Section 6
'---------------------------------------------------------------
table tbl_ccc_dcc_summary
tbl_ccc_dcc_summary.setwidth(1) 26
tbl_ccc_dcc_summary.setwidth(2) 18
tbl_ccc_dcc_summary.setwidth(3) 18
tbl_ccc_dcc_summary.setelem(1,1) Model
tbl_ccc_dcc_summary.setelem(1,2) LogLik
tbl_ccc_dcc_summary.setelem(1,3) Note
tbl_ccc_dcc_summary.setelem(2,1) CCC
tbl_ccc_dcc_summary.setelem(2,2) {ccc_eq.@logl}
tbl_ccc_dcc_summary.setelem(2,3) "H0: a=b=0"
tbl_ccc_dcc_summary.setelem(3,1) DCC
tbl_ccc_dcc_summary.setelem(3,2) {dcc_eq.@logl}
tbl_ccc_dcc_summary.setelem(3,3) "LR=2(LL_DCC-LL_CCC), df=2; see tbl_ccc_dcc"

'---------------------------------------------------------------
' 10.B) Cheung–Ng robustness already computed (m=6,12,18)
'       Build a compact summary panel of p-values
'---------------------------------------------------------------
table tbl_cn_summary
tbl_cn_summary.setwidth(1) 28
tbl_cn_summary.setwidth(2) 16
tbl_cn_summary.setwidth(3) 16
tbl_cn_summary.setwidth(4) 16
tbl_cn_summary.setelem(1,1) Direction
tbl_cn_summary.setelem(1,2) p(m=6)
tbl_cn_summary.setelem(1,3) p(m=12)
tbl_cn_summary.setelem(1,4) p(m=18)
tbl_cn_summary.setelem(2,1) Stock → Housing
tbl_cn_summary.setelem(2,2) {p_sh_06}
tbl_cn_summary.setelem(2,3) {p_sh_12}
tbl_cn_summary.setelem(2,4) {p_sh_18}
tbl_cn_summary.setelem(3,1) Housing → Stock
tbl_cn_summary.setelem(3,2) {p_hs_06}
tbl_cn_summary.setelem(3,3) {p_hs_12}
tbl_cn_summary.setelem(3,4) {p_hs_18}

'---------------------------------------------------------------
' 10.C) Student-t QML margins → DCC re-estimation → CN tests
'      (Heavy-tail robustness; leaves baseline intact)
'      NOTE: If your EViews flavor uses a different keyword for t,
'            change "(tdist)" to your local option (e.g., "(t)" or "dist=t").
'---------------------------------------------------------------

' -- GARCH(1,1) with Student-t on ARMA residuals
equation g_stock_t.arch(tdist) e_stock c arch(1) garch(1)
equation g_house_t.arch(tdist) e_house c arch(1) garch(1)

' -- Standardized residuals under t
series ht_stock = g_stock_t.@h
series ht_house = g_house_t.@h
series st_stock = @sqrt(ht_stock)
series st_house = @sqrt(ht_house)
series zt_stock = e_stock / st_stock
series zt_house = e_house / st_house

' -- DCC on t-standardized residuals (normal likelihood is fine here;
'    optionally, you can also try dist=student)
group GZT zt_stock zt_house
equation dcc_eq_t.mgarch(dcc, garch=(1,1), dist=normal, type=diagonal) zt_stock zt_house
freeze(dcc_out_t) dcc_eq_t.output

' -- (Optional) DCC with student likelihood on t-residuals:
' equation dcc_eq_t2.mgarch(dcc, garch=(1,1), dist=student, type=diagonal) zt_stock zt_house

' Extract DCC^t parameters and ρ_t^t
scalar at = @coef("dcc_eq_t","dcc_a")
scalar bt = @coef("dcc_eq_t","dcc_b")

scalar s11t = @var(zt_stock)
scalar s22t = @var(zt_house)
scalar s12t = @cov(zt_stock,zt_house)

series q11t = s11t
series q22t = s22t
series q12t = s12t
smpl @first+1 @last
q11t = (1-at-bt)*s11t + at*zt_stock(-1)^2           + bt*q11t(-1)
q22t = (1-at-bt)*s22t + at*zt_house(-1)^2           + bt*q22t(-1)
q12t = (1-at-bt)*s12t + at*zt_stock(-1)*zt_house(-1)+ bt*q12t(-1)
series rho_t_t = q12t/@sqrt(q11t*q22t)
smpl @all

scalar rho_mean_t = @mean(rho_t_t)
scalar rho_min_t  = @min(rho_t_t)
scalar rho_max_t  = @max(rho_t_t)
scalar rho_sd_t   = @stdev(rho_t_t)

table tbl_dcc_t
tbl_dcc_t.setwidth(1) 22
tbl_dcc_t.setwidth(2) 16
tbl_dcc_t.setelem(1,1) Statistic
tbl_dcc_t.setelem(1,2) Value
tbl_dcc_t.setelem(2,1) a (t-margins)
tbl_dcc_t.setelem(2,2) {at}
tbl_dcc_t.setelem(3,1) b (t-margins)
tbl_dcc_t.setelem(3,2) {bt}
tbl_dcc_t.setelem(4,1) a+b
tbl_dcc_t.setelem(4,2) {at+bt}
tbl_dcc_t.setelem(5,1) Mean(ρ_t^t)
tbl_dcc_t.setelem(5,2) {rho_mean_t}
tbl_dcc_t.setelem(6,1) Min(ρ_t^t)
tbl_dcc_t.setelem(6,2) {rho_min_t}
tbl_dcc_t.setelem(7,1) Max(ρ_t^t)
tbl_dcc_t.setelem(7,2) {rho_max_t}
tbl_dcc_t.setelem(8,1) SD(ρ_t^t)
tbl_dcc_t.setelem(8,2) {rho_sd_t}

' -- Cheung–Ng on t-standardized residuals (m=6,12,18)
vector(6)  rho_sh_06_t
vector(6)  rho_hs_06_t
vector(12) rho_sh_12_t
vector(12) rho_hs_12_t
vector(18) rho_sh_18_t
vector(18) rho_hs_18_t
scalar Q_sh_06_t  Q_hs_06_t  p_sh_06_t  p_hs_06_t
scalar Q_sh_12_t  Q_hs_12_t  p_sh_12_t  p_hs_12_t
scalar Q_sh_18_t  Q_hs_18_t  p_sh_18_t  p_hs_18_t

call cn_test(zt_stock, zt_house, 6)  Q_sh_06_t  p_sh_06_t  rho_sh_06_t
call cn_test(zt_house, zt_stock, 6)  Q_hs_06_t  p_hs_06_t  rho_hs_06_t
call cn_test(zt_stock, zt_house, 12) Q_sh_12_t  p_sh_12_t  rho_sh_12_t
call cn_test(zt_house, zt_stock, 12) Q_hs_12_t  p_hs_12_t  rho_hs_12_t
call cn_test(zt_stock, zt_house, 18) Q_sh_18_t  p_sh_18_t  rho_sh_18_t
call cn_test(zt_house, zt_stock, 18) Q_hs_18_t  p_hs_18_t  rho_hs_18_t

table tbl_cn_summary_t
tbl_cn_summary_t.setwidth(1) 28
tbl_cn_summary_t.setwidth(2) 16
tbl_cn_summary_t.setwidth(3) 16
tbl_cn_summary_t.setwidth(4) 16
tbl_cn_summary_t.setelem(1,1) Direction (t-margins)
tbl_cn_summary_t.setelem(1,2) p(m=6)
tbl_cn_summary_t.setelem(1,3) p(m=12)
tbl_cn_summary_t.setelem(1,4) p(m=18)
tbl_cn_summary_t.setelem(2,1) Stock → Housing
tbl_cn_summary_t.setelem(2,2) {p_sh_06_t}
tbl_cn_summary_t.setelem(2,3) {p_sh_12_t}
tbl_cn_summary_t.setelem(2,4) {p_sh_18_t}
tbl_cn_summary_t.setelem(3,1) Housing → Stock
tbl_cn_summary_t.setelem(3,2) {p_hs_06_t}
tbl_cn_summary_t.setelem(3,3) {p_hs_12_t}
tbl_cn_summary_t.setelem(3,4) {p_hs_18_t}

' -- Add robustness outputs to master spool and export
SP_ALL.append tbl_ccc_dcc_summary
SP_ALL.append tbl_cn_summary
SP_ALL.append dcc_out_t
SP_ALL.append tbl_dcc_t
SP_ALL.append tbl_cn_summary_t

export(table, mode=update) @runpath + "\results_tables_from_eviews.xlsx" _
  tbl_ccc_dcc_summary tbl_cn_summary tbl_dcc_t tbl_cn_summary_t

'====================[ END ROBUSTNESS BLOCK ]====================
