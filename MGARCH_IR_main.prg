'============================================================
' Reproducible EViews Program (PRG)
' Paper: Volatility Spillovers Between Housing and Stock Markets (Iran)
' Sample: 2016M04–2023M10 (monthly)
' Pipeline: ARMA → GARCH(1,1) → DCC(1,1) → Cheung–Ng (m=12) + PP
' EViews: v13+ (MGARCH DCC supported)
'============================================================

'------------------------------
' 0) Workfile & Data Import
'------------------------------
wfcreate(wf=MGARCH_IR, page=Main) m 2016m04 2023m10

' Place the Excel file next to this PRG:
' Sheet "gregorian" with columns:
'   date_gregorian (YYYY-MM), stock_index_close, tehran_price_per_m2_rial
%datapath = @runpath + "\mgarch_gregorian_unmodified.xlsx"

import(options) {%datapath} range=gregorian @date(date_gregorian,"YYYY-MM") stock_index_close tehran_price_per_m2_rial
rename stock_index_close           capidx
rename tehran_price_per_m2_rial    estp

'------------------------------
' 1) Returns (log * 100) & sample alignment
'------------------------------
series rcap = 100*(log(capidx) - log(capidx(-1)))
series rest = 100*(log(estp)  - log(estp(-1)))

smpl @all
smpl if @isna(rcap)=0 and @isna(rest)=0   ' keep non-missing overlap

'------------------------------
' 2) Descriptives, ADF, PP
'------------------------------
spool sp_desc
spool sp_adf
spool sp_pp

sp_desc.append rcap.stats
sp_desc.append rest.stats

' ADF (intercept only; level)
equation adf_rcap.adf rcap 0 0 c
equation adf_rest.adf rest 0 0 c
sp_adf.append adf_rcap.output
sp_adf.append adf_rest.output

' Phillips–Perron (level, intercept only) for SI
freeze(pp_rcap)  rcap.unitroot(pp, c)
freeze(pp_rest)  rest.unitroot(pp, c)
sp_pp.append pp_rcap
sp_pp.append pp_rest

'------------------------------
' 3) ARMA order selection (AIC; p,q ∈ {0,1,2})
'------------------------------
subroutine buildspec(string %y, scalar !p, scalar !q) string %spec
  %spec = %y + " c"
  if !p>0 then
    for !i=1 to !p
      %spec = %spec + " ar(" + @str(!i) + ")"
    next
  endif
  if !q>0 then
    for !j=1 to !q
      %spec = %spec + " ma(" + @str(!j) + ")"
    next
  endif
endsub

' RCAP
scalar !best_aic_rc = 1e+30
scalar !best_p_rc  = 0
scalar !best_q_rc  = 0
for !p=0 to 2
  for !q=0 to 2
    call buildspec("rcap", !p, !q) %s
    equation eqtmp1.ls {%s}
    if @aic < !best_aic_rc then
      !best_aic_rc = @aic
      !best_p_rc = !p
      !best_q_rc = !q
      eq_rcap_mean.copy eqtmp1
    endif
  next
next

' REST
scalar !best_aic_re = 1e+30
scalar !best_p_re  = 0
scalar !best_q_re  = 0
for !p=0 to 2
  for !q=0 to 2
    call buildspec("rest", !p, !q) %s
    equation eqtmp2.ls {%s}
    if @aic < !best_aic_re then
      !best_aic_re = @aic
      !best_p_re = !p
      !best_q_re = !q
      eq_rest_mean.copy eqtmp2
    endif
  next
next

' ARMA order summary
table tbl_arma
tbl_arma.setwidth(1) 26
tbl_arma.setwidth(2) 18
tbl_arma.setwidth(3) 18
tbl_arma.setwidth(4) 18
tbl_arma.setelem(1,1) Series
tbl_arma.setelem(1,2) Selected ARMA(p,q)
tbl_arma.setelem(1,3) AIC
tbl_arma.setelem(1,4) BIC
tbl_arma.setelem(2,1) RCAP
tbl_arma.setelem(2,2) ({!best_p_rc}, {!best_q_rc})
tbl_arma.setelem(2,3) {!best_aic_rc}
tbl_arma.setelem(2,4) {eq_rcap_mean.@schwarz}
tbl_arma.setelem(3,1) REST
tbl_arma.setelem(3,2) ({!best_p_re}, {!best_q_re})
tbl_arma.setelem(3,3) {!best_aic_re}
tbl_arma.setelem(3,4) {eq_rest_mean.@schwarz}

'------------------------------
' 4) ARMA residual diagnostics (12 lags)
'------------------------------
series e_stock = eq_rcap_mean.resid
series e_house = eq_rest_mean.resid

freeze(lb_resid_rc)  e_stock.qstat(12)      ' Ljung–Box on returns
freeze(lb_resid_re)  e_house.qstat(12)
freeze(lb_sq_rc)     (e_stock^2).qstat(12)  ' Ljung–Box on squares
freeze(lb_sq_re)     (e_house^2).qstat(12)

freeze(archlm_rc)    eq_rcap_mean.archlm(12)
freeze(archlm_re)    eq_rest_mean.archlm(12)

'------------------------------
' 5) Univariate GARCH(1,1) on ARMA residuals (Gaussian QML)
'------------------------------
equation g_stock.arch e_stock c arch(1) garch(1)
equation g_house.arch e_house c arch(1) garch(1)

' Conditional variances & standardized residuals
series h_stock = g_stock.@h
series h_house = g_house.@h
series s_stock = @sqrt(h_stock)
series s_house = @sqrt(h_house)
series z_stock = e_stock / s_stock
series z_house = e_house / s_house

' Persistence and half-life (months)
scalar pers_stock = @coef("g_stock","c(3)") + @coef("g_stock","c(4)")   ' α+β if c(3)=arch(1), c(4)=garch(1)
scalar pers_house = @coef("g_house","c(3)") + @coef("g_house","c(4)")
scalar hl_stock   = @log(0.5)/@log(pers_stock)
scalar hl_house   = @log(0.5)/@log(pers_house)

table tbl_garch_pers
tbl_garch_pers.setwidth(1) 20
tbl_garch_pers.setwidth(2) 16
tbl_garch_pers.setwidth(3) 16
tbl_garch_pers.setelem(1,1) Series
tbl_garch_pers.setelem(1,2) Alpha+Beta
tbl_garch_pers.setelem(1,3) Half-life (months)
tbl_garch_pers.setelem(2,1) RCAP
tbl_garch_pers.setelem(2,2) {pers_stock}
tbl_garch_pers.setelem(2,3) {hl_stock}
tbl_garch_pers.setelem(3,1) REST
tbl_garch_pers.setelem(3,2) {pers_house}
tbl_garch_pers.setelem(3,3) {hl_house}

'------------------------------
' 6) DCC(1,1): two-step on standardized residuals
'------------------------------
group GZ z_stock z_house

' DCC(1,1) on z's (normal)
equation dcc_eq.mgarch(dcc, garch=(1,1), dist=normal, type=diagonal) z_stock z_house
freeze(dcc_out) dcc_eq.output

' CCC vs DCC (for necessity test)
equation ccc_eq.mgarch(ccc, garch=(1,1), dist=normal, type=diagonal) z_stock z_house
scalar LR_ccc_dcc = 2*(dcc_eq.@logl - ccc_eq.@logl)
scalar p_LR = 1 - @cchisq(LR_ccc_dcc, 2)

table tbl_ccc_dcc
tbl_ccc_dcc.setwidth(1) 24
tbl_ccc_dcc.setwidth(2) 20
tbl_ccc_dcc.setwidth(3) 18
tbl_ccc_dcc.setelem(1,1) Test
tbl_ccc_dcc.setelem(1,2) Statistic
tbl_ccc_dcc.setelem(1,3) p-value
tbl_ccc_dcc.setelem(2,1) LR (CCC vs DCC)
tbl_ccc_dcc.setelem(2,2) {LR_ccc_dcc}
tbl_ccc_dcc.setelem(2,3) {p_LR}

' Extract DCC parameters
scalar a = @coef("dcc_eq","dcc_a")
scalar b = @coef("dcc_eq","dcc_b")

' Reconstruct Q_t and ρ_t (for stats and plots)
scalar s11 = @var(z_stock)
scalar s22 = @var(z_house)
scalar s12 = @cov(z_stock,z_house)

series q11 = s11
series q22 = s22
series q12 = s12
smpl @first+1 @last
q11 = (1-a-b)*s11 + a*z_stock(-1)^2            + b*q11(-1)
q22 = (1-a-b)*s22 + a*z_house(-1)^2            + b*q22(-1)
q12 = (1-a-b)*s12 + a*z_stock(-1)*z_house(-1)  + b*q12(-1)
series rho_t = q12/@sqrt(q11*q22)
smpl @all

' Summary stats for ρ_t
scalar rho_mean = @mean(rho_t)
scalar rho_min  = @min(rho_t)
scalar rho_max  = @max(rho_t)
scalar rho_sd   = @stdev(rho_t)
scalar rho_neg_share = @mean(rho_t<0)

table tbl_dcc
tbl_dcc.setwidth(1) 20
tbl_dcc.setwidth(2) 16
tbl_dcc.setelem(1,1) Statistic
tbl_dcc.setelem(1,2) Value
tbl_dcc.setelem(2,1) a
tbl_dcc.setelem(2,2) {a}
tbl_dcc.setelem(3,1) b
tbl_dcc.setelem(3,2) {b}
tbl_dcc.setelem(4,1) a+b
tbl_dcc.setelem(4,2) {a+b}
tbl_dcc.setelem(5,1) Mean(ρ_t)
tbl_dcc.setelem(5,2) {rho_mean}
tbl_dcc.setelem(6,1) Min(ρ_t)
tbl_dcc.setelem(6,2) {rho_min}
tbl_dcc.setelem(7,1) Max(ρ_t)
tbl_dcc.setelem(7,2) {rho_max}
tbl_dcc.setelem(8,1) SD(ρ_t)
tbl_dcc.setelem(8,2) {rho_sd}
tbl_dcc.setelem(9,1) Share(ρ_t<0)
tbl_dcc.setelem(9,2) {rho_neg_share}

'------------------------------
' 7) Cheung–Ng (1996) variance-causality
'    m = 12 (main), plus m = 6 & 18 (robustness)
'------------------------------
subroutine cn_test(series s_from, series s_to, scalar !m) scalar !Q, scalar !p, vector v_rho
  ' s_from^2 (lagged) → s_to^2 contemporaneous, up to m
  series u = s_from^2 - @mean(s_from^2)
  series v = s_to^2   - @mean(s_to^2)
  vector(!m) rho
  scalar Tobs = @obs(u)
  for !k=1 to !m
    smpl @first+!k @last
    series u_lag = u(-!k)
    rho(!k) = @cor(v, u_lag)
  next
  smpl @all
  scalar T_eff = Tobs - !m       ' <<< robust finite-sample effective T
  !Q = T_eff * @inner(rho, rho)
  !p = 1 - @cchisq(!Q, !m)
  v_rho.copy rho
endsub

' m = 12 (main)
vector(12) rho_sh_12
vector(12) rho_hs_12
scalar Q_sh_12
scalar Q_hs_12
scalar p_sh_12
scalar p_hs_12

call cn_test(z_stock, z_house, 12) Q_sh_12 p_sh_12 rho_sh_12
call cn_test(z_house, z_stock, 12) Q_hs_12 p_hs_12 rho_hs_12

table tbl_cn_12
tbl_cn_12.setwidth(1) 24
tbl_cn_12.setwidth(2) 20
tbl_cn_12.setwidth(3) 16
tbl_cn_12.setelem(1,1) Direction
tbl_cn_12.setelem(1,2) Q-stat (df=12)
tbl_cn_12.setelem(1,3) p-value
tbl_cn_12.setelem(2,1) Stock → Housing
tbl_cn_12.setelem(2,2) {Q_sh_12}
tbl_cn_12.setelem(2,3) {p_sh_12}
tbl_cn_12.setelem(3,1) Housing → Stock
tbl_cn_12.setelem(3,2) {Q_hs_12}
tbl_cn_12.setelem(3,3) {p_hs_12}

' Lag-wise table (m=12)
table tbl_cn_lag_12
tbl_cn_lag_12.setwidth(1) 6
tbl_cn_lag_12.setwidth(2) 28
tbl_cn_lag_12.setwidth(3) 28
tbl_cn_lag_12.setelem(1,1) Lag
tbl_cn_lag_12.setelem(1,2) Corr[stock^2→housing]
tbl_cn_lag_12.setelem(1,3) Corr[housing^2→stock]
for !k=1 to 12
  tbl_cn_lag_12.setelem(1+!k,1) {!k}
  tbl_cn_lag_12.setelem(1+!k,2) {rho_sh_12(!k)}
  tbl_cn_lag_12.setelem(1+!k,3) {rho_hs_12(!k)}
next

' Robustness m = 6 and m = 18
vector(6)  rho_sh_06
vector(6)  rho_hs_06
vector(18) rho_sh_18
vector(18) rho_hs_18
scalar Q_sh_06  Q_hs_06  p_sh_06  p_hs_06
scalar Q_sh_18  Q_hs_18  p_sh_18  p_hs_18

call cn_test(z_stock, z_house, 6)  Q_sh_06 p_sh_06 rho_sh_06
call cn_test(z_house, z_stock, 6)  Q_hs_06 p_hs_06 rho_hs_06
call cn_test(z_stock, z_house, 18) Q_sh_18 p_sh_18 rho_sh_18
call cn_test(z_house, z_stock, 18) Q_hs_18 p_hs_18 rho_hs_18

table tbl_cn_multi
tbl_cn_multi.setwidth(1) 22
tbl_cn_multi.setwidth(2) 16
tbl_cn_multi.setwidth(3) 16
tbl_cn_multi.setwidth(4) 16
tbl_cn_multi.setelem(1,1) Direction
tbl_cn_multi.setelem(1,2) Q(df=6)
tbl_cn_multi.setelem(1,3) Q(df=12)
tbl_cn_multi.setelem(1,4) Q(df=18)
tbl_cn_multi.setelem(2,1) Stock → Housing (p-values)
tbl_cn_multi.setelem(2,2) {p_sh_06}
tbl_cn_multi.setelem(2,3) {p_sh_12}
tbl_cn_multi.setelem(2,4) {p_sh_18}
tbl_cn_multi.setelem(3,1) Housing → Stock (p-values)
tbl_cn_multi.setelem(3,2) {p_hs_06}
tbl_cn_multi.setelem(3,3) {p_hs_12}
tbl_cn_multi.setelem(3,4) {p_hs_18}

'------------------------------
' 8) Figures: conditional σ_t and DCC ρ_t
'------------------------------
graph fig_vol.line s_stock s_house
fig_vol.addtext(t) "Figure 1. Conditional Volatilities (GARCH(1,1))"
fig_vol.legend on
fig_vol.axis(l) label "Sigma"

graph fig_dcc.line rho_t
fig_dcc.addtext(t) "Figure 2. DCC Dynamic Correlation: Stock vs Housing"
fig_dcc.legend off
fig_dcc.axis(l) label "Correlation"

' Save PNGs next to PRG
fig_vol.save(t=png) @runpath + "\Figure1_conditional_vols.png"
fig_dcc.save(t=png) @runpath + "\Figure2_dcc_correlation.png"

'------------------------------
' 9) Spool & Export
'------------------------------
spool SP_ALL
SP_ALL.append sp_desc
SP_ALL.append sp_adf
SP_ALL.append sp_pp
SP_ALL.append eq_rcap_mean.output
SP_ALL.append eq_rest_mean.output
SP_ALL.append tbl_arma
SP_ALL.append lb_resid_rc
SP_ALL.append lb_resid_re
SP_ALL.append lb_sq_rc
SP_ALL.append lb_sq_re
SP_ALL.append archlm_rc
SP_ALL.append archlm_re
SP_ALL.append g_stock.output
SP_ALL.append g_house.output
SP_ALL.append tbl_garch_pers
SP_ALL.append dcc_out
SP_ALL.append tbl_ccc_dcc
SP_ALL.append tbl_dcc
SP_ALL.append tbl_cn_12
SP_ALL.append tbl_cn_lag_12
SP_ALL.append tbl_cn_multi
SP_ALL.append fig_vol
SP_ALL.append fig_dcc

' Export tables/views to Excel
export(table, mode=overwrite) @runpath + "\results_tables_from_eviews.xlsx" _
  tbl_arma tbl_garch_pers tbl_dcc tbl_ccc_dcc tbl_cn_12 tbl_cn_lag_12 tbl_cn_multi pp_rcap pp_rest

' Save workfile
save @runpath + "\MGARCH_IR_EViews.wf1"

'================== End of Program ==================
