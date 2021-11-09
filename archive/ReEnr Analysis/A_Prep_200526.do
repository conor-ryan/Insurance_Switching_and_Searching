/// INITIALIZE ///

clear *
cd "G:\Shared drives\CovCAInertia"

/// PLAN PREP ///

use proc/planchars if year<2018, clear
drop if metal=="Catastrophic"
egen lprem = min(prem), by(year gra)
duplicates drop year gra, force
keep year gra prem
rename prem lplead
save proc/lpremlead, replace

/// ENROLMENT PREP ///

* Load

use "proc/ca_enroll_hh_merge", clear
sort idhh year

* Week Filter

g wk_end = substr(cov_end,-2,2)
destring wk_end, ignore(".") replace

g wk_strt = substr(cov_strt,-2,2)
destring wk_strt, ignore(".") replace

* Reenrollment, Auto Indicators

g reenr = (wk_end==. | wk_end>=45) & wk_strt[_n+1]<=8 & idhh==idhh[_n+1]
g autod = reenr==1 & act_pass=="A" // Maybe base this on auto elig?

* Filters

drop if year==2018
drop if wk_end<45 & wk_end!=. // Terminated before November

* Auto Premium Variables

replace p_post = fsize if p_post<5

g auto_p_pre = agecurve * autoprem // SETTING WRT CURRENT YEAR SINCE THAT'S WHAT'S IN LETTER
lab var auto_p_pre "Pre-Subsidy, Age-Adjusted Premium of Auto-Assigned Plan"

g auto_p_post = auto_p_pre - aptc // DITTO
replace auto_p_post = fsize if auto_p_post < fsize // ADJUST FOR APTC > PREMIUM
lab var auto_p_post "Post-Subsidy, Age-Adjusted Premium of Auto-Assigned Plan"
rename auto_p_post next_defp

g autodp = next_defp - p_post // Missing vals are contra costa, united, and 2018
lab var autodp "Change in Post-Subsidy Premium: Auto Assigned Plan re Last Year's Plan"

* Premium of Lowest Cost Plan Next Year

merge m:1 year gra using proc/lpremlead, nogen // all match

g next_lpps = lplead * fageadj
lab var next_lpps "Age Adjusted Premium of Lowest Cost Plan Next Year"

* Format Premium Variables

rename autoisschange x_switch

g delta_raw = next_defp - p_post

sum next_* delta_raw p_post

* Fixed Effects

encode metal, g(fe_mtl)
lab var fe_mtl "Metal FE"

egen fe_age = cut(fage), at(18,35,45,55,150) icodes
lab var fe_age "Age FE"

egen fe_fpl = cut(fpl), at(138,150,200,250,300,401) icodes // filters regression
lab var fe_fpl "Income FE"

* Group FE

egen fe_dem    = group(fe_age fe_fpl fsize)
egen fe_plan   = group(planid)
egen fe_planyr = group(planid year)
egen fe_grayr  = group(gra year)

* IHS

g next_defp_ihs = asinh(next_defp)
g next_lpps_ihs = asinh(next_lpps)
g delta_raw_ihs = asinh(delta_raw)

* Windsor

sum next_defp, d
g next_defp_w95 = min(next_defp, r(p95))

sum next_lpps, d
g next_lpps_w95 = min(next_lpps, r(p95))

sum delta_raw, d
g delta_raw_w95 = min(delta_raw, r(p95))
replace delta_raw_w95 = max(delta_raw_w95,r(p5))

/// REGRESSION ///

global xb i.fe_age i.fe_fpl i.fsize i.fe_mtl

reghdfe reenr next_defp_w95 next_lpps_w95 x_switch, a(fe_grayr fe_dem fe_mtl) vce(cluster gra)
reghdfe reenr delta_raw_w95 next_lpps_w95 x_switch, a(fe_grayr fe_dem fe_mtl) vce(cluster gra)

reghdfe reenr next_defp     next_lpps     x_switch, a(fe_grayr fe_dem fe_mtl) vce(cluster gra)
reghdfe reenr delta_raw     next_lpps     x_switch, a(fe_mtl fe_dem fe_grayr) vce(cluster gra)
reghdfe reenr next_defp_ihs next_lpps_ihs x_switch, a(fe_mtl fe_dem fe_grayr) vce(cluster gra)
reghdfe reenr delta_raw_ihs next_lpps_ihs x_switch, a(fe_mtl fe_dem fe_grayr) vce(cluster gra)

/* NEXT STEPS 
- p_post and auto_p_post both had some huge negative values... why???
- make sure auto referring to subsequent year
- x_switch should not be booted due to collinearity; can redefine to be bcbs/2017