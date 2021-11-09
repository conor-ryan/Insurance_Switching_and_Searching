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

use "proc\ca_enroll_hh_wi", clear // Builds off of 4_inertCat_191202

* Week Filter

g wk_end = substr(cov_end,-2,2)
destring wk_end, ignore(".") replace

g wk_strt = substr(cov_strt,-2,2)
destring wk_strt, ignore(".") replace

* Reenrollment, Auto Indicators

g reenr = (wk_end==. | wk_end>=45) & wk_strt<=8 & idhh==idhh[_n+1]
g autod = reenr==1 & act_pass=="A" // Maybe base this on auto elig?

* Filters

drop if year==2018
drop if wk_end<45 & wk_end!=. // Terminated before November

* Premium of Lowest Cost Plan Next Year

merge m:1 year gra using proc/lpremlead, nogen // all match

g next_lpps = lplead * fageadj
lab var next_lpps "Age Adjusted Premium of Lowest Cost Plan Next Year"

* Format Premium Variables

replace p_post = fsize if p_post<5

replace auto_p_post = fsize if auto_p_post<5
rename auto_p_post next_defp

rename autoisschange x_switch

g delta_raw = next_defp - p_post

sum next_* delta_raw p_post

* Fixed Effects

encode metal, g(fe_mtl)
lab var fe_mtl "Metal FE"

egen fe_age = cut(fage), at(18,35,45,55,150) 
lab var fe_age "Age FE"

egen fe_fpl = cut(fpl), at(138,150,200,250,300,401)
lab var fe_fpl "Income FE"

* Group FE

egen fe_dem    = group(fe_age fe_fpl fsize)
egen fe_plan   = group(planid)
egen fe_planyr = group(planid year)
egen fe_grayr  = group(gra year)

/// REGRESSION ///

reghdfe reenr next_defp next_lpps x_switch, a(fe_mtl fe_dem fe_grayr) vce(cluster gra)
reghdfe reenr delta_raw next_lpps x_switch, a(fe_mtl fe_dem fe_grayr) vce(cluster gra)

/* NEXT STEPS 
- p_post and auto_p_post both had some huge negative values... why???
- make sure auto referring to subsequent year
- x_switch should not be booted due to collinearity; can redefine to be bcbs/2017