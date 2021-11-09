*********************************************************************
* Combine Enrollment, Premium Data to Backout FPL and Code Inertia
*********************************************************************

/// NOTES /// 

* Available Data Files
	* 1: ca_enroll_hh: Cleaned enrollment data (hh-year)
	* 2: planchars: Plan characteristics (plan-gra-year)
	* 3: inert_xwalk: Crosswalks plan_t-1 to plan_t for auto assign (plan-gra-year)
	
/* Methods to Calc APTC_t1
	- The problem is that we don't see APTCs for HHs that left CCA the next year
	- Different methods to guess what this would be
		1. Use characteristics of all HHs in t0 to impute APTCs at t1
		2. Calculate APTCs based on non-renewing HHs' crosswalked plan, t0 dem chars, bench prem in next year 
			- Can keep the observed premiums for remaining HHs
			- However, doing this is problematic as corr of calc'd prem, actual prem is very low
			- Part of this due to changes in HH composition, income that we can't observe
	- Imputing next year's APTC thus seems like an attractive approach, but is it subject to selection bias?
		- Are non-returning enrollees systematically different than returning enrollees
		- Maybe, but they'd have to be different based on chars a year in advance of their decision to leave
		- Also, decisions to leave/return can be influenced by exogenous shocks (changes in econ)
		- Point is just to predict it as accurately as possible... 
	- Model HH's ATPC in year t+1 as a function of...
		- HH chars in year t: size, age, fpl, aptc, p_pre (or post)
		- Fixed effects for metal level, network, gra, year, maybe network#gra
		- GLM and cluster on net
	*/
	
* Roadmap: Need to determine premium hh would pay next year
	* code already attaches plan id of next year's auto plan
	* but doesn't have mechanism to figure out what post-aptc prem of said plan would be
		* assume income stays the same; cca has to do this
		* eventually will need to kick up age(s) by 1 MAKE THIS ADJUSTMENT IN 1_ENRFORMATTER
	* can mult next year's prem by ageadj to get gross prem, but need aptc
		* first, need prem of bp
		* then can just do raw calc on aptc

/// SETUP ///

clear *
cd "G:\Shared Drives\CovCAInertia"

/// MERGES: PLAN BENEFITS, AUTO ASSIGN CHARS ///

* Load

use "proc\ca_enroll_hh", clear
	drop planid16
	lab var planid "Plan ID"
	lab var csr "CSR Level"

* Merge Plan Characteristics of Selected Plan

merge m:1 year gra planid using proc/planchars, keep(master match) nogen // 6149 unmatched
	drop issuername bp
	lab var prem "Base Premium of Selected Plan"
	rename issuerid issid // check in final data set for variation

* Merge Auto Assignment Info

rename planid planid0
merge m:1 year gra planid0 using "proc/inertxwalk", keep(master match)
	rename planid0 planid
	drop metal0 netid0 issuername0 issuerid0
	tab year _merge // Nearly all 2018; Merge identifies plan to which HH xwalked t+1
	tab issuer_name gra if _merge==1 & year==2014 // 5677 Anthem enrollees 2014 gra 19; minor problem
	tab issuer_name gra if _merge==1 & year==2015 // 2 Anthem, 467 SHARP, 1 Western; too small to matter
	rename _merge inertxwalkmatch
	lab var inertxwalkmatch "Matched to Automatic Assignment X-Walk"
	rename (planid1 netid1 prem1) (autoplanid autonetid autoprem)
	rename (iplanchange inetchange) (autoplanchange autonetchange)
	
* Merge BP for Next Year

merge m:1 year gra using "proc/autobp", keep(master match) nogen // all matched

* Save

save "proc/ca_enroll_hh_merge", replace

/// CHANGE FROM SELECTED PREMIUM TO AUTO ASSIGNED PREMIUM ///

* Load

use "proc/ca_enroll_hh_merge", clear

* Gross Premiums of Current Plan, Next Period Auto Assigned Plan, and Benchmark Plans in t, t+1

g g_sp  = agecurve * prem // 6149 missing; corr(p_pre, g_sp) = 0.46, shoudl be way higher
g g_bp  = agecurve * bp
g g_asp = agecurve * autoprem // 8557 missing
g g_abp = agecurve * autobp

replace g_sp = p_pre if g_sp==. // 6149

lab var g_sp  "Gross Premium: Selected Plan"
lab var g_bp  "Gross Premium: Benchmark Plan"
lab var g_asp "Gross Premium: Auto Assigned Plan Next Year"
lab var g_abp "Gross Premium: Benchmark Plan Next Year"

* Expected Contribution Next Year

g ecpct_t1 = .
lab var ecpct_t1 "Expected Contribution Percentage t+1"

replace ecpct_t1 = 2.08 if fpl < 133 & year==2017 // https://www.thehortongroup.com/resources/affordability-percentages-will-increase-for-2019
replace ecpct_t1 = (((fpl-133)/(150-133)) * (4.03-3.02)) + 3.02 if fpl >= 133 & fpl < 150 & year==2017 
replace ecpct_t1 = (((fpl-150)/(200-150)) * (6.34-4.03)) + 4.03 if fpl >= 150 & fpl < 200 & year==2017 
replace ecpct_t1 = (((fpl-200)/(250-200)) * (8.10-6.34)) + 6.34 if fpl >= 200 & fpl < 250 & year==2017 
replace ecpct_t1 = (((fpl-250)/(300-250)) * (9.56-8.10)) + 8.10 if fpl >= 250 & fpl < 300 & year==2017 
replace ecpct_t1 = 9.56 if fpl >= 300 & fpl <= 400 & year==2017 

replace ecpct_t1 = 2.04 if fpl < 133 & year==2016 // https://www.irs.gov/pub/irs-drop/rp-16-24.pdf
replace ecpct_t1 = (((fpl-133)/(150-133)) * (4.08-3.06)) + 3.06 if fpl >= 133 & fpl < 150 & year==2016
replace ecpct_t1 = (((fpl-150)/(200-150)) * (6.43-4.08)) + 4.08 if fpl >= 150 & fpl < 200 & year==2016
replace ecpct_t1 = (((fpl-200)/(250-200)) * (8.21-6.43)) + 6.43 if fpl >= 200 & fpl < 250 & year==2016
replace ecpct_t1 = (((fpl-250)/(300-250)) * (9.69-8.21)) + 8.21 if fpl >= 250 & fpl < 300 & year==2016
replace ecpct_t1 = 9.69 if fpl >= 300 & fpl <= 400 & year==2016

replace ecpct_t1 = 2.03 if fpl < 133 & year==2015 // https://www.irs.gov/pub/irs-drop/rp-14-62.pdf
replace ecpct_t1 = (((fpl-133)/(150-133)) * (4.07-3.05)) + 3.05 if fpl >= 133 & fpl < 150 & year==2015
replace ecpct_t1 = (((fpl-150)/(200-150)) * (6.41-4.07)) + 4.07 if fpl >= 150 & fpl < 200 & year==2015
replace ecpct_t1 = (((fpl-200)/(250-200)) * (8.18-6.41)) + 6.41 if fpl >= 200 & fpl < 250 & year==2015
replace ecpct_t1 = (((fpl-250)/(300-250)) * (9.66-8.18)) + 8.18 if fpl >= 250 & fpl < 300 & year==2015
replace ecpct_t1 = 9.66 if fpl >= 300 & fpl <= 400 & year==2015

replace ecpct_t1 = 2.01 if fpl < 133 & year==2014 // https://www.irs.gov/pub/irs-drop/rp-14-37.pdf
replace ecpct_t1 = (((fpl-133)/(150-133)) * (4.02-3.02)) + 3.05 if fpl >= 133 & fpl < 150 & year==2014
replace ecpct_t1 = (((fpl-150)/(200-150)) * (6.34-4.02)) + 4.07 if fpl >= 150 & fpl < 200 & year==2014
replace ecpct_t1 = (((fpl-200)/(250-200)) * (8.10-6.34)) + 6.41 if fpl >= 200 & fpl < 250 & year==2014
replace ecpct_t1 = (((fpl-250)/(300-250)) * (9.56-8.10)) + 8.18 if fpl >= 250 & fpl < 300 & year==2014
replace ecpct_t1 = 9.56 if fpl >= 300 & fpl <= 400 & year==2014

replace ecpct_t1 = 2.00 if fpl < 133 & year==2013 // https://www.tribalselfgov.org/wp-content/uploads/2018/08/TSGAC-Brief-Applicable-Percentages-Thresholds-and-Payments-for-ACA-Provisions-2018-07-10c.pdf
replace ecpct_t1 = (((fpl-133)/(150-133)) * (4.00-3.00)) + 3.00 if fpl >= 133 & fpl < 150 & year==2013
replace ecpct_t1 = (((fpl-150)/(200-150)) * (6.30-4.00)) + 4.00 if fpl >= 150 & fpl < 200 & year==2013
replace ecpct_t1 = (((fpl-200)/(250-200)) * (8.05-6.30)) + 6.30 if fpl >= 200 & fpl < 250 & year==2013
replace ecpct_t1 = (((fpl-250)/(300-250)) * (9.50-8.05)) + 8.05 if fpl >= 250 & fpl < 300 & year==2013
replace ecpct_t1 = 9.50 if fpl >= 300 & fpl <= 400 & year==2013

g ecpct_t0 = .
lab var ecpct_t0 "Expected Contribution Percentage t"

replace ecpct_t0 = 2.08 if fpl < 133 & year==2018 // https://www.thehortongroup.com/resources/affordability-percentages-will-increase-for-2019
replace ecpct_t0 = (((fpl-133)/(150-133)) * (4.03-3.02)) + 3.02 if fpl >= 133 & fpl < 150 & year==2018 
replace ecpct_t0 = (((fpl-150)/(200-150)) * (6.34-4.03)) + 4.03 if fpl >= 150 & fpl < 200 & year==2018 
replace ecpct_t0 = (((fpl-200)/(250-200)) * (8.10-6.34)) + 6.34 if fpl >= 200 & fpl < 250 & year==2018 
replace ecpct_t0 = (((fpl-250)/(300-250)) * (9.56-8.10)) + 8.10 if fpl >= 250 & fpl < 300 & year==2018 
replace ecpct_t0 = 9.56 if fpl >= 300 & fpl <= 400 & year==2018 

replace ecpct_t0 = 2.04 if fpl < 133 & year==2017 // https://www.irs.gov/pub/irs-drop/rp-16-24.pdf
replace ecpct_t0 = (((fpl-133)/(150-133)) * (4.08-3.06)) + 3.06 if fpl >= 133 & fpl < 150 & year==2017
replace ecpct_t0 = (((fpl-150)/(200-150)) * (6.43-4.08)) + 4.08 if fpl >= 150 & fpl < 200 & year==2017
replace ecpct_t0 = (((fpl-200)/(250-200)) * (8.21-6.43)) + 6.43 if fpl >= 200 & fpl < 250 & year==2017
replace ecpct_t0 = (((fpl-250)/(300-250)) * (9.69-8.21)) + 8.21 if fpl >= 250 & fpl < 300 & year==2017
replace ecpct_t0 = 9.69 if fpl >= 300 & fpl <= 400 & year==2017

replace ecpct_t0 = 2.03 if fpl < 133 & year==2016 // https://www.irs.gov/pub/irs-drop/rp-14-62.pdf
replace ecpct_t0 = (((fpl-133)/(150-133)) * (4.07-3.05)) + 3.05 if fpl >= 133 & fpl < 150 & year==2016
replace ecpct_t0 = (((fpl-150)/(200-150)) * (6.41-4.07)) + 4.07 if fpl >= 150 & fpl < 200 & year==2016
replace ecpct_t0 = (((fpl-200)/(250-200)) * (8.18-6.41)) + 6.41 if fpl >= 200 & fpl < 250 & year==2016
replace ecpct_t0 = (((fpl-250)/(300-250)) * (9.66-8.18)) + 8.18 if fpl >= 250 & fpl < 300 & year==2016
replace ecpct_t0 = 9.66 if fpl >= 300 & fpl <= 400 & year==2016

replace ecpct_t0 = 2.01 if fpl < 133 & year==2015 // https://www.irs.gov/pub/irs-drop/rp-14-37.pdf
replace ecpct_t0 = (((fpl-133)/(150-133)) * (4.02-3.02)) + 3.05 if fpl >= 133 & fpl < 150 & year==2015
replace ecpct_t0 = (((fpl-150)/(200-150)) * (6.34-4.02)) + 4.07 if fpl >= 150 & fpl < 200 & year==2015
replace ecpct_t0 = (((fpl-200)/(250-200)) * (8.10-6.34)) + 6.41 if fpl >= 200 & fpl < 250 & year==2015
replace ecpct_t0 = (((fpl-250)/(300-250)) * (9.56-8.10)) + 8.18 if fpl >= 250 & fpl < 300 & year==2015
replace ecpct_t0 = 9.56 if fpl >= 300 & fpl <= 400 & year==2015

replace ecpct_t0 = 2.00 if fpl < 133 & year==2014 // https://www.tribalselfgov.org/wp-content/uploads/2018/08/TSGAC-Brief-Applicable-Percentages-Thresholds-and-Payments-for-ACA-Provisions-2018-07-10c.pdf
replace ecpct_t0 = (((fpl-133)/(150-133)) * (4.00-3.00)) + 3.00 if fpl >= 133 & fpl < 150 & year==2014
replace ecpct_t0 = (((fpl-150)/(200-150)) * (6.30-4.00)) + 4.00 if fpl >= 150 & fpl < 200 & year==2014
replace ecpct_t0 = (((fpl-200)/(250-200)) * (8.05-6.30)) + 6.30 if fpl >= 200 & fpl < 250 & year==2014
replace ecpct_t0 = (((fpl-250)/(300-250)) * (9.50-8.05)) + 8.05 if fpl >= 250 & fpl < 300 & year==2014
replace ecpct_t0 = 9.50 if fpl >= 300 & fpl <= 400 & year==2014

* Income 

g inc_t0 = .
lab var inc_t0 "HH Income t+1"

replace inc_t0 = (121.4 * fpl) + ((fsize-1)*4.32) if year==2017 // Need to get raw HH size?
replace inc_t0 = (120.6 * fpl) + ((fsize-1)*4.18) if year==2016
replace inc_t0 = (118.8 * fpl) + ((fsize-1)*4.14) if year==2015
replace inc_t0 = (117.7 * fpl) + ((fsize-1)*4.16) if year==2014
replace inc_t0 = (116.7 * fpl) + ((fsize-1)*4.06) if year==2013

g inc_t1 = .
lab var inc_t1 "HH Income t"

replace inc_t1 = (121.4 * fpl) + ((fsize-1)*4.32) if year==2018 // Need to get raw HH size?
replace inc_t1 = (120.6 * fpl) + ((fsize-1)*4.18) if year==2017
replace inc_t1 = (118.8 * fpl) + ((fsize-1)*4.14) if year==2016
replace inc_t1 = (117.7 * fpl) + ((fsize-1)*4.16) if year==2015
replace inc_t1 = (116.7 * fpl) + ((fsize-1)*4.06) if year==2014

* Back Out Expected Contribution

g ecamt_t1 = inc_t1 * (ecpct_t1/100)
lab var ecamt_t1 "Exp. Contr. Amt. t"

g ecamt_t0 = inc_t0 * (ecpct_t0/100)
lab var ecamt_t0 "Exp. Contr. Amt. t"

* Back Out APTC from Next Year

g aptc_t1 = g_abp - ecamt_t1
lab var aptc_t1 "APTC t+1"

g aptc_t0 = g_bp - ecamt_t0
lab var aptc_t0 "APTC t"

* Back Out Post-APTC Premium of Auto Assigned Plan

g padj_t1 = g_asp - aptc_t1
replace padj_t1 = fsize if padj_t1 < fsize // Price Floor
replace padj_t1 = g_asp if fpl>400
lab var padj_t1 "Post-Subsidy Premium of Auto Assigned Plan t1"

g padj_t0 = g_sp - aptc_t0
replace padj_t0 = fsize if padj_t0 < fsize
replace padj_t0 = g_sp if fpl>400
lab var padj_t0 "Post-Subsidy Premium of Selected Plan t0"

* Back Out Change in Net Premiums, Auto Assigned Plan t+1 from Selected Plan t

g autodp = padj_t0 - padj_t1
replace autodp = g_asp - g_sp if fpl>400 // Gross Premium Diff for those w/o APTC
lab var autodp "Difference in Premiums between Auto Assigned, Selected Plan"

* Review

sum g_* ecpct* inc* ecamt* padj* autodp if fpl<=400
corr padj_t0 p_post if fpl<=400 // 0.3
corr padj_t0 p_post if fpl >400 // 0.18

* Back up Auto Assign Characteristics One Year

* Clean

drop g_* ecpct* inc* ecamt* padj*

* Save 

save "proc/ca_enroll_hh_autodelta", replace

/// OTHER INERTIA-RELATED CHARACTERISTICS ///

* Load 

use "proc/ca_enroll_hh_autodelta", clear

* Lags

sort  idhh year
order idhh year planid

egen household = group(idhh)
tsset household year

g lagprem = p_pre[_n-1]  if household == L.household // Is APTC censored if > gross prem?
g lagaptc = aptc[_n-1]   if household == L.household
g lagplan = planid[_n-1] if household == L.household
g lagnet  = netid[_n-1]  if household == L.household
g lagiss  = issid[_n-1]  if household == L.household
g lagstat = status[_n-1] if household == L.household
replace lagstat = "New" if lagstat==""

lab var lagprem "Gross Premium of Chosen Plan from Previous Year"
lab var lagaptc "APTC from Previous Year"
lab var lagplan "Plan ID from Previous Year"
lab var lagnet  "Net ID from Previous Year"
lab var lagiss  "Iss ID from Previous Year"
lab var lagstat "Enrollment Status in Previous Year"

* Inertia Categories without Automatic Assignment

foreach i in plan net { // iss

g icat`i' = . // Use status variable? Can't have inertia if terminated in previous year

replace icat`i' = 1 if idhh!=idhh[_n-1]                                   // new enrollee
replace icat`i' = 2 if idhh==idhh[_n-1] & `i'id==`i'id[_n-1]              // stayer
replace icat`i' = 3 if idhh==idhh[_n-1] & `i'id!=`i'id[_n-1] & disc`i'==0 // switcher
replace icat`i' = 4 if idhh==idhh[_n-1] & `i'id!=`i'id[_n-1] & disc`i'==1 // discontinued

lab var icat`i' "Inertia Category (No Auto)"
lab def icat`i' 1 "New" 2 "Stayer" 3 "Switcher" 4 "Discontinued"
lab val icat`i' icat`i'

tab icat`i' year

}

* Inertia Categories with Automatic Assignment

foreach i in plan net { // iss

g acat`i' = . 

replace acat`i' = 1 if idhh!=idhh[_n-1]                          // new enr
replace acat`i' = 2 if idhh==idhh[_n-1] & `i'id==auto`i'id[_n-1] // stayer
replace acat`i' = 3 if idhh==idhh[_n-1] & `i'id!=auto`i'id[_n-1] // switcher
replace acat`i' = 4 if idhh==idhh[_n-1] & year[_n-1]==2014 & year==2015 & issid[_n-1]=="99483" // Contra Costa
replace acat`i' = 4 if idhh==idhh[_n-1] & year[_n-1]==2016 & year==2017 & issid[_n-1]=="37873" // United

lab var acat`i' "Inertia Category (Auto)"
lab def acat`i' 1 "New" 2 "Stayer" 3 "Switcher" 4 "Discontinued", replace
lab val acat`i' acat`i'

tab acat`i' year

}

* Decision Categories

g dcat = lagstat=="Enrolled" & icatplan!=4
replace dcat = 2 if lagstat=="Terminated" & icatplan!=4
lab var dcat "Decision Category"
lab def dcat 0 "New Enrollment" 1 "Continuous Enrollment" 2 "Discontinuous Enrollment"
lab val dcat dcat

* Experience

sort idhh year
g iexp = 0
forvalues y = 1/4 {
replace iexp = `y' if idhh==idhh[_n-`y']
}
lab var iexp "Inertia: Experience" 

* Tenure

foreach i in plan net { // iss

g iten`i' = 0
replace iten`i' = 1 if iexp==1
replace iten`i' = 2 if iexp==2 & `i'id[_n-1]==`i'id[_n-2]
replace iten`i' = 3 if iexp==3 & `i'id[_n-1]==`i'id[_n-2] & `i'id[_n-2]==`i'id[_n-3]
replace iten`i' = 4 if iexp==4 & `i'id[_n-1]==`i'id[_n-2] & `i'id[_n-2]==`i'id[_n-3] & `i'id[_n-3]==`i'id[_n-4]
lab var iten`i' "Inertia: `i' Tenure"

}

* Correct Issuer Names

rename  issuer_name issuername
replace issuername = "Anthem"       if issuername=="Anthem Blue Cross"
replace issuername = "Chinese"      if issuername=="Chinese Community"
replace issuername = "Contra Costa" if issuername=="Contra Costa Heal"
replace issuername = "HealthNet"    if issuername=="Health Net"
replace issuername = "Molina"       if issuername=="Molina Health Car"
replace issuername = "Oscar"        if issuername=="Oscar Health Plan"
replace issuername = "Sharp"        if issuername=="SHARP Health Plan"
replace issuername = "United"       if issuername=="UnitedHealthcare"
replace issuername = "Valley"       if issuername=="Valley Health"
replace issuername = "Western"      if issuername=="Western Health"

* Network Name

g netname = issuername + " " + plantype

* Save

save "proc\ca_enroll_hh_wi", replace // CHECK TABS OF DCAT, PSHOCK, LAGPREM, LAGPREMLAG VARIABLES
