*********************************************************************
* Combine Enrollment, Premium Data to Backout FPL and Code Inertia
*********************************************************************

/// MERGES: PLAN BENEFITS, AUTO ASSIGN CHARS ///

* Load

use "proc\ca_enroll_hh", clear
	drop planid16 
	rename (metal plantype plan_name issuer_name) (s_metal s_plantype s_plan_name s_issuer_name) // Remove fields that we want to have be determined by plan char data
	lab var planid "Plan ID"
	lab var csr "CSR Level"

* Merge Plan Characteristics of Selected Plan

merge m:1 year gra planid using proc/planchars, keep(master match) // 6149 unmatched; same as below
	drop if _merge==1 // 6149
	drop _merge s_metal s_plantype s_plan_name s_issuer_name
	lab var prem "Base Premium of Selected Plan"
	rename issuerid issid // check in final data set for variation
	
* Drop concurrent mean stats 
drop mean_prem_all mean_prem_mtl

* Merge Auto Assignment Info

rename planid planid0
merge m:1 year gra planid0 using "proc/inertxwalk", keep(master match)
	rename planid0 planid
	tab year _merge // Nearly all 2018; Merge identifies plan to which HH xwalked t+1
	tab issuername gra if _merge==1 & year==2014 // 5677 Anthem enrollees 2014 gra 19; minor problem
	tab issuername gra if _merge==1 & year==2015 // 2 Anthem, 467 SHARP, 1 Western; too small to matter
	rename _merge inertxwalkmatch
	lab var inertxwalkmatch "Matched to Automatic Assignment X-Walk"
	rename (auto_planid auto_net auto_prem auto_iss) (autoplanid autonet autoprem autoiss)  // Should just name these in code (3) but I'm lazy
	rename (iplanchange inetchange iisschange) (autoplanchange autonetchange autoisschange) // Ditto

* Save

save "proc/ca_enroll_hh_merge", replace

/// INERTIA CHARACTERISTICS ///

* Load

use "proc/ca_enroll_hh_merge", clear

* Change in Premiums, Selected t vs. Auto t+1

sort idhh year

g p_post_approx = fageadj*prem - aptc //Could potentially use this to avoid measurement error

g auto_p_pre = fageadj[_n+1] * autoprem if idhh==idhh[_n+1] /*& year[_n+1]==year+1*/
lab var auto_p_pre "Pre-Subsidy, Age-Adjusted Premium of Auto-Assigned Plan"

g auto_p_post = auto_p_pre - aptc[_n+1] if idhh==idhh[_n+1] /*& year[_n+1]==year+1*/
lab var auto_p_post "Post-Subsidy, Age-Adjusted Premium of Auto-Assigned Plan"

// g autodp = auto_p_post - p_post if idhh==idhh[_n+1] /*& year[_n+1]==year+1*/
// lab var autodp "Change in Post-Subsidy Premium: Auto Assigned Plan re Last Year's Plan"

g autodp = auto_p_post - p_post_approx if idhh==idhh[_n+1] /*& year[_n+1]==year+1*/
lab var autodp "Change in Post-Subsidy Premium: Auto Assigned Plan re Last Year's Plan"

count if autodp==. & idhh==idhh[_n+1] // 100,997
count if autodp==. & idhh==idhh[_n+1]& year[_n+1]==year+1 // 923

* Move Auto Variables Forward One Year

sort  idhh year
order idhh year planid auto* agecurve p_post

foreach y of varlist autoplanid autonet autoiss autoprem autoplanchange-autodp base_dprem mean_dprem_all mean_dprem_mtl auto_mean_prem_all auto_mean_prem_mtl { 
	local i = `i' + 1
	gen `y'_new = `y'[_n-1] if idhh==idhh[_n-1] 
	if `i'<=3 replace `y'_new = "" if year==2014 | idhh!=idhh[_n-1] | year[_n-1]!=year-1
	if `i' >3 replace `y'_new = .  if year==2014 | idhh!=idhh[_n-1] | year[_n-1]!=year-1
	replace `y' = `y'_new
	drop `y'_new
} // note: foreach used to contain autoplanchange-audodp_approx, but autodp_approx wasn't here anymore 2/25

* Lags

sort  idhh year
order idhh year planid

g daptc = aptc - aptc[_n-1] if idhh==idhh[_n-1]
lab var autodp "Change in APTC"


***** Test Network Classififcations ******

replace netname = "Anthem PPO"      if netname== "Anthem EPO"
replace product = "PPO"             if product=="EPO" & issuername=="Anthem"
replace netname = "HealthNet EPO"   if netname=="HealthNet PPO"
replace product = "EPO"             if product=="PPO" & issuername=="HealthNet"
replace netname = "Blue Shield PPO" if netname=="Blue Shield EPO" 
replace product = "PPO"             if product=="EPO" & issuername=="Blue Shield"

g lagppre   = p_pre[_n-1]        if idhh == idhh[_n-1] & year[_n-1]==year-1 // Is APTC censored if > gross prem?
g lagppost  = p_post[_n-1]       if idhh == idhh[_n-1]& year[_n-1]==year-1
g lagaptc   = aptc[_n-1]         if idhh == idhh[_n-1]& year[_n-1]==year-1
g lagplan   = planid[_n-1]       if idhh == idhh[_n-1]& year[_n-1]==year-1
g lagnet    = netname[_n-1]      if idhh == idhh[_n-1]& year[_n-1]==year-1
g lagiss    = issuername[_n-1]   if idhh == idhh[_n-1]& year[_n-1]==year-1
g lagmetal  = metal[_n-1]        if idhh == idhh[_n-1]& year[_n-1]==year-1
g lagproduct  = product[_n-1]    if idhh == idhh[_n-1]& year[_n-1]==year-1
g lagcovend = cov_end[_n-1]      if idhh == idhh[_n-1]& year[_n-1]==year-1
g lagstat   = status[_n-1]       if idhh == idhh[_n-1]& year[_n-1]==year-1
replace lagstat = "New" if lagstat==""

lab var lagppre   "Gross Premium of Chosen Plan from Previous Year"
lab var lagppost  "Net Premium of Chosen Plan from Previous Year"
lab var lagaptc   "APTC from Previous Year"
lab var lagplan   "Plan ID from Previous Year"
lab var lagnet    "Net Name from Previous Year"
lab var lagiss    "Iss Name from Previous Year"
lab var lagstat   "Enrollment Status in Previous Year"
lab var lagcovend "Coverage End Week from Previous Year"
lab var lagmetal  "Metal Level from Previous Year"
lab var lagproduct  "Product Line from Previous Year"

// * Test Code
// gen age_flag = fage>70
// egen any_age_flag = sum(age_flag), by(idhh)
// gen excess_change = autodp + daptc
// gen flag = excess_change!=. & abs(excess_change)>2000 
// egen anyflag = sum(flag), by(idhh)
// g lagfadj   = fageadj[_n-1]  if idhh == idhh[_n-1]& year[_n-1]==year-1
// g dfadj = fageadj - lagfadj

* Inertia Categories

sort idhh year
order idhh year netname

g icatplan = .
replace icatplan = 1 if idhh!=idhh[_n-1] | year[_n-1]<year-1
replace icatplan = 2 if idhh==idhh[_n-1] & netname==netname[_n-1] & metal==metal[_n-1] & product==product[_n-1] & year== year[_n-1]+1
replace icatplan = 3 if idhh==idhh[_n-1] & (netname!=netname[_n-1] | metal!=metal[_n-1] | product!=product[_n-1]) & year== year[_n-1]+1
lab var icatplan "Plan Inertia Category"
lab def icatplan 1 "New" 2 "Stayer" 3 "Switcher", replace
lab val icatplan icatplan
tab icatplan year, m

g icatnet = .
replace icatnet = 1 if idhh!=idhh[_n-1] | year[_n-1]<year-1
replace icatnet = 2 if idhh==idhh[_n-1] & netname==netname[_n-1] & year== year[_n-1]+1
replace icatnet = 3 if idhh==idhh[_n-1] & netname!=netname[_n-1] & year== year[_n-1]+1
lab var icatnet "Plan Inertia Category"
lab def icatnet 1 "New" 2 "Stayer" 3 "Switcher", replace
lab val icatnet icatnet
tab icatnet year, m

* Adjust for Enrollees not Having Inertia if Terminated Coverage

tab icatplan if lagstat=="Terminated" // 30% switchers had previously terminated, only ~10% of stayers
count if lagcovend!=".." & regexm(lagcovend,"W52")!=1 & regexm(lagcovend,"W53")!=1 & idhh==idhh[_n-1] // 209k
count if act_pass=="" & idhh==idhh[_n-1] // 196k; go w above approach as it's not dependent on a variable that seems a little questionable

rename (icatplan icatnet) (icatplan3 icatnet3)
lab var icatplan3 "Plan Inertia Category: No Adjustment for Previous Termination"
lab var icatnet3   "Net Inertia Category: No Adjustment for Previous Termination"

g icatplan = icatplan3
replace icatplan = 4 if lagcovend!=".." & regexm(lagcovend,"W52")!=1 & regexm(lagcovend,"W53")!=1 & idhh==idhh[_n-1] 
lab def icatplan 1 "New" 2 "Stayer" 3 "Switcher" 4 "Returner", replace
lab val icatplan icatplan
tab icatplan3 year, m
tab icatplan  year, m

g icatnet = icatnet3
replace icatnet  = 4 if lagcovend!=".." & regexm(lagcovend,"W52")!=1 & regexm(lagcovend,"W53")!=1 & idhh==idhh[_n-1]
lab def icatnet 1 "New" 2 "Stayer" 3 "Switcher" 4 "Returner", replace
lab val icatnet icatnet
tab icatnet3 year, m
tab icatnet  year, m

* Auto Eligibility

g autoelig = icatnet==2 | icatnet==3
lab var autoelig "Eligible for Automatic Reenrollment"
assert autoelig==0 if idhh!=idhh[_n-1]

* Drop Autoplanid if the household is not eligible
replace autoplanid="" if autoelig==0

* Categorize enrollees based on whether they stay or switch from their auto option
g icatauto = icatplan3
replace icatauto=2 if (planid==autoplanid & autoelig==1)
lab def icatauto 1 "New" 2 "Stayer" 3 "Switcher", replace
lab val icatauto icatauto

tab icatauto year if autoplanchange==1

tab icatauto autoisschange

*251 obs of staying on same plan even though auto plan is new. 
*Just dropping these households for now. (942 hh - years)
gen flag = icatplan==2&autoplanchange==1
egen hhflag = sum(flag), by(idhh)

drop if hhflag==1
drop flag hhflag

* Drop if household is autoeligible but missing autodp
*806 obs, 2,645 hh-years)
gen flag = autodp==. & autoelig==1
egen hhflag = sum(flag), by(idhh)
drop if hhflag==1
drop flag hhflag

* Cleaning

compress
drop agecurve netid

order year gra idhh ///
issuername issid netname plantype planid metal csr ///
prem p_pre p_post aptc hassub pmtlord ///
fpl fplcat fage age3 fageadj fam fsize svch ///
lag* auto* icat* inertxwalkmatch status act_pass enr_strt cov_strt cov_end

* Save Version for HH Analysis

save "proc\ca_enroll_hh_wi", replace
