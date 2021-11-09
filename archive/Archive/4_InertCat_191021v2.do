*********************************************************************
* Combine Enrollment, Premium Data to Backout FPL and Code Inertia
*********************************************************************

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

merge m:1 year gra planid using proc/planchars, keep(master match) nogen // 6149 unmatched; same as below
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

* Save

save "proc/ca_enroll_hh_merge", replace

/// INERTIA CHARACTERISTICS ///

* Load

use "proc/ca_enroll_hh_merge", clear

* Change in Premiums, Selected t vs. Auto t+1

sort idhh year

g auto_p_pre = agecurve[_n+1] * autoprem if idhh==idhh[_n+1]
lab var auto_p_pre "Pre-Subsidy, Age-Adjusted Premium of Auto-Assigned Plan"

g auto_p_post = auto_p_pre - aptc[_n+1] if idhh==idhh[_n+1]
lab var auto_p_post "Post-Subsidy, Age-Adjusted Premium of Auto-Assigned Plan"

g autodp = auto_p_post - p_post if idhh==idhh[_n+1]
lab var autodp "Change in Post-Subsidy Premium: Auto Assigned Plan re Last Year's Plan"

count if autodp==. & idhh==idhh[_n+1] // 4855

* Move Auto Variables Forward One Year

sort  idhh year
order idhh year planid autoplanid-autodp agecurve p_post

foreach y of varlist autoplanid autonetid {
	replace `y' = `y'[_n-1] if idhh==idhh[_n-1] & year==2018
	replace `y' = `y'[_n-1] if idhh==idhh[_n-1] & year==2017
	replace `y' = `y'[_n-1] if idhh==idhh[_n-1] & year==2016
	replace `y' = `y'[_n-1] if idhh==idhh[_n-1] & year==2015
	replace `y' = "" if year==2014 | idhh!=idhh[_n-1]
}

foreach y of varlist autoprem-autodp {
	replace `y' = `y'[_n-1] if idhh==idhh[_n-1] & year==2018
	replace `y' = `y'[_n-1] if idhh==idhh[_n-1] & year==2017
	replace `y' = `y'[_n-1] if idhh==idhh[_n-1] & year==2016
	replace `y' = `y'[_n-1] if idhh==idhh[_n-1] & year==2015
	replace `y' = . if year==2014 | idhh!=idhh[_n-1]
}

* Lags

sort  idhh year
order idhh year planid

g lagprem = p_pre[_n-1]  if idhh == idhh[_n-1] // Is APTC censored if > gross prem?
g lagaptc = aptc[_n-1]   if idhh == idhh[_n-1]
g lagplan = planid[_n-1] if idhh == idhh[_n-1]
g lagnet  = netid[_n-1]  if idhh == idhh[_n-1]
g lagiss  = issid[_n-1]  if idhh == idhh[_n-1]
g lagstat = status[_n-1] if idhh == idhh[_n-1]
replace lagstat = "New" if lagstat==""

lab var lagprem "Gross Premium of Chosen Plan from Previous Year"
lab var lagaptc "APTC from Previous Year"
lab var lagplan "Plan ID from Previous Year"
lab var lagnet  "Net ID from Previous Year"
lab var lagiss  "Iss ID from Previous Year"
lab var lagstat "Enrollment Status in Previous Year"

* Inertia Categories with Automatic Assignment

foreach i in plan net { // iss

g icat`i' = . 

replace icat`i' = 1 if idhh!=idhh[_n-1]                          // new enr
replace icat`i' = 2 if idhh==idhh[_n-1] & `i'id==auto`i'id[_n-1] // stayer
replace icat`i' = 3 if idhh==idhh[_n-1] & `i'id!=auto`i'id[_n-1] // switcher
*replace icat`i' = 4 if idhh==idhh[_n-1] & year[_n-1]==2014 & year==2015 & issid[_n-1]=="99483" // Contra Costa
*replace icat`i' = 4 if idhh==idhh[_n-1] & year[_n-1]==2016 & year==2017 & issid[_n-1]=="37873" // United

lab var icat`i' "Inertia Category"
lab def icat`i' 1 "New" 2 "Stayer" 3 "Switcher" 4 "Discontinued", replace
lab val icat`i' icat`i'

tab icat`i' year

}

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
lab var netname "Network Name"

* Save

save "proc\ca_enroll_hh_wi", replace
