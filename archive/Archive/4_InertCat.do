*********************************************************************
* Combine Enrollment, Premium Data to Backout FPL and Code Inertia
*********************************************************************

/// NOTES ///

* Available Data Files
	* 1: ca_enroll_hh: Cleaned enrollment data (hh-year)
	* 2: planchars: Plan characteristics (plan-gra-year)
	* 3: inert_xwalk: Crosswalks plan_t-1 to plan_t for auto assign (plan-gra-year)
	
* Roadmap
	* Discontinued Indicator: Identify if plan's previous plan discontinued (hhdisc; hh-year)
	* 

/// SETUP ///

clear *
cd "G:\Shared Drives\CovCAInertia"

/// DISCONTINUED INDICATOR ///

* Load

use year gra planid idhh using "proc\ca_enroll_hh", clear

* Merge Plan Characteristics, Define Lags

merge m:1 year gra planid using proc/planchars, keep(master match) keepusing(netid issuerid) nogen
	sort  idhh year 
	order idhh year issuerid netid planid gra
	g planidlag = planid[_n-1]   if idhh==idhh[_n-1]
	g netidlag  = netid[_n-1]    if idhh==idhh[_n-1]
	g issidlag  = issuerid[_n-1] if idhh==idhh[_n-1]
	drop planid netid issuerid
	rename (planidlag netidlag issidlag) (planid netid issuerid)
	
* See if Lags Were Available in Year t
	
merge m:1 year gra planid using proc/planchars, keep(master match) keepusing(prem)
	sort idhh year
	g discplan = _merge==1 if idhh==idhh[_n-1] // plan of previous year does not match with current year's plan list
	drop _merge
	
merge m:1 year gra netid using proc/netlist, keep(master match)
	sort idhh year
	g discnet = _merge==1 if idhh==idhh[_n-1]
	drop _merge
	
merge m:1 year gra issuerid using proc/isslist, keep(master match)
	sort idhh year
	g disciss = _merge==1 if idhh==idhh[_n-1]
	keep year idhh disc*
	
lab var discplan "Household's Plan in Previous Year Was Discontinued"
lab var discnet  "Household's Network in Previous Year Was Discontinued"
lab var disciss  "Household's Insurer in Previous Year Was Discontinued"

save "proc/hhdisc", replace

/// ADD PLAN BENEFITS, INERTIA INDICATORS TO ENROLLMENT DATA ///

* Merge Plan Characteristics, Not Discontinued Indicator

use "proc\ca_enroll_hh", clear

merge m:1 year gra planid using proc/planchars, keep(master match) nogen
drop issuername
lab var prem "Base Premium of Selected Plan"
rename issuerid issid // check in final data set for variation

merge 1:1 year idhh using proc/hhdisc, keep(master match) nogen

* Auto Assignment Chars

rename planid planid0
merge m:1 year gra planid0 using proc/inertxwalk, keep(master match)
tab year _merge // all 2018 (correct)
tab issuername gra if _merge==1 & year==2014 // 5677 Anthem enrollees 2014 gra 19; this could be a problem
tab issuername gra if _merge==1 & year==2015 // 2 Anthem 2015, 467 SHARP 2015, 1 Western 2015; too small to matter
rename _merge inertxwalkmatch
rename planid0 planid
rename (planid1 netid1 prem1) (autoplanid autonetid autoprem)
rename (iplanchange inetchange) (autoplanchange autonetchange)
drop metal0 netid0 issuername0 issuerid0

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

* Delta Net Premium for Auto Enr

g dprem = 

START HERE!!!!

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
