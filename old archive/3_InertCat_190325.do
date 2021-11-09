*********************************************************************
* Combine Enrollment, Premium Data to Backout FPL and Code Inertia
*********************************************************************

/// SETUP ///

clear *
cd "\\files.umn.edu\cla\Home\ryan0463\My Documents\CovCAInertia"

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
merge 1:1 year idhh       using proc/hhdisc,    keep(master match) nogen

sort  idhh year
order idhh year planid

****** there may be a bug with issuerid/issid variables. final data set has both, looks like only one varies with plan choice;
rename issuerid issid

* Lags
egen household = group(idhh)
tsset household year

g lagplan = planid[_n-1] if household == L.household
g lagnet  = netid[_n-1]  if household == L.household
g lagiss  = issid[_n-1]  if household == L.household
g lagstat = status[_n-1] if household == L.household
replace lagstat = "New" if lagstat==""

lab var lagplan "Plan ID from Previous Year"
lab var lagnet  "Net ID from Previous Year"
lab var lagiss  "Iss ID from Previous Year"
lab var lagstat "Enrollment Status in Previous Year"

* Inertia Categories

foreach i in plan net iss {

g icat`i' = . // Use status variable? Can't have inertia if terminated in previous year

replace icat`i' = 1 if idhh!=idhh[_n-1]                                   // new enrollee
replace icat`i' = 2 if idhh==idhh[_n-1] & `i'id==`i'id[_n-1]              // stayer
replace icat`i' = 3 if idhh==idhh[_n-1] & `i'id!=`i'id[_n-1] & disc`i'==0 // switcher
replace icat`i' = 4 if idhh==idhh[_n-1] & `i'id!=`i'id[_n-1] & disc`i'==1 // discontinued

lab var icat`i' "Inertia Category"
lab def icat`i' 1 "New" 2 "Stayer" 3 "Switcher" 4 "Discontinued"
lab val icat`i' icat`i'

tab icat`i' year

}

* Decision Categories

g dcat = lagstat=="Enrolled" & icatplan!=4
replace dcat = 2 if lagstat=="Terminated" & icatplan!=4
lab var dcat "Decision Category"
lab def dcat 0 "New Enrollment" 1 "Continuous Enrollment" 2 "Discontinuous Enrollment"
lab val dcat dcat

* Previous Plan Premium and Lag Premium

merge m:1 year gra lagplan using proc/lagplanprem, keep(master match) nogen
g dprem = lagprem - lagpremlag
egen dpremmed  = median(dprem), by(year gra)
g pshock = .
replace pshock = 0 if dprem <= dpremmed & dprem!=. & dpremmed!=.
replace pshock = 1 if dprem >  dpremmed & dprem!=. & dpremmed!=.
lab var dprem    "Change in Own-Premium of Previous Plan"
lab var dpremmed "Median Change in Own-Premium of Previous Plan"
lab var pshock   "Premium Shock (Premium Change > Median Premium Change)"

* Experience

sort idhh year
g iexp = 0
forvalues y = 1/4 {
replace iexp = `y' if idhh==idhh[_n-`y']
}
lab var iexp "Inertia: Experience" 

* Tenure

foreach i in plan net iss {

g iten`i' = 0
replace iten`i' = 1 if iexp==1
replace iten`i' = 2 if iexp==2 & `i'id[_n-1]==`i'id[_n-2]
replace iten`i' = 3 if iexp==3 & `i'id[_n-1]==`i'id[_n-2] & `i'id[_n-2]==`i'id[_n-3]
replace iten`i' = 4 if iexp==4 & `i'id[_n-1]==`i'id[_n-2] & `i'id[_n-2]==`i'id[_n-3] & `i'id[_n-3]==`i'id[_n-4]
lab var iten`i' "Inertia: `i' Tenure"

}

* Correct Issuer Names

replace issuername = "Anthem" if issuername==""    & issuer_name=="Anthem Blue Cross"
replace issuername = "HealthNet" if issuername=="" & issuer_name=="Health Net"
replace issuername = "Sharp" if issuername==""     & issuer_name=="SHARP Health Plan"
replace issuername = "Western" if issuername==""   & issuer_name=="Western Health"

* Network Name

g netname = issuername + " " + plantype

* Save

save "proc\ca_enroll_hh_wi", replace // CHECK TABS OF DCAT, PSHOCK, LAGPREM, LAGPREMLAG VARIABLES
