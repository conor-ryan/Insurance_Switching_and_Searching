*********************************************************
* Crosswalk Plans from Year to Year for Auto Assignment
*********************************************************
		
/// 2017 TO 2018 ///

* 2018 Minimum Premiums within Metal and GRA

use "proc/planchars18", clear
bys gra metal (prem): g minpmtl = _n==1
keep if minpmtl==1
unique gra metal // no ties
keep gra metal planid netid
rename (planid netid) (minpmtlid minpmtlnet) 
lab var minpmtlid  "Plan ID of Lowest Premium Plan by Metal, GRA"
lab var minpmtlnet "Net ID of Lowest Premium Plan by Metal, GRA"
save "proc/minpmtlid18", replace
		
* Connect 17-18 Mismatches
use proc/planchars if year==2018, clear
rename planid planidx
save "proc/planmerge18", replace

use proc/planchars if year==2017, clear
replace year = year+1

merge 1:1 year gra issuername plantype product metal using proc/planmerge18 // Merge identifies mismatches
	drop if _merge==2 & year!=2018
	replace year = year-1 if _merge!=2
	order issuername gra metal year
	sort  issuername gra metal year

g inertid = planidx if _merge==3
lab var inertid  "Crosswalked Plan ID"
		
tab issuername year if _merge!=3 // Determine problem insurers; Anthem and HealthNet

tab gra _merge if issuername=="Anthem"  // GRAs 1 7 10 remained; all else exited 
tab gra _merge if issuername=="Anthem" & (gra==1 | gra==7 | gra==10) // Manual matches where insurer did not exit
*replace inertid = planidx[_n+1] if inertid=="" & year==2017 & gra==1  & issuername=="Anthem" // only possible combo
replace inertid = planidx[_n+1] if inertid=="" & year==2017 & gra==10 & issuername=="Anthem" // only possible combo
*replace inertid = planidx[_n+1] if inertid=="" & year==2017 & gra==7  & issuername=="Anthem" & year[_n+1]==2018 // 2017a to 2018
*replace inertid = planidx[_n+2] if inertid=="" & year==2017 & gra==7  & issuername=="Anthem" & year[_n+2]==2018 // 2017b to 2018
assert inertid!="" if issuername=="Anthem" & _merge!=2 & (gra==1|gra==10) // HMO plans in 7 Exit? 

tab gra _merge if issuername=="HealthNet" // Only mismatch is on 3
tab gra _merge if issuername=="HealthNet" & gra==3
replace inertid = planidx[_n+1] if inertid=="" & year==2017 & gra==3 & issuername=="HealthNet"
assert inertid!="" if issuername=="HealthNet" & _merge!=2 & gra==3

drop if year==2018
tab gra issuername if inertid==""

merge m:1 gra metal using proc/minpmtlid18, nogen      // Default to lowest prem plan w/in year-gra-metal
replace inertid = minpmtlid if inertid=="" & _merge==1 // 132; Exited Anthem and HealthNet GRAs

keep gra issuername metal year planid netname product inert*

save "proc/xwalk1718", replace
		
/// 2016 TO 2017 ///

* Connect 16-17 Mismatches
use proc/planchars if year==2017, clear
rename planid planidx
save "proc/planmerge17", replace

use proc/planchars if year==2016, clear
replace year = year+1

merge 1:1 year gra issuername plantype product metal using proc/planmerge17 // Merge identifies mismatches
	drop if _merge==2 & year!=2017
	replace year = year-1 if _merge!=2
	g match = _merge==3
	order match issuername gra metal year
	sort  match issuername gra metal year
	
g inertid  = planidx if _merge==3
lab var inertid  "Crosswalked Plan ID"

tab issuername _merge // Determine problem insurers; Anthem

tab gra _merge if issuername=="Anthem" // Not an issue for GRAs 10-14

replace inertid = planidx[_n+1] if inertid=="" & year==2016 & issuername=="Anthem" // only possible combo; flip from PPO to EPO
assert inertid!="" if issuername=="Anthem" & _merge!=2 & (gra<10 | gra>14)

drop if year==2017
tab gra issuername if inertid==""

keep gra issuername metal year planid netname product inert*

save "proc/xwalk1617", replace

/// 2015 TO 2016 ///

* Connect 15-16 Matches
use proc/planchars if year==2016, clear
rename planid planidx
save "proc/planmerge16", replace

use proc/planchars if year==2015, clear
replace year = year+1

merge 1:1 year gra issuername plantype product metal using proc/planmerge16 // Merge identifies mismatches
	drop if _merge==2 & year!=2016
	replace year = year-1 if _merge!=2
	g match = _merge==3
	order issuername gra metal year
	sort  issuername gra metal year
	
g inertid  = planidx if _merge==3
lab var inertid  "Crosswalked Plan ID"

tab issuername _merge // Determine problem insurers; Blue Shield, HealthNet, Sharp

tab gra _merge if issuername=="Anthem" //Anthem EPO/PPO Switch in some markets
replace inertid = planidx[_n+1] if inertid=="" & year==2015 & issuername=="Anthem" // only possible combo; EPO to PPO

tab gra _merge if issuername=="Blue Shield"
replace inertid = planidx[_n+1] if inertid=="" & year==2015 & issuername=="Blue Shield" // only possible combo; EPO to PPO

tab gra _merge if issuername=="HealthNet" // 2 4-10 14
replace inertid = planidx[_n+1] if inertid=="" & year==2015 & issuername=="HealthNet" // only possible combo; netid already consistent

tab metal _merge if issuername=="Sharp" // Metals don't match; nobody must have enrolled 
// *Pretty sure this is an error in Sharp Plan IDs
// replace inertid = "92499CA0020008"  if planid== "92499CA0020007" 
// replace inertid = "92499CA0020010"  if planid== "92499CA0020009"
// replace inertid = "92499CA0020009"  if planid== "92499CA0020008" 

drop if year==2016
tab gra issuername if inertid==""

keep gra issuername metal year planid netname product inert*

save "proc/xwalk1516", replace

/// 2014 TO 2015 ///
use proc/planchars if year==2015, clear
rename planid planidx
save "proc/planmerge15", replace

* Connect 14-15 Matches

use proc/planchars if year==2014, clear
replace year = year+1

*rename plantype plantype0
* Many Merge to account for duplicate Anthem Silver in 2014
merge m:1 year gra issuername plantype product metal using proc/planmerge15 // Merge identifies mismatches
	drop if _merge==2 & year!=2015
	replace year = year-1 if _merge!=2
	g match = _merge==3
	rename plantype plantype1
	order issuername gra metal year
	sort  issuername gra metal year
	
g inertid = planidx if _merge==3
lab var inertid  "Crosswalked Plan ID"

tab issuername _merge // Determine problem insurers; Anthem, HealthNet, Sharp, Western.

tab gra   _merge    if issuername=="Anthem" 
tab metal _merge    if issuername=="Anthem" // More than one possible combo; can isolate by plan type
* Just the HDHP Bronze plans are problems
replace inertid = planidx[_n+1] if inertid=="" & year==2014 & year[_n+1]==2015 & issuername=="Anthem"
// replace inertid = planidx[_n+2] if inertid=="" & year==2014 & year[_n+2]==2015 & plantype0==plantype1[_n+2] & issuername=="Anthem"
// replace inertid = planidx[_n+3] if inertid=="" & year==2014 & year[_n+3]==2015 & plantype0==plantype1[_n+3] & issuername=="Anthem"
// count if inertid=="" & _merge==1 & issuername=="Anthem" // 9
// replace inertid = planidx[_n+1] if inertid=="" & year==2014 & year[_n+1]==2015 & issuername=="Anthem"
// replace inertid = planidx[_n+2] if inertid=="" & year==2014 & year[_n+2]==2015 & issuername=="Anthem"
// replace inertid = planidx[_n+3] if inertid=="" & year==2014 & year[_n+3]==2015 & issuername=="Anthem"
count if inertid=="" & _merge==1 & issuername=="Anthem" // 0; fixed bronze HSAs




tab gra   _merge if issuername=="HealthNet" 
tab metal _merge if issuername=="HealthNet" 
replace inertid = planidx[_n+1] if inertid=="" & year==2014 & year[_n+1]==2015 & issuername=="HealthNet" // only possible combo
count if inertid=="" & _merge==1 & issuername=="HealthNet" 

*tab metal _merge if issuername=="Sharp" // Metals don't match; nobody must have enrolled (CBR - Adjusted Plan IDs)

// *Pretty sure this is an error in Sharp Plan IDs
// replace inertid = "92499CA0020007"  if planid== "92499CA0020008" 
// replace inertid = "92499CA0020009" if planid== "92499CA0020010" 
// replace inertid = "92499CA0020008" if planid== "92499CA0020009"
count if inertid=="" & _merge==1 & issuername=="Sharp" 

// tab metal _merge if issuername=="Western"
// replace inertid = planid[_n+1] if inertid=="" & year==2014 & year[_n+1]==2015 & issuername=="Western" // only possible combo
count if inertid=="" & _merge==1 & issuername=="Western"

drop if year==2015
tab gra issuername if inertid==""

keep gra issuername metal year planid netname product inert*

save "proc/xwalk1415", replace

/// CONSOLIDATE /// NEED TO NOT CODE NETS AS SWITCHING WHEN IT'S JUST A CHANGE OF CODE (EG ANTHEM 2017)

* Add Networks from Next Year's Plan

forvalues x = 14/17 {

local y = `x' + 1

use proc/xwalk`x'`y', clear
replace year = year+1
rename (planid netname issuername metal product inertid) (planid0 netname0 issuername0 metal0 product0 planid)


merge m:1 year gra planid using proc/planchars, keep(master match)
	replace year = year-1
	tab issuername gra if _merge==1
	// 14 master: contra costa, 1 sharp plan
	// 15 master: 1 sharp plan
	// 16 master: united
	// 17 master: none (auto assign covers exits)

keep year gra planid* netname* issuername* metal* product* prem mean*
rename (issuername netname metal product planid prem) (auto_iss auto_net auto_mtl auto_prod auto_planid auto_prem)
rename (issuername0 netname0 metal0 product0 planid0) (issuername netname metal product planid)
order year gra issuername netname planid auto_*
	
lab var auto_planid "Crosswalked Plan ID"
lab var auto_prem   "Crosswalked Plan's Base Premium"
lab var auto_iss    "Crosswalked Issuer Name"
lab var auto_net    "Crosswalked Net Name"
lab var auto_mtl    "Crosswalked Metal Level"
lab var auto_prod    "Crosswalked Product Line"

rename (mean_prem_all mean_prem_mtl) (auto_mean_prem_all auto_mean_prem_mtl)
lab var auto_mean_prem_all    "Mean Base Premiums, Excl. Auto Plan"
lab var auto_mean_prem_mtl    "Mean Base Premiums w/in Mtl Level, Excl. Auto Plan"

// Merge in Prior Year Prems for aggregate change stats
merge m:1 year gra planid using proc/planchars, keep(master match)
	tab issuername gra if _merge==1

drop issuerid netid pmtlord plantype bp _merge mean*
	
// Compute Mean Changes of Premiums, excluding each plan.
g base_dprem = auto_prem-prem
g count = 1

egen sum_dprem_all = sum(base_dprem), by(year gra)
egen sum_dprem_mtl = sum(base_dprem), by(year gra metal)

egen count_dprem_all = sum(count), by(year gra)
egen count_dprem_mtl = sum(count), by(year gra metal)

g mean_dprem_all = (sum_dprem_all - base_dprem)/(count_dprem_all - 1)
g mean_dprem_mtl = (sum_dprem_mtl - base_dprem)/(count_dprem_mtl - 1)

lab var base_dprem 		 "Change in Base Premium"
lab var mean_dprem_all   "Mean Change in All Rival Base Premiums"
lab var mean_dprem_mtl   "Mean Change in Rival Base Premiums, Same Metal"

drop count* sum* prem
	
	
save proc/xwalk`x'`y'net, replace

}

* Append

use proc/xwalk1718net, clear
append using proc/xwalk1617net
append using proc/xwalk1516net
append using proc/xwalk1415net
compress

* Change Indicators

*g iplanchange = planid!=auto_planid
g iplanchange = !(auto_net==netname & auto_mtl==metal & auto_prod==product)
replace iplanchange = . if auto_planid==""
lab var iplanchange "Crosswalked to New Plan (Only for Use in Inattention Stage)"

g inetchange = netname!=auto_net
replace inetchange = . if auto_net==""
lab var inetchange "Crosswalked to New Network (Only for Use in Inattention Stage)"

g iisschange = issuername!=auto_iss
replace iisschange = . if auto_iss==""
lab var iisschange "Crosswalked to New Issuer (Only for Use in Inattention Stage)"

count if inetchange==1 & iplanchange==0 // 54 Network name changes

tab year iplanchange
tab year inetchange // pretty correlated...

tab year issuername if iplanchange==1 // shows issuer-year where plan changes happened
tab year issuername if inetchange==1 // shows issuer-year where net changes happened


drop metal product

* Prep For Merge and Save

drop issuername netname
rename planid planid0

save proc/inertxwalk, replace // Need to add premiums of previous plan? Determine by proceeding w code





