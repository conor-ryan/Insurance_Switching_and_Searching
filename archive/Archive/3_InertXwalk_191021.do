*********************************************************
* Crosswalk Plans from Year to Year for Auto Assignment
*********************************************************

/// SETUP ///

cd "G:\Shared Drives\CovCAInertia"
		
/// 2017 TO 2018 ///

* 2018 Minimum Premiums within Metal and GRA

use proc/planchars18, clear
bys gra metal (prem): g minpmtl = _n==1
keep if minpmtl==1
unique gra metal // no ties
keep gra metal planid netid
rename (planid netid) (minpmtlid minpmtlnet) 
lab var minpmtlid  "Plan ID of Lowest Premium Plan by Metal, GRA"
lab var minpmtlnet "Net ID of Lowest Premium Plan by Metal, GRA"
save "proc/minpmtlid18", replace
		
* Connect 17-18 Mismatches

use proc/planchars17, clear

merge 1:1 gra planid using proc/planchars18 // Merge identifies mismatches
	replace metal = "HDHP Bronze" if metal=="Bronze HDHP" | metal=="HSA Bronze"
	order issuername gra metal year
	sort  issuername gra metal year

g inertid  = planid if _merge==3
lab var inertid  "Crosswalked Plan ID"
		
tab issuername year if _merge!=3 // Determine problem insurers; Anthem and HealthNet

tab gra _merge if issuername=="Anthem"  // GRAs 1 7 10 remained; all else exited
tab gra _merge if issuername=="Anthem" & (gra==1 | gra==7 | gra==10) // Manual matches where insurer did not exit
replace inertid = planid[_n+1] if inertid=="" & year==2017 & gra==1  & issuername=="Anthem" // only possible combo
replace inertid = planid[_n+1] if inertid=="" & year==2017 & gra==10 & issuername=="Anthem" // only possible combo
replace inertid = planid[_n+1] if inertid=="" & year==2017 & gra==7  & issuername=="Anthem" & year[_n+1]==2018 // 2017a to 2018
replace inertid = planid[_n+2] if inertid=="" & year==2017 & gra==7  & issuername=="Anthem" & year[_n+2]==2018 // 2017b to 2018

tab gra _merge if issuername=="HealthNet" // Exited GRAs 1 2 4-14; remained elsewhere, only mismatch is 3
tab gra _merge if issuername=="HealthNet" & gra==3
replace inertid = planid[_n+1] if inertid=="" & year==2017 & gra==3 & issuername=="HealthNet"

drop if year==2018
tab gra issuername if inertid==""

merge m:1 gra metal using proc/minpmtlid18, nogen      // Default to lowest prem plan w/in year-gra-metal
replace inertid = minpmtlid if inertid=="" & _merge==1 // Exited Anthem and HealthNet GRAs

keep gra issuername metal year planid netid inert*

save "proc/xwalk1718", replace
		
/// 2016 TO 2017 ///

* Connect 16-17 Mismatches

use proc/planchars16, clear

merge 1:1 gra planid using proc/planchars17 // Merge identifies mismatches
	g match = _merge==3
	replace metal = "HDHP Bronze" if metal=="Bronze HDHP" | metal=="HSA Bronze"
	order match issuername gra metal year
	sort  match issuername gra metal year
	
g inertid  = planid if _merge==3
lab var inertid  "Crosswalked Plan ID"

tab issuername _merge // Determine problem insurers; Anthem

tab gra _merge if issuername=="Anthem" // Not an issue for GRAs 10-14

replace inertid = planid[_n+1] if inertid=="" & year==2016 & issuername=="Anthem" // only possible combo; flip from PPO to EPO

drop if year==2017
tab gra issuername if inertid==""

keep gra issuername metal year planid netid inert*

save "proc/xwalk1617", replace

/// 2015 TO 2016 ///

* Connect 15-16 Matches

use proc/planchars15, clear

merge 1:1 gra planid using proc/planchars16 // Merge identifies mismatches
	g match = _merge==3
	replace metal = "HDHP Bronze" if metal=="Bronze HDHP" | metal=="HSA Bronze"
	order issuername gra metal year
	sort  issuername gra metal year
	
g inertid  = planid if _merge==3
lab var inertid  "Crosswalked Plan ID"

tab issuername _merge // Determine problem insurers; Blue Shield, HealthNet, Sharp

tab gra _merge if issuername=="Blue Shield"
replace inertid = planid[_n+1] if inertid=="" & year==2015 & issuername=="Blue Shield" // only possible combo; EPO to PPO

tab gra _merge if issuername=="HealthNet" // 2 4-10 14
replace inertid = planid[_n+1] if inertid=="" & year==2015 & issuername=="HealthNet" // only possible combo; netid already consistent

tab metal _merge if issuername=="Sharp" // Metals don't match; nobody must have enrolled

drop if year==2016
tab gra issuername if inertid==""

keep gra issuername metal year planid netid inert*

save "proc/xwalk1516", replace

/// 2014 TO 2015 ///

* Connect 14-15 Matches

use proc/planchars14, clear

rename plantype plantype0

merge 1:1 gra planid using proc/planchars15 // Merge identifies mismatches
	g match = _merge==3
	replace metal = "HDHP Bronze" if metal=="Bronze HDHP" | metal=="HSA Bronze"
	rename plantype plantype1
	order issuername gra metal year
	sort  issuername gra metal year
	
g inertid = planid if _merge==3
lab var inertid  "Crosswalked Plan ID"

tab issuername _merge // Determine problem insurers; Anthem, HealthNet, Sharp, Western
tab gra   _merge    if issuername=="Anthem" 
tab metal _merge    if issuername=="Anthem" // More than one possible combo; can isolate by plan type
replace inertid = planid[_n+1] if inertid=="" & year==2014 & year[_n+1]==2015 & plantype0==plantype1[_n+1] & issuername=="Anthem"
replace inertid = planid[_n+2] if inertid=="" & year==2014 & year[_n+2]==2015 & plantype0==plantype1[_n+2] & issuername=="Anthem"
replace inertid = planid[_n+3] if inertid=="" & year==2014 & year[_n+3]==2015 & plantype0==plantype1[_n+3] & issuername=="Anthem"
count if inertid=="" & _merge==1 & issuername=="Anthem" // 9
replace inertid = planid[_n+1] if inertid=="" & year==2014 & year[_n+1]==2015 & issuername=="Anthem"
replace inertid = planid[_n+2] if inertid=="" & year==2014 & year[_n+2]==2015 & issuername=="Anthem"
replace inertid = planid[_n+3] if inertid=="" & year==2014 & year[_n+3]==2015 & issuername=="Anthem"
count if inertid=="" & _merge==1 & issuername=="Anthem" // 0; fixed bronze HSAs

tab gra   _merge if issuername=="HealthNet" 
tab metal _merge if issuername=="HealthNet" 
replace inertid = planid[_n+1] if inertid=="" & year==2014 & year[_n+1]==2015 & issuername=="HealthNet" // only possible combo
count if inertid=="" & _merge==1 & issuername=="HealthNet" 

tab metal _merge if issuername=="Sharp" // Metals don't match; nobody must have enrolled
count if inertid=="" & _merge==1 & issuername=="Sharp" 

tab metal _merge if issuername=="Western"
replace inertid = planid[_n+1] if inertid=="" & year==2014 & year[_n+1]==2015 & issuername=="Western" // only possible combo
count if inertid=="" & _merge==1 & issuername=="Western"

drop if year==2015
tab gra issuername if inertid==""

keep gra issuername metal year planid netid inert*

save "proc/xwalk1415", replace

/// CONSOLIDATE /// NEED TO NOT CODE NETS AS SWITCHING WHEN IT'S JUST A CHANGE OF CODE (EG ANTHEM 2017)

* Add Networks from Next Year's Plan

forvalues x = 14/17 {

local y = `x' + 1

use proc/xwalk`x'`y', clear
rename (planid inertid netid) (planid0 planid netid0)

merge m:1 gra planid using proc/planchars`y', keep(master match) 
	// 14 master: contra costa
	// 15 master: 1 sharp plan
	// 16 master: united
	// 17 master: none
	if `x'==15 drop if _merge==1 // could fix but just isn't going to matter

drop prodid plantype _merge
compress
order year gra issuername issuerid metal planid0 planid netid0 netid prem
rename (issuername issuerid metal) (issuername0 issuerid0 metal0)
rename (planid netid prem) (planid1 netid1 prem1)

lab var metal   "Metal Prev Plan"
lab var planid0 "Plan ID Prev Plan"
lab var netid0  "Net ID Prev Plan"
lab var netid1  "Crosswalked Net ID"
lab var prem    "Crosswalked Plan Prem (Unadj)"

save proc/xwalk`x'`y'net, replace

}

* Append

use proc/xwalk1718net, clear
append using proc/xwalk1617net
append using proc/xwalk1516net
append using proc/xwalk1415net

* Change Indicators

g iplanchange = planid0!=planid1
replace iplanchange = . if planid1==""
lab var iplanchange "Crosswalked to New Plan"

g inetchange  = netid0!=netid1
replace inetchange = . if planid1==""
lab var inetchange "Crosswalked to New Network"

assert inetchange!=1 if iplanchange==0

tab year iplanchange
tab year inetchange // pretty correlated...

tab year issuername if inetchange==1 // shows issuer-year where net changes happened

* Save

save proc/inertxwalk, replace // Need to add premiums of previous plan? Determine by proceeding w code

* LA Spot Check

use proc/inertxwalk if gra==16, clear
sort year gra issuername0 netid0 metal0

/// NOTES ///

/* Setup
	
	Question: Does HH previously enrolled in plan x (or ~x) have inertia wrt plan x 
		HH always has inertia for plan x if previously enrolled in plan x
		HH may have inertia wrt plan x if previously enrolled in plan ~x
		HH may have inertia wrt plan ~x if previously enrolled in plan x, but NOT if plan x is still offered
	
	Setup
		HH will have field for previous plan from last year
		Plan will have field for plan(s) from previous year that will be inertial to it (could be none)
		If HH matches one of those fields, HH has inertia wrt plan
		
	Scenarios
		- Plan x is same in 17, 18: Matched
		- Plan x and its insurer leave from 17 to 18: Match to lowest within-metal plan
		- Plan x leaves from 17 to 18 but its insurer doesn't: Match to lowest within-metal-ins plan
		- Plan x changes IDs from 17 to 18: Match to lowest within-metal-ins-type plan
		
	Merge 18 to 17
		1: 18 plan that does not match to 17 plan
		2: 17 plan that does not match to 18 plan
		3: 18 plan matches to 17 plan
		
	Goal: Identify inertia matches for 18 plans
		Identify mismatches between 18, 17
		Every 2017 plan (merge=2) should have a match for 2018 
		Prior to 2017, every plan should have a match if the insurer didn't exit 
		
	Dealing with Anthem, HealthNet Exits
		Identify Plan ID of minimum prem plan by metal, GRA in 2018
		Get list of Anthem plans that should crosswalk to such plans in 2018
		Set inertids for min plans as list of eligible Anthem plans
		
*/
