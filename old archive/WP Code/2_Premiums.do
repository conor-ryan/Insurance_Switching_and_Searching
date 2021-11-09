********************************************
* Format CCA Premium Data
********************************************

/// SETUP ///

cd "\\files.umn.edu\cla\Home\ryan0463\My Documents\CovCAInertia"

/// NETWORK CROSSWALKS ///

* Networks 2016-2018

forvalues x = 16/18 {
import delimited using "data\pufcaplan20`x'.csv", clear
keep standardcomponentid hiosproductid networkid plantype
rename (standardcomponentid hiosproductid networkid) (planid prodid netidorig)
g issid = substr(planid,1,5)
g netid = substr(prodid,1,5) + substr(netid,-3,3)
order netid prodid planid
keep  netid prodid planid
sort  netid prodid planid
duplicates drop planid, force
keep netid prodid 
duplicates drop prodid, force
replace prodid = subinstr(prodid," ","",.)
replace netid  = subinstr(prodid," ","",.)
save "proc/netxwalk`x'", replace
}

/// PLAN-LEVEL PREMIUMS 2014 ///

* Import and Format 2014 Premiums Spreadsheet

import excel using "data\prem2014wID.xlsx", first case(lower) cellrange(A8:E504) clear sheet("2014 Individual Rates")

rename (producttype ratingregion healthinsurancecompany monthlyprice metaltier) (plantype gra issuername prem metal)
g year = 2014
keep  year gra issuername prem metal plantype
order year gra issuername prem metal plantype

replace issuername="Anthem"       if issuername=="Anthem Blue Cross"
replace issuername="Blue Shield"  if issuername=="Blue Shield of California"
replace issuername="Chinese"      if issuername=="Chinese Community Health Plan"
replace issuername="Contra Costa" if issuername=="Contra Costa Health Plan"
replace issuername="HealthNet"    if issuername=="Health Net"
replace issuername="Kaiser"       if issuername=="Kaiser Permanente"
replace issuername="Western"      if issuername=="Western Health Advantage"

bys gra issuername metal plantype: g dups  = _N
bys gra issuername metal plantype: g dupno = _n 
sort gra issuername metal plantype prem
replace metal = "Bronze HDHP" if dups==2 & dupno==1 & metal=="Bronze"
drop dups dupno

duplicates drop gra issuername metal plantype, force // 3 Sharp plans

replace prem = prem / 1.278 // Data reflects plans for 40 y/o single; remove age curve

save "proc\prem14", replace

* Load Enrollment Data

use "proc\ca_enroll_hh" if year==2014, clear
rename (plan_name issuer_name) (planname issuername)

* Change to 14-Digit Plan ID

duplicates drop planid gra, force
keep year gra planname issuername planid

* Issuer

tab issuername
replace issuername = "Anthem"       if issuername=="Anthem Blue Cross"
replace issuername = "Chinese"      if issuername=="Chinese Community"
replace issuername = "Contra Costa" if issuername=="Contra Costa Heal"
replace issuername = "HealthNet"    if issuername=="Health Net"
replace issuername = "Molina"       if issuername=="Molina Health Car"
replace issuername = "Sharp"        if issuername=="SHARP Health Plan"
replace issuername = "Valley"       if issuername=="Valley Health"
replace issuername = "Western"      if issuername=="Western Health"

* Metal

g metal = ""
replace metal = "Catastrophic" if regexm(planname,"Minimu")==1
replace metal = "Bronze"       if regexm(planname,"Bronze")==1
replace metal = "Bronze HDHP"  if regexm(planname,"HSA")==1
replace metal = "Silver"       if regexm(planname,"Silver")==1
replace metal = "Gold"         if regexm(planname,"Gold")==1
replace metal = "Platinum"     if regexm(planname,"Platin")==1
drop if regexm(planname,"AI-AN")==1

* Plan Type

g plantype = ""
replace plantype = "HMO" if regexm(planname,"HMO")==1
replace plantype = "PPO" if regexm(planname,"PPO")==1
replace plantype = "EPO" if regexm(planname,"EPO")==1

tab planname if plantype==""
replace plantype = "HMO" if issuername=="Contra Costa"
replace plantype = "HMO" if issuername=="Kaiser"
replace plantype = "HMO" if issuername=="Molina"
replace plantype = "HMO" if issuername=="Sharp"
replace plantype = "HMO" if issuername=="Valley"
replace plantype = "HMO" if issuername=="Western"

tab planid if plantype==""
replace plantype = "PPO" if regexm(planid,"27603CA114")==1
replace plantype = "EPO" if regexm(planid,"27603CA116")==1
replace plantype = "PPO" if regexm(planid,"70285CA135")==1
replace plantype = "EPO" if regexm(planid,"70285CA136")==1

* Merge Premiums

merge m:1 gra issuername metal plantype using "proc\prem14", keep(master match) nogen
replace prem = 233.61 if prem==. // Missing bronze hdhp prem set to equivalent bronze prem

* Merge Networks

g prodid = substr(planid,1,10)
merge m:1 prodid using proc/netxwalk16, keep(master match) nogen
replace netid = "70285002" if netid=="" & issuername=="Blue Shield"  & plantype=="EPO" // Blue EPO was not continued in 2016
replace netid = "99483001" if netid=="" & issuername=="Contra Costa" & plantype=="HMO" // Contra Costa only existed with one network in 2014
replace netid = "99910002" if netid=="" & issuername=="HealthNet"    & plantype=="PPO" // This network only existed in 2014

* Save

order year gra planid issuername planname prem metal plantype
save proc/planchars14, replace

/// PREMIUMS 2015-2018 ///

* 2015

import excel using data/prem2015wID.xlsx, first case(lower) clear 
rename (ratingareaid individualrate applicant hiosid metallevel fullplanname) (gra prem issuername issuerid metal planname)
keep if age=="21"
g year = 2015
g plantype = substr(indexname,-3,3)
destring gra, ignore("Rating Area ") replace
keep  year gra issuerid planid issuername prem metal plantype
order year gra issuerid planid issuername prem metal plantype
g prodid = substr(planid,1,10)
merge m:1 prodid using proc/netxwalk16, keep(master match) nogen
replace netid = "70285002" if netid=="" & issuername=="Blue Shield" & plantype=="EPO" // Blue EPO was not continued in 2016
replace issuername = "Chinese" if issuername=="Chinese C."
replace issuername = "LA Care" if issuername=="L.A. Care"
save proc/planchars15, replace // year gra issuerid planid issuername prem metal plantype prodid netid

* 2016

import excel using data\prem2016.xlsx, first case(lower) cellrange(A8:K25354) clear
rename (ratingareaid individualrate applicant hiosid metallevel network product fullplanname) ///
(gra prem issuername issuerid metal plantype product planname)
keep if age=="21"
g year = 2016
destring gra, ignore("Rating Area ") replace
keep  year gra issuerid planid issuername prem metal plantype
order year gra issuerid planid issuername prem metal plantype
g prodid = substr(planid,1,10)
merge m:1 prodid using proc/netxwalk16, keep(master match) nogen
replace issuername = "Chinese" if issuername=="Chinese C."
replace issuername = "LA Care" if issuername=="L.A. Care"
save proc/planchars16, replace

* 2017

import excel using data\prem2017.xlsx, first case(lower) cellrange(A7:K27929) clear
rename (ratingareaid individualrate applicant hoisid metallevel network product fullplanname) ///
(gra prem issuername issuerid metal plantype product planname)
keep if age=="21"
g year = 2017
destring gra, ignore("Rating Area ") replace
keep  year gra issuerid planid issuername prem metal plantype
order year gra issuerid planid issuername prem metal plantype
g prodid = substr(planid,1,10)
merge m:1 prodid using proc/netxwalk17, keep(master match) nogen
replace issuername = "Chinese" if issuername=="Chinese Community"
replace issuername = "LA Care" if issuername=="L.A. Care"
save proc/planchars17, replace

* 2018

import excel using data\prem2018.xlsx, first case(lower) cellrange(A8:K25916) clear
rename (ratingareaid individualrate applicant hoisid metallevel network product fullplanname) ///
(gra prem issuername issuerid metal plantype product planname)
keep if age=="21"
g year = 2018
destring gra, ignore("Rating Area ") replace
keep  year gra issuerid planid issuername prem metal plantype
order year gra issuerid planid issuername prem metal plantype
g prodid = substr(planid,1,10)
merge m:1 prodid using proc/netxwalk18, keep(master match) nogen // GET THIS TO WORK W 2018
replace issuername = "Chinese"   if issuername=="Chinese C."
replace issuername = "HealthNet" if regexm(issuername,"HealthNet")==1
replace issuername = "LA Care"   if issuername=="L.A. Care"
tostring issuerid, replace
save proc/planchars18, replace

/// CREATE COMBINED FILE WITH LAG OFFER INDICATOR ///

* Append Years

use proc/planchars18, clear
append using proc/planchars17
append using proc/planchars16
append using proc/planchars15
append using proc/planchars14

sort  year gra issuerid plantype metal planid
order year gra issuerid issuername plantype metal planid prem

* Formatting

replace issuerid = substr(planid,1,5) if issuerid==""

replace metal = "HDHP Bronze" if metal=="Bronze HDHP" | metal=="HSA Bronze"

replace plantype = subinstr(plantype,"HDHP ","",.)
replace plantype = subinstr(plantype,"HSA ","",.)
replace plantype = "EPO" if regexm(plantype,"EPO")==1
replace plantype = "HMO" if regexm(plantype,"HMO")==1
replace plantype = "PPO" if regexm(plantype,"PPO")==1
replace plantype = "HSP" if plantype=="HCSP"

sort  year gra metal prem
bys year gra metal: g pmtlord = _n
replace pmtlord = 5 if pmtlord>5
order year gra metal prem pmtlord
lab var pmtlord "Premium Order (Metal-GRA-Year)"

* Clean and Save

drop prodid planname
order year gra issuerid netid metal planid
sort  year gra issuerid netid metal planid

lab var netid "Network ID"
lab var year "Year"

save proc/planchars, replace

/// LISTS ///

* Lagged Premiums 

use proc/planchars, clear
sort planid gra year
g lagpremlag = prem[_n-1] if planid==planid[_n-1] & gra==gra[_n-1] & year-1==year[_n-1]
lab var lagpremlag "Previous Plan Premium in Previous Year"
keep year gra planid prem lagpremlag
rename (planid prem) (lagplan lagprem)
lab var lagprem "Previous Plan Premium"
save proc/lagplanprem, replace

* Network List

use proc/planchars, clear
duplicates drop year gra netid, force
keep year gra netid
save proc/netlist, replace

* Issuer List

use proc/planchars, clear
duplicates drop year gra issuerid, force
keep year gra issuerid
save proc/isslist, replace
