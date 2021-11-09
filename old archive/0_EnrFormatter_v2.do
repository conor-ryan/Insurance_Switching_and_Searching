********************************************
* CovCA Formatter - Individual Data
********************************************

/// INDIVIDUAL ENROLLMENT ///

* Load

clear *
import delimited using "D:\CC PRA enrollment file_022018.csv", clear

* Rename and Label

rename (ahbx_case_id_x indv_id_x enrlee_enrlmnt_yr enrlee_status rating_area) (idhh idind year status gra)
rename (hios_id_16 active_passive_d subsidy_fpl_percent_round gross_premium_amt_round net_premium_amt_round aptc_amt_round) (planid act_pass fpl p_pre p_post aptc)
rename (zip_3 enrlee_plan_select_dt_wk cov_start_dt_wk cov_end_dt_wk) (zip3 enr_strt cov_strt cov_end)

lab var act_pass    "Acctive or Passive Enrollment"
lab var age         "Age"
lab var aptc        "Advanced Premium Tax Credit"
lab var cov_strt    "Enrollee Coverage Start Week"
lab var cov_end     "Enrollee Coverage End Week"
lab var enr_strt    "Enrollee Enrollment Week"
lab var fpl         "FPL"
lab var gender      "Gender"
lab var gra         "Rating Area"
lab var idhh        "Household ID"
lab var idind       "Individual ID"
lab var issuer_name "Issuer Name"
lab var p_post      "Post-Subsidy Premium"
lab var p_pre       "Pre-Subsidy Premium"
lab var plan_mrkt   "Plan Market"
lab var plan_name   "Plan Name"
lab var plan_type   "Plan Type"
lab var planid      "Plan ID: 16"
lab var region      "Rating Area"
lab var status      "Enrollment Status"
lab var year        "Year"
lab var zip3        "3-Digit Zip Code"

* IDs
	
replace idhh  = subinstr(idhh,"CASE_","",.)
replace idind = subinstr(idind,"INDV_","",.)

g idhh_miss  = idhh==""
g idind_miss = idin==""
lab var idhh_miss  "Missing IDHH"
lab var idind_miss "Missing IDIND"

gen id2 = _n
tostring id2, replace
replace idhh  = "ZZHH"  + id2 if idhh==""
replace idind = "ZZIND" + id2 if idind==""

drop id2

* Age

replace age = 100 if age > 100 // 329

* Merge Age Curve
		
merge m:1 age using "data\agecurve", keep(master match) nogen // 11 missing age

* APTC

replace aptc = 0 if p_pre - p_post==0 & aptc==.

* Subsidy
	
g hassub = aptc>0
lab var hassub "Enrollee has subsidy"
lab def hassub 0 "Unsubsidized" 1 "Subsidized"
lab val hassub hassub

* FPL 

tab year if fpl == . 
tab year if fpl != . & fpl < 401 & aptc==0

g fplcat = (fpl < 133 & fpl != . & fpl != 0)
replace fplcat = 2 if fpl >= 133 & fpl <= 250
replace fplcat = 3 if fpl >  250 & fpl <= 400
replace fplcat = 4 if fpl >  400
replace fplcat = 5 if fpl >= 133 & fpl <= 400 & aptc == 0
replace fplcat = 6 if fpl == 0 | fpl == .

lab var fplcat "FPL Category"
lab def fplcat 1 "<133%" 2 "133-250%" 3 "251-400%" 4 "400%+" 5 "133-400% No APTC" 6 "Unknown"
lab val fplcat fplcat

* Save

drop region plan_mrkt 
compress
save "D:\proc\ca_enroll_ind", replace

/// HOUSEHOLD ENROLLMENT ///

* Load

use "D:\proc\ca_enroll_ind", clear

* Flag for Household Selecting Multiple Plans

egen u_p_pre = mean(p_pre), by(year idhh)
g multplan = (u_p_pre!=p_pre) 
lab var multplan "HH splits plans"
drop u_p_pre

* Household Max Age, Age Curve

egen fage = max(age), by(year idhh)
lab var fage "HH Age Max"

egen fageadj = total(agecurve), by(year idhh)
lab var fageadj "HH Age Curve Total"

* Household Composition

egen fsize = count(idind), by(year idhh)
replace fsize = 5 if fsize > 5
lab var fsize "HH Size"

* Reduce to HH Level 

duplicates drop idhh year, force

* Save

save "proc\ca_enroll_hh_nofilter", replace

* Review Filters

tab year                     
bys year: count if multplan==1                            
bys year: count if gra==.                                    
bys year: count if fage<18           
bys year: count if fplcat==1                                         
bys year: count if hassub==1 & fpl>400                 
count // 1,096,430; Maintains 94.87% of original households

* Implement Filters

count                     
drop if multplan==1         // 63,137                        
drop if gra==.              // 17,226                            
drop if fage<18             // 31,685
drop if fplcat==1           // 80,293                                 
drop if hassub==1 & fpl>400 // 247,042                 
count // 1,096,430; Maintains 94.87% of original households

* Save

save "proc\ca_enroll_hh", replace

/// PLAN-LEVEL PREMIUMS 2014 ///

* Load

use "proc\ca_enroll_hh_nofilter" if year==2014 & gra!=., clear
rename (plan_name issuer_name) (planname issuername)

* Change to 14-Digit Plan ID

rename planid planid16
g planid = substr(planid,1,14)
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

merge m:1 gra issuername metal plantype using proc\prem14, keep(master match) nogen
replace prem = 233.61 if prem==. // Missing bronze hdhp prem set to equivalent bronze prem

* Save

order year gra planid issuername planname prem metal plantype
save proc/planchars14, replace

/// PLAN-LEVEL PREMIUMS 2015 ///

* Load

use "proc\ca_enroll_hh_nofilter" if year==2015 & gra!=., clear
rename (plan_name issuer_name) (planname issuername)

* Change to 14-Digit Plan ID

rename planid planid16
g planid = substr(planid,1,14)
duplicates drop planid gra, force
keep year gra planname issuername planid

* Issuer

tab issuername
replace issuername = "Anthem"       if issuername=="Anthem Blue Cross"
replace issuername = "Chinese"      if issuername=="Chinese Community"
replace issuername = "HealthNet"    if issuername=="Health Net"
replace issuername = "Molina"       if issuername=="Molina Health Car"
replace issuername = "Sharp"        if issuername=="SHARP Health Plan"
replace issuername = "Valley"       if issuername=="Valley Health"
replace issuername = "Western"      if issuername=="Western Health"

* Metal

g metal = ""
replace metal = "Catastrophic" if regexm(planname,"Minimu")==1 | regexm(planname,"Minum")==1
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
replace plantype = "HMO" if issuername=="Chinese"   & plantype==""
replace plantype = "HMO" if issuername=="HealthNet" & plantype=="" & gra>=15
replace plantype = "HMO" if issuername=="Kaiser"    & plantype==""
replace plantype = "HMO" if issuername=="Molina"    & plantype==""

tab planid if plantype==""
replace plantype = "PPO" if regexm(planid,"27603CA114")==1 & plantype==""

* Merge Premiums

merge m:1 gra issuername metal plantype using proc\prem15, keep(master match) nogen
drop if prem==. // 3 remaining ones seemed to be errors, plans not offered alongside other metals

* Save

order year gra planid issuername planname prem metal plantype
save proc/planchars15, replace

/// PREMIUMS 2016-2018 ///

* Plan Characteristics

import excel using data\prem2018.xlsx, first case(lower) cellrange(A8:K25916) clear
rename (ratingareaid individualrate applicant hoisid metallevel network product fullplanname) ///
(gra prem issuername issuerid metal plantype product planname)
keep if age=="21"
g year = 2018
destring gra, ignore("Rating Area ") replace
keep  year gra issuerid planid issuername prem metal plantype
order year gra issuerid planid issuername prem metal plantype
save proc/planchars18, replace // Not going to worry about this 'til I have 2018 enrollment

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
