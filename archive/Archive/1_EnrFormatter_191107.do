********************************************
* CovCA Formatter - Individual Data
********************************************

/// INDIVIDUAL ENROLLMENT ///

* Load

clear *
cap cd "\\files.umn.edu\cla\home\ryan0463\My Documents\CovCAInertia"
cap cd "G:\Shared drives\CovCAInertia"
import delimited using "Data\CC PRA enrollment file_082018_v2", clear

* Rename and Label

rename (ahbx_case_id_x indv_id_x enrlee_enrlmnt_yr enrlee_status rating_area) (idhh idind year status gra)
rename (hios_id_16 active_passive_d subsidy_fpl_percent_round gross_premium_amt_round net_premium_amt_round aptc_amt_round) (planid act_pass fpl p_pre p_post aptc)
rename (zip_3 enrlee_plan_select_dt_wk cov_start_dt_wk cov_end_dt_wk) (zip3 enr_strt cov_strt cov_end)
rename (service_channel metal_level_enhanced plan_network_type) (svch metal plantype)

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
lab var region      "Rating Area"
lab var plan_name   "Plan Name"
lab var plantype    "Plan Type"
lab var planid      "Plan ID: 16"
lab var status      "Enrollment Status"
lab var year        "Year"
lab var zip3        "3-Digit Zip Code"
lab var svch        "Service Channel"
lab var metal       "Metal Level"

* Age Curve

replace age = 100 if age > 100 // 956
merge m:1 age using "data\agecurve", keep(master match) nogen // 11 missing age

* APTC

count if aptc==. // 494k
replace aptc = 0 if p_pre - p_post==0 & aptc==. // 494k

* Subsidy
	
g hassub = aptc>0
lab var hassub "Enrollee Has subsidy"
lab def hassub 0 "Unsubsidized" 1 "Subsidized"
lab val hassub hassub

* FPL 

tab year if fpl == . // Huge chunk missing in 2014
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

drop region
compress
save "proc\ca_enroll_ind", replace

/// HOUSEHOLD ENROLLMENT ///

* Load

use "proc\ca_enroll_ind", clear

* Flag for Household Selecting Multiple Plans

egen u_p_pre = mean(p_pre), by(year idhh)
g multplan = (u_p_pre!=p_pre) // 2.91% of HH split
lab var multplan "HH splits plans"
drop u_p_pre

* Household Max Age, Age Curve

egen fage = max(age), by(year idhh)
lab var fage "HH Age Max"

egen fageadj = total(agecurve), by(year idhh)
lab var fageadj "HH Age Curve Total"

* Age Groups

g age3 = 0 
replace age3 = 1 if fage>=30 & fage<50
replace age3 = 2 if fage>=50
lab var age3 "HH Maximum Age"
lab def age3 0 "18-30" 1 "31-49" 2 "50-64"
lab val age3 age3

* Household Composition

egen fsize = count(idind), by(year idhh)
replace fsize = 5 if fsize > 5
lab var fsize "HH Size"

g fam = fsize>1
lab var fam "Family"

* Reduce to HH Level 

duplicates drop idhh year, force

* Save

save "proc\ca_enroll_hh_nofilter", replace

/// FILTERS ///

* Load

use "proc\ca_enroll_hh_nofilter", clear

* Create CSR Filter

g planid16_2 = substr(planid,-2,2)
g csr = 0 
replace csr = 73 if planid16_2=="04"
replace csr = 87 if planid16_2=="05"
replace csr = 94 if planid16_2=="06"
drop planid16_2

* Review Filters

tab year                     
count if multplan==1                            
count if gra==.                                    
count if fage<18     
count if csr>0 & hassub==0                       

* Implement Filters

count // 5,464,510              
  
drop if multplan==1       // 88,309                        
drop if gra==.            // 17,196                            
drop if fage<18           // 41,686        
drop if csr>0 & hassub==0 // 11,575

count // 5,305,374; Maintains 97.09% of original households -- later remove 6149 nonmatches to plan IDs... 5299225, 96.97%

* Plan ID 14

rename planid planid16
g planid = substr(planid16,1,14)

* Save

drop idind gender age zip3 multplan
order year idhh planid planid16
save "proc\ca_enroll_hh", replace
