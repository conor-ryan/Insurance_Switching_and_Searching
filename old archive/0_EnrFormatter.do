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

replace age = 100 if age >100 // 329

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

* Plan Timing

replace  cov_end = "" if regexm(cov_end,"W")==0
replace  cov_end = substr(cov_end, 6,2)
destring cov_end, replace

replace  cov_strt = substr(cov_strt,6,2)
destring cov_strt, replace

gen enr_strt2 = substr(enr_strt,6,2)
destring enr_strt2, replace
replace enr_strt2 = enr_strt2 - 53 if year==2017 & regexm(enr_strt,"2016")==1
replace enr_strt2 = enr_strt2 - 53 if year==2016 & regexm(enr_strt,"2015")==1
replace enr_strt2 = enr_strt2 - 53 if year==2015 & regexm(enr_strt,"2014")==1
drop enr_strt
rename enr_strt2 enr_strt
lab var enr_strt "Enrollee Enrollment Week"

* Drop

drop region plan_mrkt 
compress

* Save

save "D:\proc\ca_enroll_ind", replace

/// GRAB PLAN IDS FOR MATCH WITH PREMIUM DATA ///

* 2014

use "D:\proc\ca_enroll_ind" if year==2014, clear
replace planid = substr(planid,1,14)
duplicates drop planid, force
keep gra planid issuer_name plan_type plan_name

g metal = ""
replace metal = "Minimum Coverage" if regexm(plan_name,"Minimu")==1
replace metal = "Bronze"           if regexm(plan_name,"Bronze")==1
replace metal = "Silver"           if regexm(plan_name,"Silver")==1
replace metal = "Gold"             if regexm(plan_name,"Gold")==1
replace metal = "Platinum"         if regexm(plan_name,"Platin")==1
replace metal = "AI-AN Zero Cost"  if regexm(plan_name,"AI-AN")==1
drop if metal=="AI-AN Zero Cost"

g plantype = ""
replace plantype = "HMO" if regexm(plan_name,"HMO")==1
replace plantype = "PPO" if regexm(plan_name,"PPO")==1
replace plantype = "EPO" if regexm(plan_name,"EPO")==1

replace plantype = "EPO" if regexm(plan_name,"Blue Shield - Minimum Coverage EP")==1
replace plantype = "PPO" if regexm(plan_name,"Blue Shield - Minimum Coverage PP")==1
replace plantype = "PPO" if regexm(planid,"27603CA114")==1 // From CC Memb Profil 2014 09
replace plantype = "EPO" if regexm(planid,"27603CA116")==1 // From CC Memb Profil 2014 09

replace plantype = "HMO" if regexm(plan_name,"Chinese")==1
replace plantype = "HMO" if regexm(plan_name,"Contra Costa")==1
replace plantype = "HMO" if regexm(plan_name,"Kaiser")==1
replace plantype = "HMO" if regexm(plan_name,"LA Care")==1
replace plantype = "HMO" if regexm(plan_name,"Molina")==1
replace plantype = "HMO" if regexm(plan_name,"SHARP")==1
replace plantype = "HMO" if regexm(plan_name,"Sharp")==1
replace plantype = "HMO" if regexm(plan_name,"Western")==1
replace plantype = "HMO" if regexm(plan_name,"Valley")==1

rename issuer_name issuername
replace issuername = "Anthem"       if issuername=="Anthem Blue Cross"
replace issuername = "Chinese"      if issuername=="Chinese Community"
replace issuername = "Contra Costa" if issuername=="Contra Costa Heal"
replace issuername = "HealthNet"    if issuername=="Health Net"
replace issuername = "Molina"       if issuername=="Molina Health Car"
replace issuername = "Sharp"        if issuername=="SHARP Health Plan"
replace issuername = "Valley"       if issuername=="Valley Health"
replace issuername = "Western"      if issuername=="Western Health"

drop plan_type plan_name

save "proc/caenr14_planid_xwalk", replace

* 2015

use "D:\proc\ca_enroll_ind" if year==2015, clear
replace planid = substr(planid,1,14)
duplicates drop planid, force
keep gra planid issuer_name plan_type plan_name

g metal = ""
replace metal = "Minimum Coverage" if regexm(plan_name,"Minimu")==1
replace metal = "Bronze"           if regexm(plan_name,"Bronze")==1
replace metal = "Silver"           if regexm(plan_name,"Silver")==1
replace metal = "Gold"             if regexm(plan_name,"Gold")==1
replace metal = "Platinum"         if regexm(plan_name,"Platin")==1
replace metal = "AI-AN Zero Cost"  if regexm(plan_name,"AI-AN")==1
replace metal = "Minimum Coverage" if metal=="" // Not sure why this 1 isn't catching
drop if metal=="AI-AN Zero Cost"

g plantype = ""
replace plantype = "HMO" if regexm(plan_name,"HMO")==1
replace plantype = "PPO" if regexm(plan_name,"PPO")==1
replace plantype = "EPO" if regexm(plan_name,"EPO")==1
replace plantype = "HMO" if plantype=="" // Checked prod IDs against CC Memb Profile 2015 09

rename issuer_name issuername
replace issuername = "Anthem"            if issuername=="Anthem Blue Cross"
replace issuername = "HealthNet"         if issuername=="Health Net"
replace issuername = "Chinese"           if issuername=="Chinese Community"
replace issuername = "HealthNet"         if issuername=="Health Net"
replace issuername = "Molina"            if issuername=="Molina Health Car"
replace issuername = "Sharp"             if issuername=="SHARP Health Plan"
replace issuername = "Valley"            if issuername=="Valley Health"
replace issuername = "Western"           if issuername=="Western Health"

drop plan_type plan_name

save "proc/caenr15_planid_xwalk", replace

bys gra metal issuername plantype: g dup = _N // Just Sharp and HDHP Bronze

/// HOUSEHOLD CONVERSION ///

* Load Data

use "D:\proc\ca_enroll_ind", clear

* Flag for Household Selecting Multiple Plans

egen u_p_pre = mean(p_pre), by(year idhh)
g multplan = (u_p_pre!=p_pre) 
lab var multplan "HH splits plans"

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

**************************************************************************************************
 
/// HOUSEHOLD CONVERSION OLD ///

* Create 2016 File for Inertia

use idhh planid year using "D:\proc\ca_enroll_ind" if year==2016, clear

duplicates drop idhh planid, force
egen hhobs = count(idhh), by(idhh)
keep if hhobs == 1

rename planid planid_2016
keep idhh planid_2016

save "D:\proc\ca_enroll_16_inertia", replace

* Load 2017 Data

use "D:\proc\ca_enroll_ind" if year==2017, clear

* Flag for Household Selecting Multiple Plans

egen u_p_pre = mean(p_pre), by(year idhh)
g multplan = (u_p_pre!=p_pre) 
lab var multplan "HH splits plans"

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

* Merge Plan Selection

rename planid planid16
g planid = substr(planid,1,14)
merge m:1 year planid gra using "$cca\ca_rates", keepusing(planid metal) keep(master match) nogen // 3 master 

* Merge Plan Existence in 2016

replace year = 2016
merge m:1 year planid gra using "$cca\ca_rates", keepusing(planid) keep(master match)
g incplan = _merge==3
lab var incplan "Incumbent Plan (2017)"
drop _merge
replace year = 2017

* Modify Plan IDs to Standard Format

replace planid = planid + "-04" if substr(planid16,-2,2)=="04"
replace planid = planid + "-05" if substr(planid16,-2,2)=="05"
replace planid = planid + "-06" if substr(planid16,-2,2)=="06"
drop planid16

* Metal

replace metal = "Silver CSR 73" if substr(planid,-3,3)=="-04"
replace metal = "Silver CSR 87" if substr(planid,-3,3)=="-05"
replace metal = "Silver CSR 94" if substr(planid,-3,3)=="-06"

* Plan Design

merge m:1 year gra planid using "D:\proc\ca_plandesign", nogen keep(master match) // 3 master

* Inertia 

merge 1:1 idhh using "D:\proc\ca_enroll_16_inertia", keep(master match) nogen

g planid14      = substr(planid,1,14)
g planid14_2016 = substr(planid_2016,1,14)
g planidl2_2016 = substr(planid_2016,-2,2)

g hasinert = planid_2016!="" & incplan==1
lab var hasinert "Inertia"

g icat       = hasinert==0
replace icat = 2 if planid_2016!="" & incplan==0
replace icat = 3 if hasinert==1 & planid14!=planid14_2016
replace icat = 4 if hasinert==1 & planid14==planid14_2016
lab var icat "Inertia Category"
lab def icat 1 "New" 2 "Forced Switch" 3 "Switcher" 4 "Stayer"
lab val icat icat

g choice16 = planid14_2016
replace choice16 = choice16 + "-" + planidl2_2016 if planidl2_2016=="04" | planidl2_2016=="05" | planidl2_2016=="06"

drop planid14 planid14_2016 planid_2016 planidl2_2016 

* Save

save "D:\proc\ca_enroll_hh_nofilter", replace

/// APPLY FILTERS ///

* Identify Filters

use "D:\proc\ca_enroll_hh_nofilter", clear

table issuername plantype, c(freq)
table issuername plantype if metal=="Bronze", c(freq)
table issuername plantype if metal!="Bronze", c(freq) // Only epo/hmo/ppo-hsp

g socal = gra>=13
table issuername socal // Just split FE for Anthem, BCBS, and Kaiser

* Filter Counts (numbers reflect 2017, not mutually exclusive)

count                                        // 1,155,673

count if metal==""                           // Plan not in premium data; 3
count if multplan==1                         // HH splits plans; 24,377              
count if gra==.                              // Missing GRA; 1
count if fage<18                             // Not a non-elderly adult; 10,548
count if fplcat==1                           // FPL Medicaid Eligible; 21,853
count if metal=="Catastrophic" & fage>30     // Has Min Cvg. over 30; 399
count if hassub==1 & fpl>400                 // Has sub over 400; 3,213

* Filters

drop if metal==""                      
drop if multplan==1                        
drop if gra==.                            
drop if fage<18             
drop if fplcat==1                                             
drop if metal=="Catastrophic" & fage>30  
drop if hassub==1 & fpl>400                  

count // 1,096,430; Maintains 94.87% of original households

drop u_p_pre multplan 
save "D:\proc\ca_enroll_hh", replace

//////////////
/// TABLES ///
//////////////

* Table 2: Population

use "D:\proc\ca_enroll_hh_nofilter" if year==2017, clear
	table gra, c(freq)
	table issuername, c(freq)
	table issuername gra, c(freq)

* Table 2: Choices by GRA

use "D:\proc\ca_plandesign" if year==2017, clear
	drop if regexm(metal,"Enhanced")==1 | metal=="Minimum Coverage"
	table gra, c(freq)
	duplicates drop gra issuername, force
	table gra, c(freq)
	
* Table 2: Market Shares, HHI

use "D:\proc\ca_enroll_hh_nofilter" if year==2017, clear // maybe switch to limited sample w/ ca_enroll_hh?
	qui collapse (count) age, by(issuername gra)
	qui rename age issenrgra
	qui egen enrgra = total(issenrgra), by(gra)
	qui gen  issenrpctgra = (issenrgra / enrgra) * 100
	table issuername gra, c(m issenrpctgra) format(%9.2f)
	table gra issuername, c(m issenrpctgra) format(%9.2f)
		g mktshr2 = issenrpctgra^2
		collapse (sum) mktshr2, by(gra)
		rename mktshr2 hhi
		drop if gra==.
		table gra, c(m hhi)
	
* Table 2: Lowest Silver by GRA

use "D:\proc\ca_plandesign" if year==2017 & metal=="Silver", clear
	table gra, c(min prem)
	
* Table 3: Market Shares by Insurer and GRA

use "D:\proc\ca_enroll_hh_nofilter" if year==2017, clear // maybe switch to limited sample w/ ca_enroll_hh?
	qui collapse (count) age, by(issuername gra)
	qui rename age issenrgra
	qui egen enrgra = total(issenrgra), by(gra)
	qui gen  issenrpctgra = (issenrgra / enrgra) * 100
	table issuername gra, c(m issenrpctgra) format(%9.2f)
	table gra issuername, c(m issenrpctgra) format(%9.2f)
	
* Inertia

use planid year using "D:\proc\ca_plandesign" if year==2016, clear
	duplicates drop planid, force
	save "D:\proc\ca_planid16", replace

use idhh reenr planid choice16 using "D:\proc\ca_enroll_hh", clear
	merge m:1 planid using "D:\proc\ca_planid16", keep(master match)
	sum reenr // pct from 2016
	
	gen inert = reenr==1 & _merge==3
	sum inert // pct enr with inertia
	sum inert if reenr==1 // pct from 2016 with inertia
	
	gen loyal = reenr==1 & _merge==3 & planid==choice16
	sum loyal // pct loyal enrs
	sum loyal if reenr==1 // pct loyal enrs among reenrollers
	sum loyal if inert==1 // pct loyal enrs among those with inertia
	
	sum reenr inert loyal
	sum inert loyal if reenr==1
	sum loyal if inert==1
	
* Market Share by Enrollment Year for Top 5: 1/23 -- DROP PPL WHOSE PLANS WERE DISCONTINUED!!!

use "D:\proc\ca_enroll_hh_nofilter", clear
	sort idhh year
	gen choice16 = planid[_n-1]
	gen reenr    = idhh==idhh[_n-1]
	gen inert    = planid==choice16 & year==2017 
	order idhh year planid choice16 reenr inert
	
	drop if missplan==1                       
	drop if multplan==1                        
	drop if gra==.                            
	drop if fage<18 | fage>64                
	drop if fpl<133 & fpl>0                     
	drop if fpl==.                               
	drop if metal=="Minimum Coverage" & fage>30  
	drop if hassub==1 & fpl>400 
	
	bys idhh: g hhyrobs = _N
	drop if hhyrobs==1 & year==2016
	
	gen enrgrp = year==2016 
	replace enrgrp = 2 if year==2017 & hhyrobs==2
	replace enrgrp = 3 if year==2017 & hhyrobs==1

	table issuername enrgrp, c(freq)
	