*********************************************************************
* Convert HH data to HH-Choice data
*********************************************************************

/// SETUP ///

clear *
cd "\\files.umn.edu\cla\home\ryan0463\My Documents\CovCAInertia"
set seed 1234 

/// SUBSAMPLES - Need to adjust to sample by household/entering cohort ///

use "proc\ca_enroll_hh_wi", clear
sort gra idhh year

collapse (first) year, by(idhh gra)
rename year entry_year

save "proc\hh_sample_main", replace

*foreach n in 20 15 10 5 1 {
foreach n in 1 {

use "proc\hh_sample_main", clear
sample `n', by(gra)

merge 1:m idhh gra using "proc\ca_enroll_hh_wi"
keep if _merge==3
	
rename planid planchoice

save "proc\ca_enroll_hh_wi_`n'", replace

}

use "proc\ca_enroll_hh_wi", clear
	foreach n in 1 {
	merge 1:1 year idhh using "proc\ca_enroll_hh_wi_`n'", keepusing(idhh)
		g samp`n' = _merge==3
		lab var samp`n' "`n' Percent Sample"
		drop _merge
	}
*keep if samp1==1 | samp5==1
keep if samp1==1
rename planid planchoice
foreach x of varlist issuername issid netname plantype metal csr prem p_pre p_post aptc hassub pmtlord {
	rename `x' s_`x'
}
save "proc\ca_enroll_hh_wi2", replace

/// RESHAPE: HH --> HH-CHOICE ///

* Prep Plan Choice Data
	
use "proc/planchars", clear
	qui sum gra
	global n = r(N)
	di $n
	g id = _n
	save "proc\planchars_temp", replace
	
* Create Lists of Enrollees Whom May Enroll in Each Plan
	
forvalues x = 1/$n {

	di as error "`x'"
	qui {
	use "proc\planchars_temp" if id==`x', clear
	merge 1:m year gra using "proc\ca_enroll_hh_wi2", keep(match) nogen
	drop if metal=="Catastrophic" & fage>=30 // Remove catastrophic from choice set if fage>=30
	save "proc\planchoice_i_`x'", replace

	}
}
	
* Append Lists

forvalues x = 1/$n {
	di as error "`x'"
	if `x'==1 use "proc\planchoice_i_`x'", clear
	else qui append using "proc\planchoice_i_`x'"
	}
	
	save "proc\planchoice_i_app", replace
	
/// PREPARE FOR ANALYSIS ///

* Load 

use "proc\planchoice_i_app", clear
	
* Identify Chosen Plan

gen choice = (planid==planchoice)
lab var choice "Selected Plan"
tab choice 

egen choicecheck = total(choice), by(year idhh)
drop if choicecheck!=1 // <1%
drop choicecheck id
count

* Adjusted Premiums: APTC Defined as in CCA Data

gen     premadj_ccatemp = s_p_pre / s_prem if choice==1       // Determine age rating
egen    premadj_cca = max(premadj_ccatemp), by(idhh year) // Expand to other choices
replace premadj_cca = premadj_cca * prem                  // Age rate HH's options

gen     padj = premadj_cca - s_aptc                  // Deduct HH's actual APTC from all options
replace padj = 1           if padj < 1               // Correct APTC > Prem
replace padj = premadj_cca if metal=="Catastrophic"  // Correct for no subs on cat. plans

lab var padj "Post-APTC Premium"

drop premadj_* 

* Inertia Indicator 

order year gra idhh netname metal lagnet s_metal // should have lag metal

g iplan = netname == lagnet & metal == lagmetal & product==lagproduct
g inet  = netname == lagnet

lab var iplan "Plan Inertia"
lab var inet  "Network Inertia"

* Inertia: Plan Options

cap drop iopt
bys year gra idhh: g iopt = _N
lab var iopt "Inertia: Plan Options"

* Clean

drop netid
sort year gra idhh	

order year gra idhh samp1 samp5 planid planchoice choice ///
iplan inet icatplan icatnet act_pass inertxwalkmatch ///
issuerid issuername netname plantype metal padj prem pmtlord bp ///
fpl fplcat fage age3 fageadj fam fsize svch iopt status ///
enr_strt cov_strt cov_end s_* lag* auto*

lab var planid     "Plan ID: HH-Choice"
lab var planchoice "Chosen Plan ID: HH"
lab var choice     "Plan Chosen (0/1)"

* Save

save "proc\analysis_i_191112", replace
