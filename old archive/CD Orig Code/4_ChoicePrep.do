*********************************************************************
* Convert HH data to HH-Choice data
*********************************************************************

/// SETUP ///

clear *
cd "M:\My Documents\CovCAInertia"

/// SUBSAMPLES ///

use "proc\ca_enroll_hh_wi", clear

rename planid planchoice
	
set seed 1234 
sample 20, by(year gra)
save "proc\ca_enroll_hh_wi_20", replace
sample 75, by(year gra)
save "proc\ca_enroll_hh_wi_15", replace
sample 66.66, by(year gra)
save "proc\ca_enroll_hh_wi_10", replace
sample 50, by(year gra)
save "proc\ca_enroll_hh_wi_5", replace
	
use "proc\ca_enroll_hh_wi", clear
	forvalues n = 5(5)20 {
	merge 1:1 year idhh using "proc\ca_enroll_hh_wi_`n'", keepusing(idhh)
		g samp`n' = _merge==3
		lab var samp`n' "`n' Percent Sample"
		drop _merge
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
	merge 1:m year gra using "proc\ca_enroll_hh_wi_5", keep(match) nogen
	drop if metal=="Catastrophic" & fage>30 // Remove catastrophic from choice set if >30
	replace metal = "Silver: CSR 73" if metal=="Silver" & csr==73 // Change silver metal depending on CSR availability
	replace metal = "Silver: CSR 87" if metal=="Silver" & csr==87
	replace metal = "Silver: CSR 94" if metal=="Silver" & csr==94
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

gen     premadj_ccatemp = p_pre / prem if choice==1       // Determine age rating
egen    premadj_cca = max(premadj_ccatemp), by(idhh year) // Expand to other choices
replace premadj_cca = premadj_cca * prem                  // Age rate HH's options

gen     padj = premadj_cca - aptc                    // Deduct HH's actual APTC from all options
replace padj = 1           if padj < 1               // Correct APTC > Prem
replace padj = premadj_cca if metal=="Catastrophic"  // Correct for no subs on cat. plans

lab var padj "Post-APTC Premium"

drop premadj_* 

* Age Groups

g age3 = 0 
replace age3 = 1 if fage>=30 & fage<50
replace age3 = 2 if fage>=50
lab var age3 "HH Maximum Age"
lab def age3 0 "18-30" 1 "31-49" 2 "50-64"
lab val age3 age3

* Family Indicator

g fam = fsize>1
lab var fam "Family"

* Inertia Indicator

g iplan = planid==lagplan
g inet  = planid==lagnet
g iiss  = planid==lagiss

lab var iplan "Plan Inertia"
lab var inet  "Network Inertia"
lab var iiss  "Issuer Inertia"

* Inertia: Plan Options

cap drop iopt
bys year gra idhh: g iopt = _N
lab var iopt "Inertia: Plan Options

* Save

sort  year gra idhh	
order year gra idhh planid planchoice choice padj metal
	
save "proc\analysis_i", replace
	