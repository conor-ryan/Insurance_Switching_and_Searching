/// SETUP ///

clear *
cd "\\files.umn.edu\cla\home\ryan0463\My Documents\CovCAInertia"

/// PREP ///

* Load

use "proc\analysis_i_191112", clear 
keep if gra==16
rename product prodline

* Groups for Fixed Effects

egen issfe    = group(issuername)
egen netfe    = group(netname)
egen defissfe = group(autoiss)
egen defnetfe = group(autonet)
egen mtlfe    = group(metal)
egen hh_id     = group(idhh)
egen product  = group(planid)
egen hh_year_id   = group(idhh year)

* Default Product

g temp_def_product = product if iplan==1
egen def_product = max(temp_def_product), by (hh_year_id)

* Active Choice

gen active = 1
replace active = 0 if act_pass == "A"
lab var active "Active Chooser"

* Fake Active Choice

gen draw = uniform()
gen fake_active = 0 
gen switcher = icatplan==3 // Switchers
replace fake_active = 1 if (switcher==1 & active==1) // Control for discontinued plans?
replace fake_active = 1 if (draw<.695 & fake_active==0)

* Metal FE

g mtl_brz = metal=="Bronze"
g mtl_cat = metal=="Catastrophic"
g mtl_gld = metal=="Gold"
g mtl_hdp = metal=="HDHP Bronze"
g mtl_plt = metal=="Platinum"
g mtl_slv = metal=="Silver"

* Default Metal FE

g def_mtl_brz = lagmetal=="Bronze"
g def_mtl_cat = lagmetal=="Catastrophic"
g def_mtl_gld = lagmetal=="Gold"
g def_mtl_hdp = lagmetal=="HDHP Bronze"
g def_mtl_plt = lagmetal=="Platinum"
g def_mtl_slv = lagmetal=="Silver"

* Issuer FE

levelsof issfe, local(issfe)
foreach i in `issfe' { // 7 for LA
	g issfe_`i' = issfe==`i'
}

* Default Issuer FE

// levelsof autoiss, local(autoiss)
// foreach i in `autoiss' {
// 	g def_issfe = autoiss==`i'
// }

* Net FE

levelsof netfe, local(netfe)
forvalues i = 1/13 { // 13 for LA
	g netfe_`i' = netfe==`i'
}

* Default Network FE

// levelsof autoiss, local(autoiss)
// foreach i in `autoiss' {
// 	g def_issfe = autoiss==`i'
// }

* Age FE

forvalues i = 0/2 {
	g agefe_`i' = age3==`i'
}

* Year FE

forvalues i = 2014/2018 {
	g year_`i' = year==`i'
}

* Scaling 

replace padj = sqrt(padj/100)

*Drop Products that are never chosen

egen choicetotal = sum(choice), by(product year)
drop if choicetotal==0
drop choicetotal

* Inertial Plan Choice

gen hasi = icatplan!=1
lab var hasi "Has Inertia"

* Inertial Network Choice

* Drop observations missing change in premium

drop if autodp==. & hasi==1 // 30 total

* Drop observations where the auto assigned plan and iplan are different
* ~87 households in 1% sample in LA 
gen autoplan = autoplanid==planid
egen hasautoplan = sum(autoplan), by(hh_year_id)
gen diff = hasautoplan==1 & iplan==1 & autoplan==0
egen diffsum = sum(diff), by(idhh)
drop if diffsum==1
drop diffsum diff autoplan
* Still have observations where people have autoassigned plans, but those plans don't appear in the choice set. 


gen autochange = hasautoplan==1 & autoplan==1 & iplan==0


save "proc\analysis_i2", replace

keep year gra hh_id hh_year_id planid product padj choice metal mtlfe issuername issfe issfe_* netname netfe dpremmed /// YOU WILL WANT TO CHANGE THESE
pmtlord hassub fam age3 iplan inet iiss itenplan iexp iopt icatplan dcat pshock mtl_* netfe_* agefe_* def_* year_* ///
dprem active lag_disc fake_active iplan_* inet_* samp1 // samp5

outsheet using "Output\analysis_i2.csv", comma replace nolabel

/* CD Notes for CR
	- There is no indicator anymore for discontinued plans, but there were only ~800 total that didn't get reassigned over 5 yrs, none in LA
	- Default price given by auto_p_post
	- No more observations with missing padj */

/// OUTLINE FOR ATTENTION AND CHOICE MODELS ///

* Set Globals for Interactions (Still need to do this)

* Attention Stage

logit latent_attention autodp agefe_1 agefe_2 fam s_hassub lagmetal6 grafe_* yearfe_*
	/* lagmetal6 captures metal level of default plan */
	/* demographic characteristics are: agefe_1 agefe_2 fam s_hassub */ 

* Choice Stage

clogit   choice padj padj#dem iplan iplan#dem inet inet#dem mtlfe_* netfe_*, group(hhidyr)
mixlogit choice padj padj#dem iplan iplan#dem inet inet#dem mtlfe_* netfe_*, group(hhidyr) rand(iss_kaiser iss_blue iss_anthem iss_hnet)
