/// SETUP ///

clear *
cd "\\files.umn.edu\cla\home\ryan0463\My Documents\CovCAInertia"

/// PREP ///

* Load
foreach n in 1 2 {


use "proc\analysis_i_191112", clear 
keep if samp5==1

if (`n'==1) {
keep if act_pass!="A"
}
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


tostring gra, gen(gra_str)

gen metal_gra = metal + gra_str
gen issuername_gra = issuername + gra_str
gen netname_gra = netname + gra_str
gen iss_net_gra = netname + gra_str

egen issfe_gra    = group(issuername_gra)
egen netfe_gra    = group(netname_gra)

* Default Product

g temp_def_product = product if iplan==1
egen def_product = max(temp_def_product), by (hh_year_id)

* Active Choice

gen active = 1
replace active = 0 if act_pass == "A"
lab var active "Active Chooser"


* Incorporate CSR Silver Variants
egen hh_csr = max(s_csr),by(hh_year_id)
replace hh_csr = 94 if fpl>=100 & fpl<=150
replace hh_csr = 87 if fpl>150 & fpl<=200
replace hh_csr = 73 if fpl>200 & fpl<=250

replace metal = "Silver73" if metal=="Silver" & hh_csr==73
replace metal = "Silver87" if metal=="Silver" & hh_csr==87
replace metal = "Silver94" if metal=="Silver" & hh_csr==94

* Metal FE

g mtl_brz = metal=="Bronze"
g mtl_cat = metal=="Catastrophic"
g mtl_gld = metal=="Gold"
g mtl_hdp = metal=="HDHP Bronze"
g mtl_plt = metal=="Platinum"
g mtl_slv = metal=="Silver" & hh_csr==0
g mtl_s73 = metal=="Silver73" 
g mtl_s87 = metal=="Silver87" 
g mtl_s94 = metal=="Silver94"

* Default Metal FE

g def_mtl_brz = lagmetal=="Bronze"
g def_mtl_cat = lagmetal=="Catastrophic"
g def_mtl_gld = lagmetal=="Gold"
g def_mtl_hdp = lagmetal=="HDHP Bronze"
g def_mtl_plt = lagmetal=="Platinum"
g def_mtl_slv = lagmetal=="Silver" & hh_csr==0
g def_mtl_s73 = lagmetal=="Silver" & hh_csr==73
g def_mtl_s87 = lagmetal=="Silver" & hh_csr==87
g def_mtl_s94 = lagmetal=="Silver" & hh_csr==94


* Rating Area FE

forvalues i = 1/19 { // 13 for LA
	g grafe_`i' = gra==`i'
}

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

replace padj = padj/100


*Drop Products that are never chosen

egen choicetotal = sum(choice), by(product year)
drop if choicetotal==0
drop choicetotal

*Drop Households with missing prices

egen missingprem = sum(padj==.), by(hh_id)
drop if missingprem>0
drop missingprem

* Drop Households with an active selection inconsistency 
gen inactive_switch_flag = iplan==1 &choice==0 & active==0
egen inactive_switch = sum(inactive_switch_flag), by(hh_id)
tab inactive_switch
drop if inactive_switch>0
drop inactive_switch_flag inactive_switch

* Inertial Plan Choice
egen has_inertia_plan = sum(iplan), by(idhh year)
*gen hasi2 = icatplan!=1
gen hasi = has_inertia_plan>0
lab var hasi "Has Inertia"
replace active = 1 if hasi==0



* Subsidy Flag
gen hassub = s_aptc>0

*Kaiser Flag
gen lag_kaiser = lagiss=="Kaiser"

*Kaiser Flag
gen inet_Kaiser = inet*lag_kaiser


* Fixed Effects for change in auto-premium
gen autodp_disc_1 = autodp<=-50
gen autodp_disc_2 = autodp >-50 & autodp<=-10
gen autodp_disc_3 = autodp >-10 & autodp<=0
gen autodp_disc_4 = autodp >0   & autodp<=10
gen autodp_disc_5 = autodp >10  & autodp<=50
gen autodp_disc_6 = autodp >50


if (`n'==1){
save "proc\analysis_i2_active", replace
}
if (`n'==2){
save "proc\analysis_i2", replace
}


keep year gra hh_id hh_year_id planid product padj choice metal issuername netname autodp autoplanchange /// YOU WILL WANT TO CHANGE THESE
hassub fam age3 iplan inet inet_Kaiser iiss iopt icatplan agefe_* def_* year_* grafe_* issfe_* ///
active autoelig samp1 samp5 metal_gra iss_net_gra autodp_disc*

*netfe_*, issfe_*

sort hh_id year product

if (`n'==1){
outsheet using "Output\analysis_i2_active.csv", comma replace nolabel
}
if (`n'==2){
outsheet using "Output\analysis_i2.csv", comma replace nolabel
}
}

/* CD Notes for CR
	- There is no indicator anymore for discontinued plans, but there were only ~800 total that didn't get reassigned over 5 yrs, none in LA
	- Default price given by auto_p_post
	- No more observations with missing padj */

/// OUTLINE FOR ATTENTION AND CHOICE MODELS ///

* Set Globals for Interactions (Still need to do this)

* Attention Stage
/*
logit latent_attention autodp agefe_1 agefe_2 fam s_hassub lagmetal6 grafe_* yearfe_*
	/* lagmetal6 captures metal level of default plan */
	/* demographic characteristics are: agefe_1 agefe_2 fam s_hassub */ 

* Choice Stage

clogit   choice padj padj#dem iplan iplan#dem inet inet#dem mtlfe_* netfe_*, group(hhidyr)
mixlogit choice padj padj#dem iplan iplan#dem inet inet#dem mtlfe_* netfe_*, group(hhidyr) rand(iss_kaiser iss_blue iss_anthem iss_hnet)
