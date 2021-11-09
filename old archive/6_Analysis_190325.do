/// SETUP ///

clear *
cd "\\files.umn.edu\cla\Home\ryan0463\My Documents\CovCAInertia"

/// PREP ///

* Load

use "proc\analysis_i" if year<2018, clear // Need full 2018 data
*keep if gra==16 // LA Only for now. 

* Groups 
replace netname = issuername + " " + plantype
egen issfe = group(issuername)
egen netfe = group(netname)
egen mtlfe = group(metal)
egen household = group(idhh)
egen product = group(planid)
egen hh_year_id  = group(household year)
lab var netname "Network Name"


* Default Product
g temp_def_product = product if iplan==1
egen def_product = max(temp_def_product), by (hh_year_id)

* Active Choice
gen active = 1
replace active = 0 if act_pass == "A"

*Quadratic and Threshold in plan options
gen ioptsq = iopt^2

* Metal FE

g mtl_brz = metal=="Bronze"
g mtl_cat = metal=="Catastrophic"
g mtl_gld = metal=="Gold"
g mtl_hdp = metal=="HDHP Bronze"
g mtl_plt = metal=="Platinum"
g mtl_s73 = metal=="Silver: CSR 73"
g mtl_s87 = metal=="Silver: CSR 87"
g mtl_s94 = metal=="Silver: CSR 94"

* Default Metal FE
g def_metal = metal if iplan==1
sort hh_year_id def_metal
bysort hh_year_id: replace def_metal = def_metal[_N]


g def_mtl_brz = def_metal=="Bronze"
g def_mtl_cat = def_metal=="Catastrophic"
g def_mtl_gld = def_metal=="Gold"
g def_mtl_hdp = def_metal=="HDHP Bronze"
g def_mtl_plt = def_metal=="Platinum"
g def_mtl_s73 = def_metal=="Silver: CSR 73"
g def_mtl_s87 = def_metal=="Silver: CSR 87"
g def_mtl_s94 = def_metal=="Silver: CSR 94"

// foreach var in mtl_brz mtl_cat mtl_gld mtl_hdp mtl_plt mtl_s73 mtl_s87 mtl_s94 {
// g temp_def_`var' = `var' if iplan==1
// egen def_`var' = max(temp_def_`var'), by(hh_year_id)
// }
// drop temp_def*

* Issuer FE

forvalues i = 2/13 {
	g issfe_`i' = issfe==`i'
}

* Default Issuer FE
g temp_def_issfe = issfe if iplan==1
egen def_issfe=max(temp_def_issfe), by(hh_year_id)

* Net FE

forvalues i = 1/23 {
	if `i'==14 continue
	g netfe_`i' = netfe==`i'
}
* Default Issuer FE
g temp_def_netfe = netfe if iplan==1
egen def_netfe = max(temp_def_netfe), by(hh_year_id)


* Age FE

forvalues i = 0/2 {
	g agefe_`i' = age3==`i'
}

* Scaling 

replace padj = sqrt(padj / 100)
replace iopt = iopt / 10

* Default Price
g temp_def_padj = padj if iplan==1
egen def_padj = max(temp_def_padj), by(hh_year_id)

drop temp_def*

*Discontinued
gen lag_disc = icatplan==4



/* Numeric Product and HH ID
gen product = regexr(planid,"CA","")
destring product, replace

encode idhh, gen(household)
*/
* Cleaning 

keep year gra idhh hh_year_id household planid product padj choice metal mtlfe issuername issfe issfe_* netname netfe ///
pmtlord hassub fam age3 iplan inet iiss itenplan iexp iopt icatplan dcat pshock mtl_* netfe_* agefe_* def_* dprem active lag_disc
rename (itenplan icatplan) (iten icat)
order year gra idhh planid household product hh_year_id choice iplan padj choice 
sort household year gra 

* Drop Missing Prices (Why are these here?)
drop if padj==.

*Drop Products that are never chosen
egen choicetotal = sum(choice),by(product year)
drop if choicetotal==0
drop choicetotal

* Has Inertia

egen hasi  = sum(iplan), by(hh_year_id)
lab var hasi "Has Inertia"

* Drop observations missing dprem
drop if lag_disc==0 & dprem==. & hasi==1

* Save

save "proc\analysis_i2", replace
outsheet using "Output\analysis_i2.csv", comma replace nolabel

/// CHOICE MODELS ///

* Load 

use "proc\analysis_i2", clear
global inert c.iplan#c.iten c.iplan#c.iexp c.iplan#c.iopt // add pshock?
global dem_int c.padj#c.fam c.padj#c.agefe_1 c.padj#c.agefe_2 c.padj#c.hassub ///
		c.iplan#c.fam c.iplan#c.agefe_1 c.iplan#c.agefe_2 c.iplan#c.hassub 
global conv  tolerance(0.001) ltolerance(0.001) nrtolerance(0.001) diff


* Conditional Logit
clogit choice c.padj c.iplan c.inet c.iiss i.mtlfe i.netfe if hasi==1 & dcat==1, group(hh_year_id)
outreg using "CL_results.doc", replace
clogit choice c.padj c.iplan $dem_int c.inet c.iiss i.mtlfe i.netfe if hasi==1 & dcat==1, group(hh_year_id)
outreg using "CL_results.doc", merge


/*
clogit choice c.padj c.padj#c.iplan c.padj#c.dprem#c.iplan c.padj#c.active#c.iplan c.iplan c.iplan#c.dprem c.iplan#c.active ///
 i.mtlfe i.netfe if hasi==1 & dcat==1, group(hh_year_id)
outreg using "CL_results.doc", merge
clogit choice c.padj c.padj#c.iplan c.padj#c.dprem#c.iplan c.padj#c.active#c.iplan c.iplan c.iplan#c.dprem c.iplan#c.active ///
$dem_int i.mtlfe i.netfe if hasi==1 & dcat==1, group(hh_year_id)
outreg using "CL_results.doc", merge

clogit choice c.padj c.iplan c.padj#c.iplan ///
i.mtlfe i.mtlfe#c.iplan i.netfe if hasi==1 & dcat==1, group(hh_year_id)
outreg using "CL_results.doc", merge
clogit choice c.padj c.iplan c.padj#c.iplan ///
$dem_int  i.mtlfe i.mtlfe#c.iplan i.netfe if hasi==1 & dcat==1, group(hh_year_id)
outreg using "CL_results.doc", merge
*/

* DSC Logit
alogit choice padj iplan inet iiss $dem_int mtl_* netfe_* if hasi==1 & dcat==1, /// might also test inertia for population that terminated
group(hh_year_id) def(iplan) avars(agefe* fam hassub active dprem def_padj def_mtl_*) model(dsc) diff
outreg using "DSC_results.doc", replace
alogit choice padj iplan inet iiss  mtl_* netfe_* if hasi==1 & dcat==1, /// might also test inertia for population that terminated
group(hh_year_id) def(iplan) avars(agefe* fam hassub active dprem def_padj def_mtl_*) model(dsc) diff
outreg using "DSC_results.doc", merge

alogit choice padj iplan inet iiss $dem_int mtl_* netfe_* if hasi==1 & dcat==1, /// might also test inertia for population that terminated
group(hh_year_id) def(iplan) avars(agefe* fam hassub dprem def_padj def_mtl_*) model(dsc) diff
outreg using "DSC_results.doc", merge
alogit choice padj iplan inet iiss  mtl_* netfe_* if hasi==1 & dcat==1, /// might also test inertia for population that terminated
group(hh_year_id) def(iplan) avars(agefe* fam hassub dprem def_padj def_mtl_*) model(dsc) diff
outreg using "DSC_results.doc", merge

// di "$S_TIME" // status, age
// alogit choice padj iplan mtl_* netfe_* if hasi==1 & dcat==1, /// might also test inertia for population that terminated
// group(hh_year_id) def(iplan) exclude(padj iplan mtl_* netfe_*) model(dsc) diff ///
// b0(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,00,0,0,0,0,0,0,0) d0(0,0)
//
// di "$S_TIME"
// clogit choice padj iplan $inert mtl_* netfe_* if dcat!=1 | (dcat==1 & hasi==0), group(hhid) diff
//
// di "$S_TIME"



// * Conditional Logit without Auto Enrollment
//
// clogit choice c.padj c.iplan $inert i.mtlfe i.netfe if gra==16, group(hhid)
// est sto pcnt
// est save "proc/pcnt", replace
// clogit choice i.pmtlord c.iplan $inert i.mtlfe i.netfe if gra==16, group(hhid)
// est sto pord
// est save "proc/pord", replace
// esttab pcnt pord, drop(*netfe) b(%9.3f) nostar nobase wide stats(N ll r2_p)

* Mixed Logit
eststo mlpcnt: asmixlogit choice $dem_int  i.mtlfe i.netfe if hasi==1 & dcat==1, case(hh_year_id) alt(planid) nocons random(padj iplan)
outreg using "ML_results.doc", replace

eststo mlpcnt: asmixlogit choice c.padj#c.iplan c.padj#c.dprem#c.iplan c.padj#c.active#c.iplan c.iplan#c.dprem c.iplan#c.active $dem_int  i.mtlfe i.netfe if hasi==1 & dcat==1, case(hh_year_id) alt(planid) nocons random(padj iplan)
outreg using "ML_results.doc", merge

// eststo mlpcnt: asmixlogit choice padj  $inert i.mtlfe i.netfe if gra==16, case(hhid) alt(planid) nocons random(iplan)
// est save "proc/mlpcnt", replace
// eststo mlpord: asmixlogit choice i.pmtlord  $inert i.mtlfe i.netfe if gra==16, case(hhid) alt(planid) nocons random(iplan)
// est save "proc/mlpord", replace
//
// esttab ml*, drop(*netfe) b(%9.3f) nostar nobase wide stats(N ll r2_p)
