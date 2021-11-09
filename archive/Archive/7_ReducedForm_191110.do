/// SETUP ///

clear *
cd "G:\Shared Drives\CovCAInertia"

/// PREP ///

* Load

use "proc\analysis_i", clear // Need full 2018 data

* Groups 
egen issfe      = group(issuername)
egen netfe      = group(netname)
egen mtlfe      = group(metal6)
egen household  = group(idhh)
egen product    = group(planid)
egen hh_year_id = group(household year)
egen hh_id      = group(household)

* Default Product
g temp_def_product = product if iplan==1
egen def_product = max(temp_def_product), by (hh_year_id)

* Active Choice
gen active = 1
replace active = 0 if act_pass == "A"

* Fake Active Choice
gen draw = uniform()
gen fake_active = 0 
gen switcher = icatplan==3
replace fake_active = 1 if (switcher==1 & active==1) // Control for discontinued plans?
replace fake_active =1 if (draw<.695 & fake_active==0)

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

forvalues i = 1/13 {
	g issfe_`i' = issfe==`i'
}

* Default Issuer FE
g temp_def_issfe = issfe if iplan==1
egen def_issfe=max(temp_def_issfe), by(hh_year_id)

forvalues i = 1/13 {
	g def_issfe_`i' = def_issfe==`i'
}


* Net FE

forvalues i = 1/23 {
	if `i'==14 continue
	g netfe_`i' = netfe==`i'
}
* Default Network FE
g temp_def_netfe = netfe if iplan==1
egen def_netfe = max(temp_def_netfe), by(hh_year_id)

forvalues i = 1/23 {
	if `i'==14 continue
	g def_netfe_`i' = def_netfe==`i'
}

* Age FE

forvalues i = 0/2 {
	g agefe_`i' = age3==`i'
}

* Year FE

forvalues i = 2014/2018 {
	g year_`i' = year==`i'
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

*** Kaiser Switching
gen iplan_kais = iplan*(issuername=="Kaiser")
gen inet_kais = inet*(issuername=="Kaiser")




/* Numeric Product and HH ID
gen product = regexr(planid,"CA","")
destring product, replace

encode idhh, gen(household)
*/
* Cleaning 


rename (itenplan icatplan) (iten icat)
order gra hh_id year planid household product hh_year_id choice iplan padj choice 
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



save "proc\analysis_i3", replace


/// REDUCED FORM MOTIVATION ///
use "proc\analysis_i3", clear // Need full 2018 data

// egen hh_year_id  = group(idhh year)
//
// * Active Choice
// gen active = 1
// replace active = 0 if act_pass == "A"
//
// egen hasi  = sum(iplan), by(hh_year_id)
// lab var hasi "Has Inertia"
//
// *Discontinued
// gen lag_disc = icatplan==4
// * Drop observations missing dprem
// drop if lag_disc==0 & dprem==. & hasi==1

* Keep Only Returning Consumers
keep if hasi==1


collapse (sum) choice (mean) padj (mean) dprem, by(hh_year_id iplan year gra active dpremmed dpremmean  age3 fam agefe_* hassub def_*)

gen switch = iplan==choice
replace dprem = dprem/100
replace dpremmed = dpremmed/100
replace dpremmean = dpremmean/100

gen padj_rival_temp = -100
replace padj_rival_temp = padj if choice==0
egen padj_rival = max(padj_rival_temp), by(hh_year_id)
drop padj_rival_temp

** Regression, price level
reg active padj padj_rival agefe_1 agefe_2 fam hassub i.year if choice==1
reg switch padj padj_rival agefe_1 agefe_2 fam hassub i.year if choice==1

** Regression, price change
reg active dprem dpremmean agefe_1 agefe_2 fam hassub i.year if choice==1
reg switch dprem dpremmean agefe_1 agefe_2 fam hassub i.year if choice==1


