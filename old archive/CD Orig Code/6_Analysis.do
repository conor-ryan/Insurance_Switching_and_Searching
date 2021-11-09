/// SETUP ///

clear *
cd "M:\My Documents\CovCAInertia"

/// PREP ///

* Load

use "proc\analysis_i" if year<2018, clear // Need full 2018 data

* Groups 
replace netname = issuername + " " + plantype
egen issfe = group(issuername)
egen netfe = group(netname)
egen mtlfe = group(metal)
egen hhid  = group(year idhh)
lab var netname "Network Name"

* Has Inertia

egen hasi  = total(iplan), by(year idhh)
lab var hasi "Has Inertia"

* Metal FE

g mtl_brz = metal=="Bronze"
g mtl_cat = metal=="Catastrophic"
g mtl_gld = metal=="Gold"
g mtl_hdp = metal=="HDHP Bronze"
g mtl_plt = metal=="Platinum"
g mtl_s73 = metal=="Silver: CSR 73"
g mtl_s87 = metal=="Silver: CSR 87"
g mtl_s94 = metal=="Silver: CSR 94"

* Issuer FE

forvalues i = 2/13 {
	g issfe_`i' = issfe==`i'
}

* Net FE

forvalues i = 1/23 {
	if `i'==14 continue
	g netfe_`i' = netfe==`i'
}

* Scaling 

replace padj = sqrt(padj / 100)
replace iopt = iopt / 10

* Cleaning 

keep year gra hhid idhh planid padj choice metal mtlfe issuername issfe issfe_* netname netfe ///
pmtlord hassub fam age3 hasi iplan itenplan iexp iopt icatplan dcat pshock mtl_* netfe_*
rename (itenplan icatplan) (iten icat)
order year gra idhh planid choice iplan padj choice 
sort gra idhh year 

* Save

save "proc\analysis_i2", replace

/// CHOICE MODELS ///

* Load 

use "proc\analysis_i2", clear
global inert c.iplan#c.iten c.iplan#c.iexp c.iplan#c.iopt // add pshock?
global conv  tolerance(0.001) ltolerance(0.001) nrtolerance(0.001) diff

* DSC Logit

di "$S_TIME" // status, age
alogit choice padj iplan $inert mtl_* netfe_* if hasi==1 & dcat==1, /// might also test inertia for population that terminated
group(hhid) def(iplan) avars(padj pshock iten iexp mtl_* netfe_*) model(dsc) diff

di "$S_TIME"
clogit choice padj iplan $inert mtl_* netfe_* if dcat!=1 | (dcat==1 & hasi==0), group(hhid) diff

di "$S_TIME"

* Conditional Logit

clogit choice c.padj c.iplan $inert i.mtlfe i.netfe if gra==16, group(hhid)
est sto pcnt
est save "proc/pcnt", replace
clogit choice i.pmtlord c.iplan $inert i.mtlfe i.netfe if gra==16, group(hhid)
est sto pord
est save "proc/pord", replace
esttab pcnt pord, drop(*netfe) b(%9.3f) nostar nobase wide stats(N ll r2_p)
unique hhid if gra==16

* Conditional Logit without Auto Enrollment

clogit choice c.padj c.iplan $inert i.mtlfe i.netfe if gra==16, group(hhid)
est sto pcnt
est save "proc/pcnt", replace
clogit choice i.pmtlord c.iplan $inert i.mtlfe i.netfe if gra==16, group(hhid)
est sto pord
est save "proc/pord", replace
esttab pcnt pord, drop(*netfe) b(%9.3f) nostar nobase wide stats(N ll r2_p)

* Mixed Logit

eststo mlpcnt: asmixlogit choice padj       $inert i.mtlfe i.netfe if gra==16, case(hhid) alt(planid) nocons random(iplan)
est save "proc/mlpcnt", replace
eststo mlpord: asmixlogit choice i.pmtlord  $inert i.mtlfe i.netfe if gra==16, case(hhid) alt(planid) nocons random(iplan)
est save "proc/mlpord", replace

esttab ml*, drop(*netfe) b(%9.3f) nostar nobase wide stats(N ll r2_p)
