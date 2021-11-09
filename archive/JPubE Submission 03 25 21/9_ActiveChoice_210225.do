/// FINAL PREP ///

use "proc\analysis_i2", clear
keep active autodp agefe_* fam s_hassub auto_mtl gra choice year idhh issuername padj iplan inet mtlfe netfe planid

encode auto_mtl, g(mtlauto)
egen hhid = group(idhh)
drop auto_mtl idhh

g iss_anth = issuername=="Anthem"
g iss_blue = issuername=="Blue Shield"
g iss_kais = issuername=="Kaiser"
g iss_hnet = issuername=="HealthNet"

g autodp_cens = autodp
sum autodp, d
replace autodp_cens = r(p5)  if autodp < r(p5)
replace autodp_cens = r(p95) if autodp > r(p95)

g autodp_disc = .
replace autodp_disc = 1 if autodp<=-50
replace autodp_disc = 2 if autodp >-50 & autodp<=-10
replace autodp_disc = 3 if autodp >-10 & autodp<=0
replace autodp_disc = 4 if autodp >0   & autodp<=10
replace autodp_disc = 5 if autodp >10  & autodp<=50
replace autodp_disc = 6 if autodp >50


gen autodp_disc_1 = autodp<=-50
gen autodp_disc_2 = autodp >-50 & autodp<=-10
gen autodp_disc_3 = autodp >-10 & autodp<=0
gen autodp_disc_4 = autodp >0   & autodp<=10
gen autodp_disc_5 = autodp >10  & autodp<=50
gen autodp_disc_6 = autodp >50


compress

save "proc\analysis_i3", replace
outsheet using "Output\analysis_i3.csv", comma replace nolabel

/// ACTIVE SELECTION ///

use "proc\analysis_i3", clear
reghdfe  active ib3.autodp_disc agefe_1 agefe_2 fam s_hassub, a(mtlauto gra year)
reghdfe  active c.autodp agefe_1 agefe_2 fam s_hassub, a(mtlauto gra year)
logistic active c.autodp agefe_1 agefe_2 fam s_hassub i.mtlauto i.gra i.year
margins, dydx(autodp)

logistic active ib3.autodp_disc agefe_1 agefe_2 fam s_hassub i.mtlauto i.gra i.year
margins, dydx( ib3.autodp_disc agefe_1 agefe_2 fam s_hassub)

/// PLAN CHOICE ///

* Load

use "analysis_i3_720" if active==1, clear
drop autodp autodp_cens mtlauto active agefe_0

* Prep

encode issuername, g(issfe)
g issfeml = issfe
replace issfeml = 8 if issfe==1 | issfe==2 | issfe==5 | issfe==6 // Setting Molina as reference for big 4

cap egen hhidyr = group(hhid year)

foreach y in padj iplan inet {
foreach x in agefe_1 agefe_2 fam s_hassub {
	g `y'_`x' = `y' * `x'
}
}

gsort hhid year mtlfe netfe issuername planid
order hhidyr hhid year gra choice planid issuername netfe mtlfe padj iplan inet agefe_1 agefe_2 fam s_hassub iss_* padj_* iplan_* inet_*

* Cond'l Logit

clogit choice padj* iplan* inet* i.mtlfe i.netfe i.issfe, group(hhidyr)

* Mixed Logit
/*
cmset hhid year planid // chooser year alternative

cmxtmixlogit choice padj* iplan* inet* i.mtlfe i.netfe i.issfeml ///
if gra==15, random(iss_anth iss_blue iss_kais iss_hnet)
