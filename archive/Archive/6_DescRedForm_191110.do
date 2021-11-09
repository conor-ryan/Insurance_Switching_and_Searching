/// SETUP ///

clear *
cd "G:\Shared Drives\CovCAInertia"

/// TABLE I: ENROLLMENT, PLAN OFFERINGS, PREMIUMS, AND MARKET SHARES ///

* Number of Enrollees

use "proc\ca_enroll_ind", clear
tab year // # enrollees

* Number HHs

use "proc\ca_enroll_hh", clear
tab year // # HHs

* Inertia Enrollment Chars (Plan)

use "proc\ca_enroll_hh_wi", clear
g icatplan_1 = icatplan==1 // New
g icatplan_2 = icatplan==2 // Stayer
g icatplan_3 = icatplan==3 // Switcher
tabstat icatplan_1 icatplan_2 icatplan_3, by(year) s(me) f(%9.4fc)

* Choice Characteristics of Returning Households

use "proc\ca_enroll_hh_wi", clear

g act_pass_bin = 1 if act_pass=="A" // 1 auto, 0 manual, else missing
replace act_pass_bin = 0 if act_pass=="M"

g sameplanid = planid==lagplan
g sameautoid = planid==autoplanid & act_pass_bin==0
g diffautoid = planid!=autoplanid & act_pass_bin==0

tabstat sameplanid sameautoid, by(act_pass_bin) // suggests A is the auto enrollers; 97% of autos pick auto assn'd plan!

tabstat act_pass_bin sameautoid diffautoid if idhh==idhh[_n-1], by(year) f(%9.4fc)

* Health Plan Offerings

use "proc\ca_enroll_hh_wi", clear
duplicates drop year planid, force
tab year
duplicates drop year netname, force
tab year
duplicates drop year issuer, force
tab year

* Market Shares

use "proc\ca_enroll_hh_wi", clear
replace issuername = "Other Insurers" if issuername!="Anthem" & issuername!="Blue Shield" & issuername!="HealthNet" & issuername!="Kaiser"
tab issuername year, col nofreq

* Median Monthly Post-APTC Premiums

use "proc\analysis_i_191112", clear
tabstat prem padj s_aptc, s(med p25 p75) by(year) c(s) f(%9.0fc)

/// TABLE II: CHARS OF PLAN SELECTIONS ACROSS ENROLLMENT COHORTS ///

use using "proc\ca_enroll_hh_wi", clear
keep year gra idhh prem pmtlord metal

egen cohort = min(year), by(idhh)

replace prem = prem * 1.024 if year < 2015 // medical CPI
replace prem = prem * 1.026 if year < 2016
replace prem = prem * 1.038 if year < 2017
replace prem = prem * 1.025 if year < 2018

g mtlmin = pmtlord==1
replace mtlmin = 100*mtlmin

g minprem_a = prem if mtlmin==100
egen minprem_b = max(minprem_a), by(year gra metal) 

g mindiff = prem - minprem_b

table cohort year, c(m prem)    f(%9.0f) row
table cohort year, c(m mtlmin)  f(%9.2f) row
table cohort year, c(m mindiff) f(%9.2f) row

/// FIGURE I: INSURER MARKET SHARES AMONG NEW AND RETURNING HOUSEHOLDS ///

* Prep

use "proc\ca_enroll_hh_wi", clear

g returnhh = idhh==idhh[_n-1]

g issanth = issuername=="Anthem"
g issbcbs = issuername=="Blue Shield"
g isshltn = issuername=="HealthNet"
g isskais = issuername=="Kaiser"
g issoth  = issanth==0 & issbcbs==0 & isskais==0 & isshltn==0

foreach i of varlist issanth issbcbs isshltn isskais issoth {
	replace `i' = `i' * 100
}

collapse issanth issbcbs isshltn isskais issoth, by(returnhh year) fast

g istr = "Inertia" 
replace istr = "No Inertia" if returnhh==0

* Graph

set scheme plotplain

graph bar issanth issbcbs isshltn isskais issoth if istr=="No Inertia" & year!=2014, over(year) stack ///
saving(proc/ino.gph, replace) title("New Households") ytitle("Market Share (%)") legend(title("Insurer") rows(1) pos(6) ///
label(1 "Anthem") label(2 "BCBS") label(3 "HealthNet") label(4 "Kaiser") label(5 "Other"))

graph bar issanth issbcbs isshltn isskais issoth if istr=="Inertia" & year!=2014, over(year) stack ///
saving(proc/iyes.gph, replace) title("Returning Households") ytitle("Market Share (%)") legend(title("Insurer") rows(1) pos(6) ///
label(1 "Anthem") label(2 "BCBS") label(3 "HealthNet") label(4 "Kaiser") label(5 "Other"))

grc1leg proc/ino.gph proc/iyes.gph

/// TABLE III: SENSITIVITY TO DEFAULT AND RIVAL PLAN PREMIUMS ///

* Load

use "proc\analysis_i_191112", clear 
keep if icatplan!=1 // No New Enrollees

* Prep

egen grafe = group(gra)

egen temp = mean(padj), by(year idhh)
g undp = (iopt*temp - auto_p_post) / (iopt - 1)
lab var undp "Mean of Non-Default Premiums"
drop temp

g switch = (netname!=autonet | metal!=lagmetal) & planchoice!=autoplanid
lab var switch "HH Switched Plans"

g active = .
replace active = 1 if act_pass=="M"
replace active = 0 if act_pass!="M" // or =="M"?
lab var active "HH Made Active Selection"

keep if choice==1

replace auto_p_post = auto_p_post / 100
replace undp        = undp        / 100
replace autodp      = autodp      / 100

replace switch = switch * 100
replace active = active * 100

egen mtlfe = group(lagmetal)
egen netfe = group(lagnet)

* Regression

global dem i.age3 fam s_hassub
global fe i.mtlfe i.year // i.netfe

eststo s1: qui areg switch auto_p_post undp $dem $fe, a(gra) cluster(idhh)
eststo s2: qui areg switch autodp      undp $dem $fe, a(gra) cluster(idhh)

eststo a1: qui areg active auto_p_post undp $dem $fe, a(gra) cluster(idhh)
eststo a2: qui areg active autodp      undp $dem $fe, a(gra) cluster(idhh)

esttab s* a*, order(auto_p_post autodp undp 1.age3 2.age3 fam s_hassub) stats(r2 N) drop(0.age3 *.mtlfe *.year) 

/// APPENDIX TABLE AII: CHAMBERLAIN TESTS ///

* Enrollment Shares

use planid gra year metal using "proc\ca_enroll_hh_wi", clear
bys year gra planid: g planct = _N
bys year gra: g gract = _N
duplicates drop planid gra year, force
g s  = planct / gract
g ls = log(s) 
save proc/planshares, replace

* Merge with Plan Data

use proc/planchars, clear

replace prem = prem * 1.024 if year < 2015 // medical CPI
replace prem = prem * 1.026 if year < 2016
replace prem = prem * 1.038 if year < 2017
replace prem = prem * 1.025 if year < 2018
replace prem = prem / 100 // scaling

egen mtlfe = group(metal) // FE
egen netfe = group(netname)
egen issfe = group(issuerid)

sort  gra planid year
g premlag = prem[_n-1] if planid==planid[_n-1] // Lagged Premium
g premorig = premlag                           // Original Premium
replace premorig = prem[_n-2] if planid==planid[_n-2]
replace premorig = prem[_n-3] if planid==planid[_n-3]
replace premorig = prem[_n-4] if planid==planid[_n-4]
order gra planid year prem premlag premorig
corr prem premlag // .9734, .9042
corr prem premorig

sort gra planid year
g pmtlordlag  = pmtlord[_n-1] if planid==planid[_n-1]
g pmtlordorig = pmtlordlag
replace pmtlordorig = pmtlord[_n-2] if planid==planid[_n-2]
replace pmtlordorig = pmtlord[_n-3] if planid==planid[_n-3]
replace pmtlordorig = pmtlord[_n-4] if planid==planid[_n-4]
order year gra metal prem pmtlord*
spearman pmtlord pmtlordlag, stats(rho p)
bys year: spearman pmtlord pmtlordlag, stats(rho p)
corr pmtlord pmtlordlag  // 0.74
corr pmtlord pmtlordorig // 0.57

order planid gra year
sort  planid gra year
g planage = 0
forvalues y = 1/4 {
	replace planage = `y' if planid==planid[_n-`y'] & gra==gra[_n-`y']
}
order planid gra year planage

merge 1:1 year gra planid using proc/planshares, keep(match) nogen // 20 / 2650 mismatched
save proc/chambtest, replace

* Chamberlain Test

eststo tchamb1: qui glm s i.pmtlord i.mtlfe i.netfe  i.year if pmtlordlag!=., cluster(issfe) link(log) family(gamma)
	qui margins, eydx(*pmtlord) post
	eststo tmarg1
	
eststo tchamb2: qui glm s i.pmtlord i.mtlfe i.netfe##i.year if pmtlordlag!=., cluster(issfe) link(log) family(gamma)
	qui margins, eydx(*pmtlord) post
	eststo tmarg2

eststo tchamb3: qui glm s i.pmtlord i.mtlfe i.netfe  i.year i.pmtlordlag , cluster(issfe) link(log) family(gamma)
	qui margins, eydx(*pmtlord *pmtlordlag) post
	eststo tmarg3
	
eststo tchamb4: qui glm s i.pmtlord i.mtlfe i.netfe##i.year i.pmtlordlag , cluster(issfe) link(log) family(gamma)
	qui margins, eydx(*pmtlord *pmtlordlag) post
	eststo tmarg4

esttab tchamb*, keep(*pmtlord *pmtlordlag) drop(1.pmtlord 1.pmtlordlag) stats(N) star(^ 0.1 * 0.05 ** 0.01 *** 0.001)
esttab tmarg*,  keep(*pmtlord *pmtlordlag) drop(1.pmtlord 1.pmtlordlag) stats(N) star(^ 0.1 * 0.05 ** 0.01 *** 0.001) tr(1*@)
