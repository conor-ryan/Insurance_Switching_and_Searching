/// SETUP ///

clear *
cd "M:\My Documents\CovCAInertia"

/// INERTIA CHECKS ///

* Raw File Checks

use "proc\ca_enroll_ind", clear
order idhh idind year act_pass enr_strt cov_strt cov_end planid
sort idhh year idind // act_pass is consistent across HH/years

* HH Checks

use "proc\ca_enroll_hh_wi", clear
sort idhh year

replace act_pass = "X" if act_pass=="" | year==2014
g new = idhh!=idhh[_n-1]
g new100 = new*100
g disenr = cov_end!=".."

order idhh year act_pass disenr enr_strt cov_strt cov_end planid

count if disenr==0 & idhh==idhh[_n-1] & act_pass=="X" & year!=2014 // Returning enrollees that didn't disenroll with missing active/passive values

table year act_pass, c(m new100) f(%9.1f) // A, M categories are WAY too small

tab year act_pass
tab act_pass if icatiss==1

/// HH LEVEL DESCRIPTIVES ///

* Table 1: Desc Stats of CCA Plans

use "proc\ca_enroll_ind", clear
tab year // # enrollees

use "proc/planchars", clear

replace prem = prem * 1.024 if year < 2015 // medical CPI
replace prem = prem * 1.026 if year < 2016
replace prem = prem * 1.038 if year < 2017
replace prem = prem * 1.025 if year < 2018
table year, c(m prem) f(%9.2f)

forvalues y = 2014/2018 {
	unique planid if year==`y'
	unique netid if year==`y'
	unique issuerid if year==`y'
}

use "proc\ca_enroll_hh_wi", clear
tab issuername year

* Table 2: Mean HH Chars

use "proc\analysis_i", clear
duplicates drop year idhh, force
table year, c(m iopt sd iopt)

use "proc\ca_enroll_hh_wi", clear

g icatnew  = icatplan==1
g icatstay = icatplan==2
g icatswch = icatplan==3
g icatdisc = icatplan==4
g hasinert = icatplan==2 | icatplan==3

foreach x of varlist prem p_post p_pre {
	replace `x' = `x' * 1.024 if year < 2015 // medical CPI
	replace `x' = `x' * 1.026 if year < 2016
	replace `x' = `x' * 1.038 if year < 2017
	replace `x' = `x' * 1.025 if year < 2018
}

table year, c(m fage m fpl m fsize) f(%9.2fc) row
table year, c(m icatnew m icatstay m icatswch m icatdisc m hasinert) f(%9.4f) row
table year, c(m iexp m itenplan m p_pre m p_post) f(%9.2f) row

table year, c(m prem m p_post) f(%9.2f) row 
table year if hasinert==0, c(m prem m p_post) f(%9.2f) row 
table year if hasinert==1, c(m prem m p_post) f(%9.2f) row 

collapse p_post, by(year hasinert)
lab var hasinert "Has Inertia"
lab def hasinert 0 "No" 1 "Yes"
lab val hasinert hasinert
set scheme plottig
twoway line p_post year if hasinert==1 || line p_post year if hasinert==0, ytitle("Post-Subsidy Monthly Premiums ($)", size(med)) ///
xtitle("Year", size(med)) legend(pos(6) title("Has Inertia") label(1 "No") label(2 "Yes") rows(1))

* Table 3: Premium Sensitivity by Cohort

use idhh year prem pmtlord metal using "proc\ca_enroll_hh_wi", clear
replace prem = prem * 1.024 if year < 2015 // medical CPI
replace prem = prem * 1.026 if year < 2016
replace prem = prem * 1.038 if year < 2017
replace prem = prem * 1.025 if year < 2018
egen cohort = min(year), by(idhh)
g mtlmin = pmtlord==1
replace mtlmin = mtlmin*100
table cohort year, c(m prem)   f(%9.2f)
table cohort year, c(m mtlmin) f(%9.2f)
keep if regexm(metal,"Silver")==1
table cohort year, c(m mtlmin) f(%9.2f)

* Appendix Tables 1 and 2: Market Shares

use "proc\ca_enroll_hh_wi", clear

tab issuername year
tab netname year

* Figure: Premium Densities by Inertia

use "proc\ca_enroll_hh_wi", clear

kdensity prem if (icatplan==2 | icatplan==3) || ///
kdensity prem if (icatplan==1 | icatplan==4)

* Figure: Market Shares of Big 4

use "proc\ca_enroll_hh_wi", clear

cap g hasinert = icatplan==2 | icatplan==3

g issanth = issuername=="Anthem"
g issbcbs = issuername=="Blue Shield"
g isshltn = issuername=="HealthNet"
g isskais = issuername=="Kaiser"
g issoth  = issanth==0 & issbcbs==0 & isskais==0 & isshltn==0

collapse issanth issbcbs isshltn isskais issoth, by(hasinert year)

g istr = "Inertia" 
replace istr = "No Inertia" if hasinert==0

set scheme plottig

graph bar issanth issbcbs isshltn isskais issoth if istr=="No Inertia", over(year) stack ///
saving(proc/ino.gph, replace) title("Enrollees without Inertia") ytitle("Market Share") legend(title("Insurer") rows(1) pos(6) ///
label(1 "Anthem") label(2 "BCBS") label(3 "HealthNet") label(4 "Kaiser") label(5 "Other"))

graph bar issanth issbcbs isshltn isskais issoth if istr=="Inertia", over(year) stack ///
saving(proc/iyes.gph, replace) title("Enrollees with Inertia") ytitle("Market Share") legend(title("Insurer") rows(1) pos(6) ///
label(1 "Anthem") label(2 "BCBS") label(3 "HealthNet") label(4 "Kaiser") label(5 "Other"))

grc1leg proc/ino.gph proc/iyes.gph

/// CHAMBERLAIN AND HARVEST TESTS ///

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
egen netfe = group(netid)
egen issfe = group(issuerid)
egen typfe = group(plantype)

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
	qui margins, dydx(*pmtlord) post
	eststo tmarg1
	
eststo tchamb2: qui glm s i.pmtlord i.mtlfe i.netfe##i.year if pmtlordlag!=., cluster(issfe) link(log) family(gamma)
	qui margins, dydx(*pmtlord) post
	eststo tmarg2

eststo tchamb3: qui glm s i.pmtlord i.mtlfe i.netfe  i.year i.pmtlordlag , cluster(issfe) link(log) family(gamma)
	qui margins, dydx(*pmtlord *pmtlordlag) post
	eststo tmarg3
	
eststo tchamb4: qui glm s i.pmtlord i.mtlfe i.netfe##i.year i.pmtlordlag , cluster(issfe) link(log) family(gamma)
	qui margins, dydx(*pmtlord *pmtlordlag) post
	eststo tmarg4

esttab tchamb*, keep(*pmtlord *pmtlordlag) drop(1.pmtlord 1.pmtlordlag) stats(r2 N) star(^ 0.1 * 0.05 ** 0.01 *** 0.001)
esttab tmarg*,  keep(*pmtlord *pmtlordlag) drop(1.pmtlord 1.pmtlordlag) stats(r2 N) star(^ 0.1 * 0.05 ** 0.01 *** 0.001) 

* Harvest Test

cap g lp = log(prem)

eststo harv_typ_nw: qui reg lp i.planage i.mtlfe i.typfe i.year, cluster(gra) 
eststo harv_iss_nw: qui reg lp i.planage i.mtlfe i.issfe i.year, cluster(gra) 
eststo harv_net_nw: qui reg lp i.planage i.mtlfe i.netfe i.year, cluster(gra)

eststo harv_typ_w: qui reg lp i.planage i.mtlfe i.typfe i.year [aw=planct], cluster(gra)
eststo harv_iss_w: qui reg lp i.planage i.mtlfe i.issfe i.year [aw=planct], cluster(gra)
eststo harv_net_w: qui reg lp i.planage i.mtlfe i.netfe i.year [aw=planct], cluster(gra) 

esttab harv_typ_nw harv_iss_nw harv_net_nw harv_typ_w harv_iss_w harv_net_w, keep(*planage) drop(0.planage) stats(r2 N)
