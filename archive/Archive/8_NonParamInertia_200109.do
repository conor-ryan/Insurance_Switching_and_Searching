*********************************************************************
* Non-Parametrically Identify Inertia
*********************************************************************

/// SETUP ///

clear *
cd "\\files.umn.edu\cla\home\ryan0463\My Documents\CovCAInertia"
set seed 1234 

/// Read in straight panel data ///

use "proc\ca_enroll_hh_wi", clear
sort gra idhh year


/// Define Demographic Types ///
gen age_cat = floor(fage/10)*10

gen mem_cat = fsize
replace mem_cat = 2 if mem_cat>2

gen fpl_cat = 0
replace fpl_cat = 1 if fpl<=400 & fpl>250
replace fpl_cat = 2 if fpl<=250 & fpl>=100

/// Use First Three Years of continuous HH Enrollment ///
egen enter_year = min(year), by(idhh)

drop if year>(enter_year+2)

gen cnt = 1
egen cov_yrs = sum(cnt), by(idhh)

keep if cov_yrs==3


/// Initial Plan - assume for now ID doesn't change... /// 
gen prod_name = netname + "-" + product + "-" + metal
egen prodnum = group(prod_name)

egen prod_cnt = sum(cnt), by(prod_name year gra)

gen cnt_2014_tmp = prod_cnt if year==2014
gen cnt_2015_tmp = prod_cnt if year==2015
gen cnt_2016_tmp = prod_cnt if year==2016
gen cnt_2017_tmp = prod_cnt if year==2017
gen cnt_2018_tmp = prod_cnt if year==2018

egen cnt_2014 = max(cnt_2014_tmp), by(prod_name gra)
egen cnt_2015 = max(cnt_2015_tmp), by(prod_name gra)
egen cnt_2016 = max(cnt_2016_tmp), by(prod_name gra)
egen cnt_2017 = max(cnt_2017_tmp), by(prod_name gra)
egen cnt_2018 = max(cnt_2018_tmp), by(prod_name gra)


gen init_prod_tmp = prodnum if year==enter_year
replace init_prod_tmp = 0 if year>enter_year
egen init_prod = sum(init_prod_tmp),by(idhh)
drop init_prod_tmp

// //Drop Anthem Enrollees in 2016, as the product exits//
// drop if init_prod<=10 & enter_year==2016
//
// //Drop Contra Costa, as the product exits//
// drop if (init_prod>=26 & init_prod<=29)
//
// //Drop United Initial Buyers which is small and not offered every year//
// drop if (init_prod>=75 & init_prod<=80)
// tab issuername enter_year
// tab issuername year

sort idhh year
gen switch = 0
replace switch = 1 if prodnum!=prodnum[_n-1] & idhh == idhh[_n-1]

tab switch icatplan3

gen init_choice = prodnum == init_prod
gen mkt_tenure = year-enter_year

egen init_cnt = sum(init_choice),by(idhh)
egen switch_cnt = sum(switch),by(idhh)

/// Switching Sequence 
gen switch_seq = "Unknown"
replace switch_seq = "Init - Init - Init" if init_cnt==3 & switch_cnt==0
replace switch_seq = "Init - Init - Switch" if init_cnt==2 & switch_cnt==1
replace switch_seq = "Init - Switch - Init" if init_cnt==2 & switch_cnt==2
replace switch_seq = "Init - Switch - Stay" if init_cnt==1 & switch_cnt==1
replace switch_seq = "Init - Switch - Switch" if init_cnt==1 & switch_cnt==2


save "proc/non_param_intert.dta", replace

/// Compare Frequencies ///
use "proc/non_param_intert.dta", clear

keep if year == enter_year

//Products are continually offered
drop if enter_year==2016 & (cnt_2018==. | cnt_2017==.)
drop if enter_year==2015 & (cnt_2017==. | cnt_2016==.)
drop if enter_year==2014 & (cnt_2016==. | cnt_2015==.)

//Drop Sharp Enrollees. Something Odd with the HDHP/Catastrophic plans...
drop if issuername=="Sharp"


collapse (sum) cnt, by (prod_name switch_seq age_cat mem_cat fpl_cat gra)

egen type_total = sum(cnt), by(prod_name age_cat mem_cat fpl_cat gra)
sort age_cat mem_cat gra prod_name switch_seq


gen seq_share = cnt/type_total
gen I_I_I_share = seq_share if switch_seq == "Init - Init - Init"
gen I_I_S_share = seq_share if switch_seq == "Init - Init - Switch"
gen I_S_I_share = seq_share if switch_seq == "Init - Switch - Init"
gen I_S_S_share = seq_share if switch_seq == "Init - Switch - Switch"
gen I_S_St_share = seq_share if switch_seq == "Init - Switch - Stay"


collapse (max) I_*, by(prod_name age_cat mem_cat fpl_cat gra type_total)

foreach v of varlist I_I_I_share-I_S_St_share{
	replace `v'=0 if `v'==.
}

egen test_sum = sum(type_total)
keep if type_total>=50
egen test_sum_2 = sum(type_total)

gen prob_ratio = I_S_I_share/I_I_S_share


foreach v of varlist I_I_I_share-prob_ratio{
	pctile pct_`v' = `v', nq(100)
}

gen pctile = _n

keep if _n==10 | _n==25 | _n==50 | _n==75 | _n==90 | _n==99
keep pct*

 
