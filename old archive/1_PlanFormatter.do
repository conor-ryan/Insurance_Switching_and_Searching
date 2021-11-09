* Directory

cd "H:\My Drive\Dissertation\Paper3"

/// CCA PLAN CHARACTERISTICS ///

* Master List of Issuer Names from 2016 File

/*
Anthem
Blue Shield
Chinese
Contra Costa
HealthNet
Kaiser
LA Care
Molina
Oscar
Sharp
United
Valley
Western 
*/

* Networks 2016-2017

import delimited using data\pufcaplan2017.csv, clear
keep standardcomponentid hiosproductid networkid plantype
rename (standardcomponentid hiosproductid networkid) (planid prodid netidorig)
g issid = substr(planid,1,5)
g netid = substr(prodid,1,5) + substr(netid,-3,3)
order netid prodid planid
keep  netid prodid planid
sort  netid prodid planid
duplicates drop planid, force
keep netid prodid
duplicates drop prodid, force
save proc/netxwalk17, replace

import delimited using data\pufcaplan2016.csv, clear
keep standardcomponentid hiosproductid networkid
rename (standardcomponentid hiosproductid networkid) (planid prodid netidorig)
g issid = substr(planid,1,5)
g netid = substr(prodid,1,5) + substr(netid,-3,3)
order netid prodid planid
keep  netid prodid planid
sort  netid prodid planid
duplicates drop planid, force
keep netid prodid
duplicates drop prodid, force
save proc/netxwalk16, replace

* Plan Characteristics from HIX Compare 2015

use st planmarket year planid area carrier planname metal plantype premi27 using ///
"H:\My Drive\Data\HIX Compare\plans_2015" if st=="CA" & planmarket!=2, clear
destring area, g(gra) ignore("CA")
g prem = premi27 / 1.048
rename carrier issuername
replace metal = "Minimum Coverage" if metal=="Catastrophic"
replace metal = "Bronze HDHP" if regexm(planname,"HSA")==1
rename  plantype plantype_num
g plantype = ""
replace plantype = "PPO" if plantype_num==1
replace plantype = "HMO" if plantype_num==2
replace plantype = "POS" if plantype_num==3
replace plantype = "EPO" if plantype_num==4
drop st planmarket area premi27 plantype_num
replace issuername = "Anthem"      if issuername=="Anthem Blue Cross"
replace issuername = "Blue Shield" if issuername=="Blue Shield of California"
replace issuername = "Chinese"     if issuername=="Chinese Community Health Plan"
replace issuername = "HealthNet"   if issuername=="Health Net of California, Inc"
replace issuername = "Kaiser"      if issuername=="Kaiser Permanente"
replace issuername = "LA Care"     if issuername=="L,A, Care Health Plan"
replace issuername = "Molina"      if issuername=="MOLINA HEALTHCARE OF CALIFORNIA"
replace issuername = "Moda"        if issuername=="Moda Health Insurance" // What is this?
replace issuername = "Valley"      if issuername=="Valley Health Plan"
replace issuername = "Western"     if issuername=="Western Health Advantage"
save proc/hix15ca, replace // Bad way to do this by itself... could be used to augment

* Plan Characteristics

import excel using data\prem2018.xlsx, first case(lower) cellrange(A8:K25916) clear
rename (ratingareaid individualrate applicant hoisid metallevel network product fullplanname) ///
(gra prem issuername issuerid metal plantype product planname)
keep if age=="21"
g year = 2018
destring gra, ignore("Rating Area ") replace
keep  year gra issuerid planid issuername prem metal plantype
order year gra issuerid planid issuername prem metal plantype
save proc/planchars18, replace // Not going to worry about this 'til I have 2018 enrollment

import excel using data\prem2017.xlsx, first case(lower) cellrange(A7:K27929) clear
rename (ratingareaid individualrate applicant hoisid metallevel network product fullplanname) ///
(gra prem issuername issuerid metal plantype product planname)
keep if age=="21"
g year = 2017
destring gra, ignore("Rating Area ") replace
keep  year gra issuerid planid issuername prem metal plantype
order year gra issuerid planid issuername prem metal plantype
g prodid = substr(planid,1,10)
merge m:1 prodid using proc/netxwalk17, keep(master match) nogen
replace issuername = "Chinese" if issuername=="Chinese Community"
replace issuername = "LA Care" if issuername=="L.A. Care"
save proc/planchars17, replace

import excel using data\prem2016.xlsx, first case(lower) cellrange(A8:K25354) clear
rename (ratingareaid individualrate applicant hiosid metallevel network product fullplanname) ///
(gra prem issuername issuerid metal plantype product planname)
keep if age=="21"
g year = 2016
destring gra, ignore("Rating Area ") replace
keep  year gra issuerid planid issuername prem metal plantype
order year gra issuerid planid issuername prem metal plantype
g prodid = substr(planid,1,10)
merge m:1 prodid using proc/netxwalk16, keep(master match) nogen
replace issuername = "Chinese" if issuername=="Chinese C."
replace issuername = "LA Care" if issuername=="L.A. Care"
save proc/planchars16, replace

* 2015

import excel using data\prem2015.xlsx, first case(lower) cellrange(A7:G22593) clear

rename (ratingareaid individualrate metallevel healthinsurancecompany) (gra prem metal issuername)
keep if age=="21"
g year = 2015
destring gra, ignore("Rating Area ") replace
g plantype   = substr(planname,-3,3)
keep  year gra issuername prem metal plantype
order year gra issuername prem metal plantype

replace issuername = "Chinese" if issuername=="Chinese Community"
replace issuername = "LA Care" if issuername=="L.A. Care"

replace metal = "Bronze HDHP" if metal=="HSA Bronze"

replace plantype = "HMO" if issuername=="Kaiser"
replace plantype = "HMO" if issuername=="HealthNet" & gra>=15
replace plantype = "EPO" if issuername=="Anthem" & plantype=="HSA" & (gra==4 | gra==15 | gra==16 | gra==18 | gra==19)
replace plantype = "PPO" if issuername=="Anthem" & plantype=="HSA" & (gra!=4 | gra!=15 | gra!=16 | gra!=18 | gra!=19)

/*
replace plantype = "PPO" if issuername=="Anthem" & plantype=="HSA" & gra==1
replace plantype = "PPO" if issuername=="Anthem" & plantype=="HSA" & gra==2
replace plantype = "HMO" if issuername=="Anthem" & plantype=="HSA" & gra==3
replace plantype = "EPO" if issuername=="Anthem" & plantype=="HSA" & gra==4
replace plantype = "PPO" if issuername=="Anthem" & plantype=="HSA" & gra==5
replace plantype = "PPO" if issuername=="Anthem" & plantype=="HSA" & gra==6
replace plantype = "HMO" if issuername=="Anthem" & plantype=="HSA" & gra==7
replace plantype = "PPO" if issuername=="Anthem" & plantype=="HSA" & gra==8
replace plantype = "PPO" if issuername=="Anthem" & plantype=="HSA" & gra==9
replace plantype = "PPO" if issuername=="Anthem" & plantype=="HSA" & gra==10
replace plantype = "HMO" if issuername=="Anthem" & plantype=="HSA" & gra==11
replace plantype = "PPO" if issuername=="Anthem" & plantype=="HSA" & gra==12
replace plantype = "PPO" if issuername=="Anthem" & plantype=="HSA" & gra==13
replace plantype = "PPO" if issuername=="Anthem" & plantype=="HSA" & gra==14
replace plantype = "HMO" if issuername=="Anthem" & plantype=="HSA" & gra==15
replace plantype = "HMO" if issuername=="Anthem" & plantype=="HSA" & gra==16
replace plantype = "HMO" if issuername=="Anthem" & plantype=="HSA" & gra==17
replace plantype = "HMO" if issuername=="Anthem" & plantype=="HSA" & gra==18
replace plantype = "HMO" if issuername=="Anthem" & plantype=="HSA" & gra==19
*/
duplicates drop issuername gra plantype metal, force // dropping 3 sharp plans

save proc/prem15, replace

* 2014

import excel using data\prem2014.xlsx, first case(lower) cellrange(A8:E504) clear

rename (producttype ratingregion healthinsurancecompany monthlyprice metaltier) (plantype gra issuername prem metal)
g year = 2014
keep  year gra issuername prem metal plantype
order year gra issuername prem metal plantype

replace issuername="Anthem"       if issuername=="Anthem Blue Cross"
replace issuername="Blue Shield"  if issuername=="Blue Shield of California"
replace issuername="Chinese"      if issuername=="Chinese Community Health Plan"
replace issuername="Contra Costa" if issuername=="Contra Costa Health Plan"
replace issuername="HealthNet"    if issuername=="Health Net"
replace issuername="Kaiser"       if issuername=="Kaiser Permanente"
replace issuername="Western"      if issuername=="Western Health Advantage"

bys gra issuername metal plantype: g dups  = _N
bys gra issuername metal plantype: g dupno = _n 
sort gra issuername metal plantype prem
replace metal = "Bronze HDHP" if dups==2 & dupno==1 & metal=="Bronze"
drop dups dupno

duplicates drop gra issuername metal plantype, force // 3 Sharp plans

save proc\prem14, replace
