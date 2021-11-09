/* SET WORKING DIRECTORY */

clear *
cap cd "G:\Shared Drives\CovCAInertia"
cap cd "\\files.umn.edu\cla\home\ryan0463\My Documents\CovCAInertia"

/* SETTINGS */

set scheme plotplain

/* SCRIPTS */

do "code\1_EnrFormatter_210225"
do "code\2_Premiums_210225"
do "code\3_InertXwalk_210811"
do "code\4_InertCat_210811"
do "code\5_ChoicePrep_210225"
do "code\6_DescRedForm_210225"
do "code\7_AnalysisPrep_210812"