using DataFrames
using CSV
########################################################################
#################### Loading and Cleaning Data #########################
########################################################################
# Load the data
df = CSV.read("$(homedir())/Documents/Research/CovCAInertia/Output/analysis_i2.csv";copycols=true)

# No constant
df[:constant] = ones(size(df, 1))
df[:autodp][ismissing.(df[:autodp])] .= 0.0
for key in [:def_mtl_brz,:def_mtl_cat,:def_mtl_gld,:def_mtl_plt,:def_mtl_slv,:def_mtl_hdp]
    df[key][df[:autoelig].==0] .= 0.0
end
