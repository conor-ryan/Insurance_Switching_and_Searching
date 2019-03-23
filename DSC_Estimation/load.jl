using DataFrames
using CSV
########################################################################
#################### Loading and Cleaning Data #########################
########################################################################
# Load the data
df = CSV.read("$(homedir())/Documents/Research/CovCAInertia/Output/analysis_i2.csv")

# No constant
df[:constant] = ones(size(df, 1))
df[:dprem][ismissing.(df[:dprem])] .= 0.0
for key in [:lag_disc,:def_padj,:def_mtl_brz,:def_mtl_cat,:def_mtl_gld,:def_mtl_plt,:def_mtl_s73,:def_mtl_s87,:def_mtl_s94]
    df[key][df[:hasi].==0] .= 0.0
end
