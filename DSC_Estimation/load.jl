using DataFrames
using CSV
########################################################################
#################### Loading and Cleaning Data #########################
########################################################################
# Load the data
df = CSV.read("$home/Output/analysis_i2.csv";copycols=true)

# No constant
df[:constant] = ones(size(df, 1))
df[:autodp][ismissing.(df[:autodp])] .= 0.0
df[:choiceset_exceed_thresh_95][ismissing.(df[:choiceset_exceed_thresh_95])] .= 0.0
df[:choiceset_exceed_thresh_98][ismissing.(df[:choiceset_exceed_thresh_98])] .= 0.0
for key in [:def_mtl_brz,:def_mtl_cat,:def_mtl_gld,:def_mtl_plt,:def_mtl_slv,:def_mtl_hdp]
    df[key][df[:autoelig].==0] .= 0.0
end

df_active = CSV.read("$home/Output/analysis_i2_active.csv";copycols=true)

# No constant
df_active[:constant] = ones(size(df_active, 1))
df_active[:autodp][ismissing.(df_active[:autodp])] .= 0.0
for key in [:def_mtl_brz,:def_mtl_cat,:def_mtl_gld,:def_mtl_plt,:def_mtl_slv,:def_mtl_hdp]
    df_active[key][df_active[:autoelig].==0] .= 0.0
end
