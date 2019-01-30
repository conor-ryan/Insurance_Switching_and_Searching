using DataFrames
using CSV
########################################################################
#################### Loading and Cleaning Data #########################
########################################################################
# Load the data
df = CSV.read("$(homedir())/Documents/Research/CovCAInertia/Output/analysis_i2.csv")

# No constant
df[:constant] = ones(size(df, 1))
