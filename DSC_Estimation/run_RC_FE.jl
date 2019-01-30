using BenchmarkTools
using JLD2
using CSV
using Random
using Dates
using LinearAlgebra
using Statistics


# Data Structure
include("InsChoiceData.jl")

#Halton Draws
include("Halton.jl")

# Random Coefficients MLE
include("RandomCoefficients.jl")
include("RandomCoefficients_der.jl")
include("Log_Likehood.jl")
include("Estimate_Basic.jl")
include("utility.jl")
println("Code Loaded")

# Load the Data
include("load.jl")
println("Data Loaded")


df_small = df[df[:gra].==16,:]
df = 0.0

# Structre the data
c = ChoiceData(df_small;
    person = [:hh_year_id],
    product = [:product],
    choice = [:choice],
    prodchars=Vector{Symbol}(undef,0),
    prodchars_0=Vector{Symbol}(undef,0),
    demoRaw=Vector{Symbol}(undef,0),
    fixedEffects=[:product],
    wgt=[:constant])

# Fit into model
m = InsuranceLogit(c,1)
println("Data Loaded")

#γ0start = rand(1)-.5
# γstart = rand(m.parLength[:γ])/10 .-.05

β0start = rand(m.parLength[:β])/10 .-.05
βstart = rand(m.parLength[:γ])/10 .- .05
σstart = rand(m.parLength[:σ])/10 .- .05
FEstart = rand(m.parLength[:FE])/100 .-.005

p0 = vcat(β0start,βstart,σstart,FEstart)
par0 = parDict(m,p0)
#
individual_values!(m,par0)
individual_shares(m,par0)

ll = log_likelihood(m,par0)

# println("Gradient Test")
# grad_2 = Vector{Float64}(undef,length(p0))
# hess_2 = Matrix{Float64}(undef,length(p0),length(p0))
# ll = log_likelihood!(hess_2,grad_2,m,par0)
#
#
# f_obj(x) = log_likelihood(m,x)
# grad_1 = Vector{Float64}(undef,length(p0))
# hess_1 = Matrix{Float64}(undef,length(p0),length(p0))
# fval_old = f_obj(p0)
# # # #
# println("Grad")
# ForwardDiff.gradient!(grad_1,f_obj, p0)#, cfg)
# println("Hessian")
# cfg = ForwardDiff.HessianConfig(f_obj, p0, ForwardDiff.Chunk{3}())
# ForwardDiff.hessian!(hess_1,f_obj, p0, cfg)
#
# println(fval_old-ll)
# println(maximum(abs.(grad_1-grad_2)))
# println(maximum(abs.(hess_1-hess_2)))



est = newton_raphson_ll(m,p0)




# using Profile
# Profile.init(n=10^8,delay=.001)
# Profile.clear()
# #Juno.@profile add_obs_mat!(hess,grad,hess_obs,grad_obs,Pop)
# # Juno.@profile log_likelihood!(thD_2,hess_2,grad_2,m,par0)
# Juno.@profile res = GMM_objective(m,p0,W)
# Juno.profiletree()
# Juno.profiler()
