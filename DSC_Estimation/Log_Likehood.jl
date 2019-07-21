# using NLopt
using ForwardDiff

# Calculate Log Likelihood
function log_likelihood(d::InsuranceLogit,p::parDict{T}) where T
    individual_values!(d,p)
    individual_shares(d,p)

    ll = mean(p.L_i)
    return ll
end


function log_likelihood!(grad::Vector{Float64},
                            d::InsuranceLogit,p::parDict{T}) where T
    Q = d.parLength[:All]
    N = size(d.draws,1)
    grad[:] .= 0.0
    ll = 0.0
    Pop = length(p.L_i)


    individual_values!(d,p)
    # individual_shares(d,p)

    # ll = mean(p.L_i)

    for app in eachperson(d.data)
        ll+= ll_obs!(grad,app,d,p)
    end
    if isnan(ll)
        ll = -1e20
    end
    for q in 1:Q
        grad[q]=grad[q]/Pop
    end
    return ll/Pop
end



function log_likelihood!(hess::Matrix{Float64},grad::Vector{Float64},
                            d::InsuranceLogit,p::parDict{T}) where T
    Q = d.parLength[:All]
    N = size(d.draws,1)
    hess[:] .= 0.0
    grad[:] .= 0.0
    ll = 0.0
    Pop = length(p.L_i)
    # Pop =sum(weight(d.data).*choice(d.data))

    individual_values!(d,p)
    # individual_shares(d,p)

    for app in eachperson(d.data)
        ll_obs = ll_obs!(hess,grad,app,d,p)
        ll+=ll_obs
    end

    if isnan(ll)
        ll = -1e20
    end

    for q in 1:Q
        grad[q]=grad[q]/Pop
        for r in 1:Q
            hess[q,r]=hess[q,r]/Pop
        end
    end

    return ll/Pop
end


function log_likelihood!(thD::Array{Float64,3},
                            hess::Matrix{Float64},grad::Vector{Float64},
                            d::InsuranceLogit,p::parDict{T};cont_flag::Bool=false) where T
    Q = d.parLength[:All]
    N = size(d.draws,1)

    thD[:] .= 0.0
    hess[:] .= 0.0
    grad[:] .= 0.0
    ll = 0.0
    Pop =sum(weight(d.data).*choice(d.data))

    individual_values!(d,p)
    individual_shares(d,p)

    for app in eachperson(d.data)
        ll_obs,pars_relevant = ll_obs!(thD,hess,grad,app,d,p)
        ll+=ll_obs
    end

    if isnan(ll)
        ll = -1e20
    end

    for q in 1:Q
        grad[q]=grad[q]/Pop
        for r in 1:Q
            hess[q,r]=hess[q,r]/Pop
            for t in 1:Q
                thD[q,r,t]=thD[q,r,t]/Pop
            end
        end
    end

    return ll/Pop
end

function log_likelihood(d::InsuranceLogit,p::Array{T,1},x::Vector{Float64}) where T
    p_vec = Vector{T}(undef,length(x))
    p_vec[:] = x[:]
    p_vec[1:length(p)] = p[:]
    params = parDict(d,p_vec)
    ll = log_likelihood(d,params)
    return ll
end

function log_likelihood(d::InsuranceLogit,p::Array{T}) where T
    params = parDict(d,p)
    ll = log_likelihood(d,params)
    return ll
end

function log_likelihood!(grad::Vector{Float64},d::InsuranceLogit,p::Array{T}) where T
    params = parDict(d,p)
    ll = log_likelihood!(grad,d,params)
    return ll
end

function log_likelihood!(hess::Matrix{Float64},grad::Vector{Float64},
                            d::InsuranceLogit,p::Array{T}) where T
    params = parDict(d,p)
    ll = log_likelihood!(hess,grad,d,params)
    return ll
end

function log_likelihood!(thD::Array{Float64,3},
                            hess::Matrix{Float64},grad::Vector{Float64},
                            d::InsuranceLogit,p::Array{T}) where T
    params = parDict(d,p)
    ll = log_likelihood!(thD,hess,grad,d,params)
    return ll
end



# Calculate Standard Errors
# Hiyashi, p. 491
function calc_Avar(d::InsuranceLogit,p::parDict{T}) where T
    # Calculate μ_ij, which depends only on parameters
    individual_values!(d,p)
    individual_shares(d,p)

    Σ = zeros(d.parLength[:All],d.parLength[:All])
    Pop =sum(weight(d.data).*choice(d.data))
    grad_obs = Vector{Float64}(undef,d.parLength[:All])

    for app in eachperson(d.data)
        grad_obs[:].=0
        ll_obs,pars_relevant = ll_obs!(grad_obs,app,d,p)
        S_n = grad_obs*grad_obs'
        Σ+= S_n
    end

    Σ = Σ./Pop
    # This last line is correct
    AsVar = inv(Σ)
    return AsVar
end
