function calc_derSums_xy!(n::Int64,s_n::Vector{T},
                    X_mat::Matrix{Float64},Y::Matrix{Float64},
                    μ_ij::Matrix{T},
                    μ_ij_sums_n::T) where T

    dμ_ij_x_sums = 0.0
    dμ_ij_y_sums = 0.0
    dμ_ij_xy_sums = 0.0

    K = length(s_n)
    for k in 1:K
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[n,k]
        @inbounds @fastmath y = Y[n,k]

        @inbounds @fastmath s_n[k] = u/μ_ij_sums_n
        @inbounds @fastmath dμ_ij_x_sums+= u*x
        @inbounds @fastmath dμ_ij_y_sums+= u*y
        @inbounds @fastmath dμ_ij_xy_sums+= u*x*y
    end

    return dμ_ij_x_sums, dμ_ij_y_sums, dμ_ij_xy_sums
end

function calc_derSums_xy!(n::Int64,s_n::Vector{T},
                    X_mat::Matrix{Float64},Y::Vector{Float64},
                    μ_ij::Matrix{T},
                    μ_ij_sums_n::T) where T

    dμ_ij_x_sums = 0.0
    dμ_ij_y_sums = 0.0
    dμ_ij_xy_sums = 0.0

    K = length(s_n)
    for k in 1:K
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[n,k]
        @inbounds @fastmath y = Y[k]

        @inbounds @fastmath s_n[k] = u/μ_ij_sums_n
        @inbounds @fastmath dμ_ij_x_sums+= u*x
        @inbounds @fastmath dμ_ij_y_sums+= u*y
        @inbounds @fastmath dμ_ij_xy_sums+= u*x*y
    end

    return dμ_ij_x_sums, dμ_ij_y_sums, dμ_ij_xy_sums
end

function calc_derSums_xy!(n::Int64,s_n::Vector{T},
                    X_mat::Matrix{Float64},y::Float64,
                    μ_ij::Matrix{T},
                    μ_ij_sums_n::T) where T

    dμ_ij_x_sums = 0.0
    dμ_ij_y_sums = 0.0
    dμ_ij_xy_sums = 0.0

    K = length(s_n)
    for k in 1:K
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[n,k]

        @inbounds @fastmath s_n[k] = u/μ_ij_sums_n
        @inbounds @fastmath dμ_ij_x_sums+= u*x
        @inbounds @fastmath dμ_ij_y_sums+= u*y
        @inbounds @fastmath dμ_ij_xy_sums+= u*x*y
    end

    return dμ_ij_x_sums, dμ_ij_y_sums, dμ_ij_xy_sums
end

function calc_derSums_xy!(n::Int64,s_n::Vector{T},
                    X_mat::Vector{Float64},Y::Matrix{Float64},
                    μ_ij::Matrix{T},
                    μ_ij_sums_n::T) where T

    dμ_ij_x_sums = 0.0
    dμ_ij_y_sums = 0.0
    dμ_ij_xy_sums = 0.0

    K = length(s_n)
    for k in 1:K
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[k]
        @inbounds @fastmath y = Y[n,k]

        @inbounds @fastmath s_n[k] = u/μ_ij_sums_n
        @inbounds @fastmath dμ_ij_x_sums+= u*x
        @inbounds @fastmath dμ_ij_y_sums+= u*y
        @inbounds @fastmath dμ_ij_xy_sums+= u*x*y
    end

    return dμ_ij_x_sums, dμ_ij_y_sums, dμ_ij_xy_sums
end

function calc_derSums_xy!(n::Int64,s_n::Vector{T},
                    X_mat::Vector{Float64},Y::Vector{Float64},
                    μ_ij::Matrix{T},
                    μ_ij_sums_n::T) where T

    dμ_ij_x_sums = 0.0
    dμ_ij_y_sums = 0.0
    dμ_ij_xy_sums = 0.0

    K = length(s_n)
    for k in 1:K
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[k]
        @inbounds @fastmath y = Y[k]

        @inbounds @fastmath s_n[k] = u/μ_ij_sums_n
        @inbounds @fastmath dμ_ij_x_sums+= u*x
        @inbounds @fastmath dμ_ij_y_sums+= u*y
        @inbounds @fastmath dμ_ij_xy_sums+= u*x*y
    end

    return dμ_ij_x_sums, dμ_ij_y_sums, dμ_ij_xy_sums
end

function calc_derSums_xy!(n::Int64,s_n::Vector{T},
                    X_mat::Vector{Float64},y::Float64,
                    μ_ij::Matrix{T},
                    μ_ij_sums_n::T) where T

    dμ_ij_x_sums = 0.0
    dμ_ij_y_sums = 0.0
    dμ_ij_xy_sums = 0.0

    K = length(s_n)
    for k in 1:K
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[k]

        @inbounds @fastmath s_n[k] = u/μ_ij_sums_n
        @inbounds @fastmath dμ_ij_x_sums+= u*x
        @inbounds @fastmath dμ_ij_y_sums+= u*y
        @inbounds @fastmath dμ_ij_xy_sums+= u*x*y
    end

    return dμ_ij_x_sums, dμ_ij_y_sums, dμ_ij_xy_sums
end

function calc_derSums_xy!(n::Int64,s_n::Vector{T},
                    x::Float64,Y::Matrix{Float64},
                    μ_ij::Matrix{T},
                    μ_ij_sums_n::T) where T

    dμ_ij_x_sums = 0.0
    dμ_ij_y_sums = 0.0
    dμ_ij_xy_sums = 0.0

    K = length(s_n)
    for k in 1:K
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath y = Y[n,k]

        @inbounds @fastmath s_n[k] = u/μ_ij_sums_n
        @inbounds @fastmath dμ_ij_x_sums+= u*x
        @inbounds @fastmath dμ_ij_y_sums+= u*y
        @inbounds @fastmath dμ_ij_xy_sums+= u*x*y
    end

    return dμ_ij_x_sums, dμ_ij_y_sums, dμ_ij_xy_sums
end

function calc_derSums_xy!(n::Int64,s_n::Vector{T},
                    x::Float64,Y::Vector{Float64},
                    μ_ij::Matrix{T},
                    μ_ij_sums_n::T) where T

    dμ_ij_x_sums = 0.0
    dμ_ij_y_sums = 0.0
    dμ_ij_xy_sums = 0.0

    K = length(s_n)
    for k in 1:K
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath y = Y[k]

        @inbounds @fastmath s_n[k] = u/μ_ij_sums_n
        @inbounds @fastmath dμ_ij_x_sums+= u*x
        @inbounds @fastmath dμ_ij_y_sums+= u*y
        @inbounds @fastmath dμ_ij_xy_sums+= u*x*y
    end

    return dμ_ij_x_sums, dμ_ij_y_sums, dμ_ij_xy_sums
end

function calc_derSums_xy!(n::Int64,s_n::Vector{T},
                    x::Float64,y::Float64,
                    μ_ij::Matrix{T},
                    μ_ij_sums_n::T) where T

    dμ_ij_x_sums = 0.0
    dμ_ij_y_sums = 0.0
    dμ_ij_xy_sums = 0.0

    K = length(s_n)
    for k in 1:K
        @inbounds @fastmath u = μ_ij[n,k]

        @inbounds @fastmath s_n[k] = u/μ_ij_sums_n
        @inbounds @fastmath dμ_ij_x_sums+= u*x
        @inbounds @fastmath dμ_ij_y_sums+= u*y
        @inbounds @fastmath dμ_ij_xy_sums+= u*x*y
    end

    return dμ_ij_x_sums, dμ_ij_y_sums, dμ_ij_xy_sums
end


function calc_derSums_x!(n::Int64,
                    dμ_ij_x_sums::Vector{T},
                    X_mat::Matrix{Float64},
                    μ_ij::Matrix{T},
                    yr_next::Vector{Int}) where T

    dμ_ij_x_sums[:].= 0.0

    K = size(μ_ij,2)
    yr_ind = 1
    for k in 1:K
        if k in yr_next
            yr_ind +=1
        end
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[n,k]
        @inbounds @fastmath dμ_ij_x_sums[yr_ind] += u*x
    end

    return nothing
end

function calc_derSums_x!(n::Int64,
                    dμ_ij_x_sums::Vector{T},
                    X_mat::Vector{Float64},
                    μ_ij::Matrix{T},
                    yr_next::Vector{Int}) where T

    dμ_ij_x_sums[:].= 0.0

    K = size(μ_ij,2)
    yr_ind = 1
    ind_max = length(yr_next)
    for k in 1:K
        if (yr_ind<=ind_max)
            if (k == yr_next[yr_ind])
                yr_ind +=1
            end
        end
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[k]
        @inbounds @fastmath dμ_ij_x_sums[yr_ind] += u*x
    end

    return nothing
end

function calc_derSums_x!(n::Int64,
                    dμ_ij_x_sums::Vector{T},
                    x::Float64,
                    μ_ij::Matrix{T},
                    yr_next::Vector{Int}) where T

    dμ_ij_x_sums[:].= 0.0

    K = size(μ_ij,2)
    yr_ind = 1
    for k in 1:K
        if k in yr_next
            yr_ind +=1
        end
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath dμ_ij_x_sums[yr_ind] += u*x
    end

    return nothing
end



function calc_dll(dS_x::Matrix{T},n::Int64,X_mat::Matrix{Float64},Γ_x::Vector{T},ω::Vector{T},
                    s_n::Matrix{T},s_hat::Matrix{T},choice_ind::Vector{Int64}) where T
    dlnLdβ = 0.0
    yr_ind = 1
    for (y,k) in enumerate(choice_ind)
        s = (X_mat[n,k]-Γ_x[yr_ind])*s_n[k,n]
        dS_x[y,n] = s
        dlnLdβ += ω[k]*s./s_hat[k,n]
        yr_ind+=1
    end
    # dlnLdβ = sum(dlnSdβ[choice_ind])
    # ll_n = prod(s_hat[choice_ind,n])
    return dlnLdβ
end

function calc_dll(dS_x::Matrix{T},n::Int64,X_mat::Vector{Float64},Γ_x::Vector{T},ω::Vector{T},
                    s_n::Matrix{T},s_hat::Matrix{T},choice_ind::Vector{Int64}) where T
    dlnLdβ = 0.0
    yr_ind = 1
    for (y,k) in enumerate(choice_ind)
        s = (X_mat[k]-Γ_x[yr_ind])*s_n[k,n]
        dS_x[y,n] = s
        dlnLdβ += ω[k]*s./s_hat[k,n]
        yr_ind+=1
    end
    # dlnLdβ = sum(dlnSdβ[choice_ind])
    # ll_n = prod(s_hat[choice_ind,n])
    return dlnLdβ
end

function calc_dll(dS_x::Matrix{T},n::Int64,x::Float64,Γ_x::Vector{T},ω::Vector{T},
                    s_n::Matrix{T},s_hat::Matrix{T},choice_ind::Vector{Int64}) where T
    dlnLdβ = 0.0
    yr_ind = 1
    for (y,k) in enumerate(choice_ind)
        s = (x-Γ_x[yr_ind])*s_n[k,n]
        dS_x[y,n] = s
        dlnLdβ += ω[k]*s./s_hat[k,n]
        yr_ind+=1
    end
    # dlnLdβ = sum(dlnSdβ[choice_ind])
    # ll_n = prod(s_hat[choice_ind,n])
    return dlnLdβ
end

function calc_dll_ω(dS_x::Matrix{T},n::Int64,X_mat::Matrix{Float64},y_last::Vector{Float64},ω::Vector{T},
                    s_n::Matrix{T},s_hat::Matrix{T},choice_ind::Vector{Int64}) where T
    dlnLdβ = 0.0
    for (y,k) in enumerate(choice_ind)
        s = X_mat[n,k]*ω[k]*(1 - ω[k])
        dS_x[y,n] = s
        dlnLdβ += s*(s_n[k,n] - y_last[k])/s_hat[k,n]
    end
    # ll_n = prod(s_hat[choice_ind,n])
    return dlnLdβ
end

function calc_dll_ω(dS_x::Matrix{T},n::Int64,X_mat::Vector{Float64},y_last::Vector{Float64},ω::Vector{T},
                    s_n::Matrix{T},s_hat::Matrix{T},choice_ind::Vector{Int64}) where T
    dlnLdβ = 0.0
    for (y,k) in enumerate(choice_ind)
        s = X_mat[k]*ω[k]*(1 - ω[k])
        dS_x[y,n] = s
        dlnLdβ += s*(s_n[k,n] - y_last[k])/s_hat[k,n]
    end
    # ll_n = prod(s_hat[choice_ind,n])
    return dlnLdβ
end

function calc_dll_ω(dS_x::Matrix{T},n::Int64,x::Float64,y_last::Vector{Float64},ω::Vector{T},
                    s_n::Matrix{T},s_hat::Matrix{T},choice_ind::Vector{Int64}) where T
    dlnLdβ = 0.0
    for (y,k) in enumerate(choice_ind)
        s = x*ω[k]*(1 - ω[k])
        dS_x[y,n] = s
        dlnLdβ += s*(s_n[k,n] - y_last[k])/s_hat[k,n]
    end
    # ll_n = prod(s_hat[choice_ind,n])
    return dlnLdβ
end




























function calc_prodTerms_xy!(n::Int64,dS_xy::Vector{Float64},
                        X_mat::Matrix{Float64},Y::Matrix{Float64},
                        s_n::Vector{Float64},
                        Γ_x::Float64,Γ_y::Float64,Γ_xy::Float64)

    K = length(dS_xy)
    t_last = Γ_x*Γ_y -Γ_xy
    for k in 1:K
        @inbounds @fastmath txy = s_n[k]*((X_mat[n,k] - Γ_x)*(Y[n,k] - Γ_y)+t_last)
        @inbounds @fastmath dS_xy[k]+= txy
    end
end

function calc_prodTerms_xy!(n::Int64,dS_xy::Vector{Float64},
                        X_mat::Matrix{Float64},Y::Vector{Float64},
                        s_n::Vector{Float64},
                        Γ_x::Float64,Γ_y::Float64,Γ_xy::Float64)

    K = length(dS_xy)
    t_last = Γ_x*Γ_y -Γ_xy
    for k in 1:K
        @inbounds @fastmath txy = s_n[k]*((X_mat[n,k] - Γ_x)*(Y[k] - Γ_y)+t_last)
        @inbounds @fastmath dS_xy[k]+= txy
    end
end

function calc_prodTerms_xy!(n::Int64,dS_xy::Vector{Float64},
                        X_mat::Matrix{Float64},Y::Float64,
                        s_n::Vector{Float64},
                        Γ_x::Float64,Γ_y::Float64,Γ_xy::Float64)

    K = length(dS_xy)
    t_last = Γ_x*Γ_y -Γ_xy
    for k in 1:K
        @inbounds @fastmath txy = s_n[k]*((X_mat[n,k] - Γ_x)*(Y - Γ_y)+t_last)
        @inbounds @fastmath dS_xy[k]+= txy
    end
end



function calc_prodTerms_xy!(n::Int64,dS_xy::Vector{Float64},
                        X_mat::Vector{Float64},Y::Matrix{Float64},
                        s_n::Vector{Float64},
                        Γ_x::Float64,Γ_y::Float64,Γ_xy::Float64)

    K = length(dS_xy)
    t_last = Γ_x*Γ_y -Γ_xy
    for k in 1:K
        @inbounds @fastmath txy = s_n[k]*((X_mat[k] - Γ_x)*(Y[n,k] - Γ_y)+t_last)
        @inbounds @fastmath dS_xy[k]+= txy
    end
end

function calc_prodTerms_xy!(n::Int64,dS_xy::Vector{Float64},
                        X_mat::Vector{Float64},Y::Vector{Float64},
                        s_n::Vector{Float64},
                        Γ_x::Float64,Γ_y::Float64,Γ_xy::Float64)

    K = length(dS_xy)
    t_last = Γ_x*Γ_y -Γ_xy
    for k in 1:K
        @inbounds @fastmath txy = s_n[k]*((X_mat[k] - Γ_x)*(Y[k] - Γ_y)+t_last)
        @inbounds @fastmath dS_xy[k]+= txy
    end
end

function calc_prodTerms_xy!(n::Int64,dS_xy::Vector{Float64},
                        X_mat::Vector{Float64},Y::Float64,
                        s_n::Vector{Float64},
                        Γ_x::Float64,Γ_y::Float64,Γ_xy::Float64)

    K = length(dS_xy)
    t_last = Γ_x*Γ_y -Γ_xy
    for k in 1:K
        @inbounds @fastmath txy = s_n[k]*((X_mat[k] - Γ_x)*(Y - Γ_y)+t_last)
        @inbounds @fastmath dS_xy[k]+= txy
    end
end


function calc_prodTerms_xy!(n::Int64,dS_xy::Vector{Float64},
                        X_mat::Float64,Y::Matrix{Float64},
                        s_n::Vector{Float64},
                        Γ_x::Float64,Γ_y::Float64,Γ_xy::Float64)

    K = length(dS_xy)
    t_last = Γ_x*Γ_y -Γ_xy
    for k in 1:K
        @inbounds @fastmath txy = s_n[k]*((X_mat - Γ_x)*(Y[n,k] - Γ_y)+t_last)
        @inbounds @fastmath dS_xy[k]+= txy
    end
end

function calc_prodTerms_xy!(n::Int64,dS_xy::Vector{Float64},
                        X_mat::Float64,Y::Vector{Float64},
                        s_n::Vector{Float64},
                        Γ_x::Float64,Γ_y::Float64,Γ_xy::Float64)

    K = length(dS_xy)
    t_last = Γ_x*Γ_y -Γ_xy
    for k in 1:K
        @inbounds @fastmath txy = s_n[k]*((X_mat - Γ_x)*(Y[k] - Γ_y)+t_last)
        @inbounds @fastmath dS_xy[k]+= txy
    end
end

function calc_prodTerms_xy!(n::Int64,dS_xy::Vector{Float64},
                        X_mat::Float64,Y::Float64,
                        s_n::Vector{Float64},
                        Γ_x::Float64,Γ_y::Float64,Γ_xy::Float64)

    K = length(dS_xy)
    t_last = Γ_x*Γ_y -Γ_xy
    for k in 1:K
        @inbounds @fastmath txy = s_n[k]*((X_mat - Γ_x)*(Y - Γ_y)+t_last)
        @inbounds @fastmath dS_xy[k]+= txy
    end
end

function calc_prodTerms_x!(n::Int64,
                        dS_x::Vector{T},
                        X_mat::Matrix{Float64},
                        s_n::Vector{T},
                        Γ_x::T) where T

    K = length(dS_x)
    for k in 1:K
        @inbounds @fastmath tx = s_n[k]*(X_mat[n,k] - Γ_x)
        @inbounds @fastmath dS_x[k]+= tx
    end
end

function calc_prodTerms_x!(n::Int64,
                        dS_x::Vector{T},
                        X_mat::Vector{Float64},
                        s_n::Vector{T},
                        Γ_x::T) where T

    K = length(dS_x)
    for k in 1:K
        @inbounds @fastmath tx = s_n[k]*(X_mat[k] - Γ_x)
        @inbounds @fastmath dS_x[k]+= tx
    end
end

function calc_prodTerms_x!(n::Int64,
                        dS_x::Vector{T},
                        X_mat::Float64,
                        s_n::Vector{T},
                        Γ_x::T) where T

    K = length(dS_x)
    for k in 1:K
        @inbounds @fastmath tx = s_n[k]*(X_mat - Γ_x)
        @inbounds @fastmath dS_x[k]+= tx
    end
end



function hess_calc_ω!(dS_xy::Vector{T},s_uncond::Vector{T},
                    ω::T,
                    X_mat::Matrix{Float64},Y_mat::Matrix{Float64},
                    y_last::Vector{Float64}) where T

    dS_xy[:] .= 0.0
    K = length(dS_xy)
    for k in 1:K
        dS_xy[k] = (X_mat[1,k]*Y_mat[1,k]*ω*(1-ω)*(1-2*ω))*(s_uncond[k] - y_last[k])
    end
    return nothing
end

function hess_calc_ω!(dS_xy::Vector{T},s_uncond::Vector{T},
                    ω::T,
                    X_mat::Matrix{Float64},Y_mat::Vector{Float64},
                    y_last::Vector{Float64}) where T

    dS_xy[:] .= 0.0
    K = length(dS_xy)
    for k in 1:K
        dS_xy[k] = (X_mat[1,k]*Y_mat[k]*ω*(1-ω)*(1-2*ω))*(s_uncond[k] - y_last[k])
    end
    return nothing
end

function hess_calc_ω!(dS_xy::Vector{T},s_uncond::Vector{T},
                    ω::T,
                    X_mat::Matrix{Float64},Y_mat::Float64,
                    y_last::Vector{Float64}) where T

    dS_xy[:] .= 0.0
    K = length(dS_xy)
    for k in 1:K
        dS_xy[k] = (X_mat[1,k]*Y_mat*ω*(1-ω)*(1-2*ω))*(s_uncond[k] - y_last[k])
    end
    return nothing
end


function hess_calc_ω!(dS_xy::Vector{T},s_uncond::Vector{T},
                    ω::T,
                    X_mat::Vector{Float64},Y_mat::Matrix{Float64},
                    y_last::Vector{Float64}) where T

    dS_xy[:] .= 0.0
    K = length(dS_xy)
    for k in 1:K
        dS_xy[k] = (X_mat[k]*Y_mat[1,k]*ω*(1-ω)*(1-2*ω))*(s_uncond[k] - y_last[k])
    end
    return nothing
end

function hess_calc_ω!(dS_xy::Vector{T},s_uncond::Vector{T},
                    ω::T,
                    X_mat::Vector{Float64},Y_mat::Vector{Float64},
                    y_last::Vector{Float64}) where T

    dS_xy[:] .= 0.0
    K = length(dS_xy)
    for k in 1:K
        dS_xy[k] = (X_mat[k]*Y_mat[k]*ω*(1-ω)*(1-2*ω))*(s_uncond[k] - y_last[k])
    end
    return nothing
end

function hess_calc_ω!(dS_xy::Vector{T},s_uncond::Vector{T},
                    ω::T,
                    X_mat::Vector{Float64},Y_mat::Float64,
                    y_last::Vector{Float64}) where T

    dS_xy[:] .= 0.0
    K = length(dS_xy)
    for k in 1:K
        dS_xy[k] = (X_mat[k]*Y_mat*ω*(1-ω)*(1-2*ω))*(s_uncond[k] - y_last[k])
    end
    return nothing
end

function hess_calc_ω!(dS_xy::Vector{T},s_uncond::Vector{T},
                    ω::T,
                    X_mat::Float64,Y_mat::Matrix{Float64},
                    y_last::Vector{Float64}) where T

    dS_xy[:] .= 0.0
    K = length(dS_xy)
    for k in 1:K
        dS_xy[k] = (X_mat*Y_mat[1,k]*ω*(1-ω)*(1-2*ω))*(s_uncond[k] - y_last[k])
    end
    return nothing
end

function hess_calc_ω!(dS_xy::Vector{T},s_uncond::Vector{T},
                    ω::T,
                    X_mat::Float64,Y_mat::Vector{Float64},
                    y_last::Vector{Float64}) where T

    dS_xy[:] .= 0.0
    K = length(dS_xy)
    for k in 1:K
        dS_xy[k] = (X_mat*Y_mat[k]*ω*(1-ω)*(1-2*ω))*(s_uncond[k] - y_last[k])
    end
    return nothing
end

function hess_calc_ω!(dS_xy::Vector{T},s_uncond::Vector{T},
                    ω::T,
                    X_mat::Float64,Y_mat::Float64,
                    y_last::Vector{Float64}) where T

    dS_xy[:] .= 0.0
    K = length(dS_xy)
    for k in 1:K
        dS_xy[k] = (X_mat*Y_mat*ω*(1-ω)*(1-2*ω))*(s_uncond[k] - y_last[k])
    end
    return nothing
end

function grad_calc_ω!(dω_x::Vector{T},
                    ω::T,
                    X_mat::Matrix{Float64}) where T

    dω_x[:] .= 0.0
    K = length(dω_x)
    for k in 1:K
        dω_x[k] = (X_mat[1,k]*ω*(1-ω))#*(s_uncond[k] - y_last[k])
    end
    return nothing
end

function grad_calc_ω!(dω_x::Vector{T},
                    ω::T,
                    X_mat::Vector{Float64}) where T

    dω_x[:] .= 0.0
    K = length(dω_x)
    for k in 1:K
        dω_x[k] = (X_mat[k]*ω*(1-ω))#*(s_uncond[k] - y_last[k])
    end
    return nothing
end

function grad_calc_ω!(dω_x::Vector{T},
                    ω::T,
                    X_mat::Float64) where T

    dω_x[:] .= 0.0
    K = length(dω_x)
    val = (X_mat*ω*(1-ω))
    for k in 1:K
        dω_x[k] = val#*(s_uncond[k] - y_last[k])
    end
    return nothing
end


function getIndex_broad(A::Float64,n::Int64,k::Int64)
    return A
end

function getIndex_broad(A::Array{Float64,1},n::Int64,k::Int64)
    return A[k]
end

function getIndex_broad(A::Array{Float64,2},n::Int64,k::Int64)
    return A[n,k]
end
