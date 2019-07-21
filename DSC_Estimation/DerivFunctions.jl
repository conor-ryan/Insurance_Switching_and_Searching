function calc_derSums_xy!(n::Int64,
                    dμ_ij_x_sums::Vector{T},
                    dμ_ij_y_sums::Vector{T},
                    dμ_ij_xy_sums::Vector{T},
                    X_mat::Matrix{Float64},
                    Y_mat::Matrix{Float64},
                    μ_ij::Matrix{T},
                    yr_next::Vector{Int}) where T

    dμ_ij_x_sums[:].= 0.0
    dμ_ij_y_sums[:].= 0.0
    dμ_ij_xy_sums[:].= 0.0

    K = size(μ_ij,2)
    yr_ind = 1
    for k in 1:K
        if k in yr_next
            yr_ind +=1
        end
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[n,k]
        @inbounds @fastmath y = Y_mat[n,k]

        @inbounds @fastmath dμ_ij_x_sums[yr_ind] += u*x
        @inbounds @fastmath dμ_ij_y_sums[yr_ind] += u*y
        @inbounds @fastmath dμ_ij_xy_sums[yr_ind] += u*x*y
    end

    return nothing
end

function calc_derSums_xy!(n::Int64,
                    dμ_ij_x_sums::Vector{T},
                    dμ_ij_y_sums::Vector{T},
                    dμ_ij_xy_sums::Vector{T},
                    X_mat::Matrix{Float64},
                    Y_mat::Vector{Float64},
                    μ_ij::Matrix{T},
                    yr_next::Vector{Int}) where T

    dμ_ij_x_sums[:].= 0.0
    dμ_ij_y_sums[:].= 0.0
    dμ_ij_xy_sums[:].= 0.0

    K = size(μ_ij,2)
    yr_ind = 1
    for k in 1:K
        if k in yr_next
            yr_ind +=1
        end
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[n,k]
        @inbounds @fastmath y = Y_mat[k]

        @inbounds @fastmath dμ_ij_x_sums[yr_ind] += u*x
        @inbounds @fastmath dμ_ij_y_sums[yr_ind] += u*y
        @inbounds @fastmath dμ_ij_xy_sums[yr_ind] += u*x*y
    end

    return nothing
end

function calc_derSums_xy!(n::Int64,
                    dμ_ij_x_sums::Vector{T},
                    dμ_ij_y_sums::Vector{T},
                    dμ_ij_xy_sums::Vector{T},
                    X_mat::Vector{Float64},
                    Y_mat::Matrix{Float64},
                    μ_ij::Matrix{T},
                    yr_next::Vector{Int}) where T

    dμ_ij_x_sums[:].= 0.0
    dμ_ij_y_sums[:].= 0.0
    dμ_ij_xy_sums[:].= 0.0

    K = size(μ_ij,2)
    yr_ind = 1
    for k in 1:K
        if k in yr_next
            yr_ind +=1
        end
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[k]
        @inbounds @fastmath y = Y_mat[n,k]

        @inbounds @fastmath dμ_ij_x_sums[yr_ind] += u*x
        @inbounds @fastmath dμ_ij_y_sums[yr_ind] += u*y
        @inbounds @fastmath dμ_ij_xy_sums[yr_ind] += u*x*y
    end

    return nothing
end

function calc_derSums_xy!(n::Int64,
                    dμ_ij_x_sums::Vector{T},
                    dμ_ij_y_sums::Vector{T},
                    dμ_ij_xy_sums::Vector{T},
                    X_mat::Vector{Float64},
                    Y_mat::Vector{Float64},
                    μ_ij::Matrix{T},
                    yr_next::Vector{Int}) where T

    dμ_ij_x_sums[:].= 0.0
    dμ_ij_y_sums[:].= 0.0
    dμ_ij_xy_sums[:].= 0.0

    K = size(μ_ij,2)
    yr_ind = 1
    for k in 1:K
        if k in yr_next
            yr_ind +=1
        end
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[k]
        @inbounds @fastmath y = Y_mat[k]

        @inbounds @fastmath dμ_ij_x_sums[yr_ind] += u*x
        @inbounds @fastmath dμ_ij_y_sums[yr_ind] += u*y
        @inbounds @fastmath dμ_ij_xy_sums[yr_ind] += u*x*y
    end

    return nothing
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

function calc_d2S(dS_xy::Matrix{T},n::Int64,
                    X_mat::Matrix{Float64},Y_mat::Matrix{Float64},
                    Γ_x::Vector{T},Γ_y::Vector{T},Γ_xy::Vector{T},
                    ω::Vector{T},s_n::Matrix{T},choice_ind::Vector{Int64}) where T
    for (c,k) in enumerate(choice_ind)
        s = ((X_mat[n,k]-Γ_x[c])*(Y_mat[n,k]-Γ_y[c])-Γ_xy[c] + Γ_x[c]*Γ_y[c])*s_n[k,n]
        dS_xy[c,n] = ω[k]*s
    end
    return nothing
end

function calc_d2S(dS_xy::Matrix{T},n::Int64,
                    X_mat::Matrix{Float64},Y_mat::Vector{Float64},
                    Γ_x::Vector{T},Γ_y::Vector{T},Γ_xy::Vector{T},
                    ω::Vector{T},s_n::Matrix{T},choice_ind::Vector{Int64}) where T
    for (c,k) in enumerate(choice_ind)
        s = ((X_mat[n,k]-Γ_x[c])*(Y_mat[k]-Γ_y[c])-Γ_xy[c] + Γ_x[c]*Γ_y[c])*s_n[k,n]
        dS_xy[c,n] = ω[k]*s
    end
    return nothing
end

function calc_d2S(dS_xy::Matrix{T},n::Int64,
                    X_mat::Vector{Float64},Y_mat::Matrix{Float64},
                    Γ_x::Vector{T},Γ_y::Vector{T},Γ_xy::Vector{T},
                    ω::Vector{T},s_n::Matrix{T},choice_ind::Vector{Int64}) where T
    for (c,k) in enumerate(choice_ind)
        s = ((X_mat[k]-Γ_x[c])*(Y_mat[n,k]-Γ_y[c])-Γ_xy[c] + Γ_x[c]*Γ_y[c])*s_n[k,n]
        dS_xy[c,n] = ω[k]*s
    end
    return nothing
end

function calc_d2S(dS_xy::Matrix{T},n::Int64,
                    X_mat::Vector{Float64},Y_mat::Vector{Float64},
                    Γ_x::Vector{T},Γ_y::Vector{T},Γ_xy::Vector{T},
                    ω::Vector{T},s_n::Matrix{T},choice_ind::Vector{Int64}) where T
    for (c,k) in enumerate(choice_ind)
        s = ((X_mat[k]-Γ_x[c])*(Y_mat[k]-Γ_y[c])-Γ_xy[c] + Γ_x[c]*Γ_y[c])*s_n[k,n]
        dS_xy[c,n] = ω[k]*s
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

function calc_d2S_ω(dS_xy::Matrix{T},n::Int64,X_mat::Vector{Float64},Y_mat::Vector{Float64},
                    y_last::Vector{Float64},ω::Vector{T},
                    s_n::Matrix{T},choice_ind::Vector{Int64}) where T
    for (c,k) in enumerate(choice_ind)
        s = X_mat[k]*Y_mat[k]*ω[k]*(1 - ω[k])*(1-2*ω[k])*(s_n[k,n] - y_last[k])
        dS_xy[c,n] = s
    end
    return nothing
end


function calc_d2ll(n::Int64,dS_xy::Matrix{T},dS_x::Matrix{T},dS_y::Matrix{T},s_hat::Matrix{T},choice_ind::Vector{Int64}) where T
    dlnLdβ = 0.0
    for (c,k) in enumerate(choice_ind)
        dlnLdβ += dS_xy[c,n]/s_hat[k,n] - (dS_x[c,n]*dS_y[c,n])/(s_hat[k,n])^2
    end
    return dlnLdβ
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
