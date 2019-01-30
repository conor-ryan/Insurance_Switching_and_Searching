function unPackChars(app::ChoiceData,d::InsuranceLogit)
    ind = person(app)[1]
    S_ij = transpose(choice(app))
    wgt = transpose(weight(app))
    idxitr = d.data._personDict[ind]

    X_t = prodchars(app)
    X_0_t = prodchars0(app)
    Z = demoRaw(app)[:,1]
    #F_t = fixedEffects(app)
    F_t = fixedEffects(app,idxitr)

    return ind, S_ij, wgt,idxitr, X_t, X_0_t, Z, F_t
end

function unPackParChars(p::parDict{T},idxitr::UnitRange{Int}) where T
    μ_ij = p.μ_ij[:,idxitr]
    #Get Market Shares
    s_hat = p.s_hat[idxitr]

    return μ_ij,s_hat
end

function sumShares!(s_hat::Vector{T}) where T
    # Fix possible computational error
    for k in eachindex(s_hat)
        if abs(s_hat[k])<=1e-300
            s_hat[k]=1e-15
            #println("Hit Share Constraint for person $ind, product $k")
        end
    end
    # s_insured = sum(s_hat)
    # if s_insured>=(1-1e-300)
    #     s_insured= 1 - 1e-15
    #     #println("Hit insured constraint for person $ind")
    # end

    return nothing
end

function relPar(app::ChoiceData,d::InsuranceLogit,F_t::SubArray,ind::Float64)
    Q = d.parLength[:All]
    Q_0 = Q - size(F_t,1)
    ## Relevant Parameters for this observation
    pars_relevant = vcat(1:Q_0,Q_0 .+ app._rel_fe_Dict[ind])

    return pars_relevant
end

function preCalcμ(μ_ij::Matrix{T}) where T
    μ_ij_sums = sum(μ_ij,dims=2)
    return μ_ij_sums[:]
end

function ll_Terms(wgt::Array{Float64,N},S_ij::Array{Float64,N},
                    s_hat::Vector{T}) where {N,T}
    K = length(S_ij)
    ll_obs = 0.0
    gll_t1 = Vector{T}(undef,K)
    gll_t2 = Vector{T}(undef,K)
    gll_t3 = Vector{T}(undef,K)
    # gll_t4 = Vector{T}(undef,K)
    # gll_t5 = Vector{T}(undef,K)
    # gll_t6 = Vector{T}(undef,K)
    for k in 1:K
        #Gradient Terms
        gll_t1[k] = wgt[k]*S_ij[k]*(1/s_hat[k])
        gll_t2[k] = wgt[k]*S_ij[k]*(1/s_hat[k]^2)
        gll_t3[k] = wgt[k]*S_ij[k]*2*(1/s_hat[k]^3)
        # gll_t4[k] = wgt[k]*S_ij[k]*urate[k]*(1/(s_insured) + 1/(1-s_insured))
        # gll_t5[k] = wgt[k]*S_ij[k]*urate[k]*(1/(s_insured^2) - 1/((1-s_insured)^2))
        # gll_t6[k] = wgt[k]*S_ij[k]*urate[k]*2*(1/(s_insured^3) + 1/((1-s_insured)^3))
        # Log Likelihood
        # ll_obs+=wgt[k]*S_ij[k]*(log(s_hat[k]) -
        #                 urate[k]*(log(s_insured)-log(1-s_insured)))
        ll_obs+=wgt[k]*S_ij[k]*log(s_hat[k])
    end

    return gll_t1, gll_t2, gll_t3, ll_obs
end


function ll_obs!(thD::Array{Float64,3},hess::Matrix{Float64},grad::Vector{Float64},
                            app::ChoiceData,d::InsuranceLogit,p::parDict{T}) where T

        ind, S_ij, wgt, idxitr, X_t, X_0_t, Z, F_t = unPackChars(app,d)
        wgt = convert(Array{Float64,2},wgt)
        S_ij = convert(Array{Float64,2},S_ij)

        prodidx = Int.(product(app))

        draws = d.draws

        # Get Utility and derivative of Utility
        μ_ij,s_hat = unPackParChars(p,idxitr)

        sumShares!(s_hat)

        (N,K) = size(μ_ij)

        # Initialize Gradient
        #(Q,N,K) = size(dμ_ij)
        pars_relevant = relPar(app,d,F_t,ind)

        # Pre-Calculate Squares
        μ_ij_sums = preCalcμ(μ_ij)

        # Pre-Calculate Log-Likelihood Terms for Gradient
        # Also Calculate Log-Likelihood itself
        gll_t1, gll_t2, gll_t3, ll_obs = ll_Terms(wgt,S_ij,s_hat)

        #hess = zeros(Q,Q)
        #hess[:] = 0.0
        #grad[:] = 0.0
        X_mat = Array{Float64}(undef,N,K)
        Y_list = Array{Array{Float64,2},1}(undef,length(pars_relevant))
        #Y_mat = Array{Float64}(N,K)

        # Allocate Memory
        dS_xyz = Vector{Float64}(undef,K)
        dS_xy = Vector{Float64}(undef,K)
        dS_x = Vector{Float64}(undef,K)


        dS_x_list = Vector{Vector{Float64}}(undef,length(pars_relevant))
        dS_xy_list = Array{Vector{Float64},2}(undef,length(pars_relevant),length(pars_relevant))

        s_n = Vector{Float64}(undef,K)


        #γlen = 1 + d.parLength[:γ]
        γlen = 0
        β0len = γlen + d.parLength[:β]
        βlen = β0len + d.parLength[:γ]
        σlen = βlen + d.parLength[:σ]
        FElen = σlen + d.parLength[:FE]


        for (q_i,q) in enumerate(pars_relevant)
            returnParameter!(q,X_mat,
                            Z,X_0_t,X_t,draws,F_t,
                            γlen,β0len,βlen,σlen)
            @inbounds Y_list[q_i] = X_mat[:,:]

            grad_calc!(dS_x,s_n,
                        μ_ij,X_mat,
                        μ_ij_sums)

            # Save Answers
            dS_x_list[q_i] = dS_x[:]

            ## Calculate Gradient
            grad[q]+= combine_grad(N,gll_t1,dS_x)


            for (r_i,r) in enumerate(pars_relevant)
                if r>q
                    continue
                end
                Y_mat = Y_list[r_i]


                dS_y = dS_x_list[r_i]

                hess_calc!(dS_xy,s_n,
                            μ_ij,X_mat,Y_mat,
                            μ_ij_sums)

                dS_xy_list[q_i,r_i] = dS_xy[:]

                hess_obs = combine_hess(N,gll_t1,gll_t2,
                            dS_xy,dS_x,dS_y)

                hess[q,r]+= hess_obs
                if (q!=r)
                    hess[r,q]+= hess_obs
                end


                for (t_i,t) in enumerate(pars_relevant)
                    if t>r
                        continue
                    end
                    Z_mat = Y_list[t_i]

                    dS_z = dS_x_list[t_i]

                    dS_xz = dS_xy_list[q_i,t_i]

                    dS_yz = dS_xy_list[r_i,t_i]


                    thD_calc!(dS_xyz,s_n,
                                μ_ij,X_mat,Y_mat,Z_mat,
                                μ_ij_sums)

                    thD_obs = combine_thD(N,gll_t1,gll_t2,gll_t3,
                                dS_xyz,dS_xy,dS_xz,dS_yz,
                                dS_x,dS_y,dS_z)

                    thD[q,r,t]+= thD_obs
                    if (r!=t)
                        thD[q,t,r]+= thD_obs
                        thD[t,r,q]+= thD_obs
                        if (r!=q)
                            thD[r,t,q]+= thD_obs
                            thD[t,q,r]+= thD_obs
                            thD[r,q,t]+= thD_obs
                        end
                    elseif (r!=q)
                        thD[t,r,q]+= thD_obs
                        thD[r,q,t]+= thD_obs
                    end

                end

            end
        end

    return ll_obs,pars_relevant
end


function ll_obs!(hess::Matrix{Float64},grad::Vector{Float64},
                    app::ChoiceData,d::InsuranceLogit,p::parDict{T}) where T

        ind, S_ij, wgt, idxitr, X_t, X_0_t, Z, F_t = unPackChars(app,d)
        wgt = convert(Array{Float64,2},wgt)
        S_ij = convert(Array{Float64,2},S_ij)

        prodidx = Int.(product(app))

        draws = d.draws

        # Get Utility and derivative of Utility
        μ_ij,s_hat = unPackParChars(p,idxitr)

        sumShares!(s_hat)

        (N,K) = size(μ_ij)

        # Initialize Gradient
        #(Q,N,K) = size(dμ_ij)
        pars_relevant = relPar(app,d,F_t,ind)

        # Pre-Calculate Squares
        μ_ij_sums = preCalcμ(μ_ij)

        # Pre-Calculate Log-Likelihood Terms for Gradient
        # Also Calculate Log-Likelihood itself
        gll_t1, gll_t2, gll_t3, ll_obs = ll_Terms(wgt,S_ij,s_hat)

        #hess = zeros(Q,Q)
        #hess[:] = 0.0
        #grad[:] = 0.0
        X_mat = Array{Float64}(undef,N,K)
        Y_list = Array{Array{Float64,2},1}(undef,length(pars_relevant))
        #Y_mat = Array{Float64}(N,K)

        # Allocate Memory
        # dS_xyz = Vector{Float64}(undef,K)
        dS_xy = Vector{Float64}(undef,K)
        dS_x = Vector{Float64}(undef,K)


        dS_x_list = Vector{Vector{Float64}}(undef,length(pars_relevant))
        # dS_xy_list = Array{Vector{Float64},2}(undef,length(pars_relevant),length(pars_relevant))

        s_n = Vector{Float64}(undef,K)


        #γlen = 1 + d.parLength[:γ]
        γlen = 0
        β0len = γlen + d.parLength[:β]
        βlen = β0len + d.parLength[:γ]
        σlen = βlen + d.parLength[:σ]
        FElen = σlen + d.parLength[:FE]


        for (q_i,q) in enumerate(pars_relevant)
            returnParameter!(q,X_mat,
                            Z,X_0_t,X_t,draws,F_t,
                            γlen,β0len,βlen,σlen)
            @inbounds Y_list[q_i] = X_mat[:,:]

            grad_calc!(dS_x,s_n,
                        μ_ij,X_mat,
                        μ_ij_sums)

            # Save Answers
            dS_x_list[q_i] = dS_x[:]

            ## Calculate Gradient
            grad[q]+= combine_grad(N,gll_t1,dS_x)


            for (r_i,r) in enumerate(pars_relevant)
                if r>q
                    continue
                end
                Y_mat = Y_list[r_i]


                dS_y = dS_x_list[r_i]

                hess_calc!(dS_xy,s_n,
                            μ_ij,X_mat,Y_mat,
                            μ_ij_sums)

                # dS_xy_list[q_i,r_i] = dS_xy[:]

                hess_obs = combine_hess(N,gll_t1,gll_t2,
                            dS_xy,dS_x,dS_y)

                hess[q,r]+= hess_obs
                if (q!=r)
                    hess[r,q]+= hess_obs
                end
            end
        end

    return ll_obs,pars_relevant
end



function ll_obs!(grad::Vector{Float64},
                    app::ChoiceData,d::InsuranceLogit,p::parDict{T}) where T

        ind, S_ij, wgt, idxitr, X_t, X_0_t, Z, F_t = unPackChars(app,d)
        wgt = convert(Array{Float64,2},wgt)
        S_ij = convert(Array{Float64,2},S_ij)

        prodidx = Int.(product(app))

        draws = d.draws

        # Get Utility and derivative of Utility
        μ_ij,s_hat = unPackParChars(p,idxitr)

        sumShares!(s_hat)

        (N,K) = size(μ_ij)

        # Initialize Gradient
        #(Q,N,K) = size(dμ_ij)
        pars_relevant = relPar(app,d,F_t,ind)

        # Pre-Calculate Squares
        μ_ij_sums = preCalcμ(μ_ij)

        # Pre-Calculate Log-Likelihood Terms for Gradient
        # Also Calculate Log-Likelihood itself
        gll_t1, gll_t2, gll_t3, ll_obs = ll_Terms(wgt,S_ij,s_hat)

        #hess = zeros(Q,Q)
        #hess[:] = 0.0
        #grad[:] = 0.0
        X_mat = Array{Float64}(undef,N,K)
        # Y_list = Array{Array{Float64,2},1}(undef,length(pars_relevant))
        #Y_mat = Array{Float64}(N,K)

        # Allocate Memory
        # dS_xyz = Vector{Float64}(undef,K)
        # dS_xy = Vector{Float64}(undef,K)
        dS_x = Vector{Float64}(undef,K)


        # dS_x_list = Vector{Vector{Float64}}(undef,length(pars_relevant))
        # dS_xy_list = Array{Vector{Float64},2}(undef,length(pars_relevant),length(pars_relevant))

        s_n = Vector{Float64}(undef,K)


        #γlen = 1 + d.parLength[:γ]
        γlen = 0
        β0len = γlen + d.parLength[:β]
        βlen = β0len + d.parLength[:γ]
        σlen = βlen + d.parLength[:σ]
        FElen = σlen + d.parLength[:FE]


        for (q_i,q) in enumerate(pars_relevant)
            returnParameter!(q,X_mat,
                            Z,X_0_t,X_t,draws,F_t,
                            γlen,β0len,βlen,σlen)
            # @inbounds Y_list[q_i] = X_mat[:,:]

            grad_calc!(dS_x,s_n,
                        μ_ij,X_mat,
                        μ_ij_sums)

            # Save Answers
            # dS_x_list[q_i] = dS_x[:]

            ## Calculate Gradient
            grad[q]+= combine_grad(N,gll_t1,dS_x)

        end

    return ll_obs,pars_relevant
end


function combine_thD(N::Int64,
                    gll_t1::Vector{T},gll_t2::Vector{T},
                    gll_t3::Vector{T},
                    dS_xyz::Vector{T},
                    dS_xy::Vector{T},dS_xz::Vector{T},dS_yz::Vector{T},
                    dS_x::Vector{T},dS_y::Vector{T},dS_z::Vector{T}) where T
    thD_obs = 0.0
    K = length(dS_xyz)

    @inbounds @fastmath @simd for k in 1:K
        thD_obs += gll_t1[k]*dS_xyz[k]/N -
                        gll_t2[k]*((dS_x[k]/N)*(dS_yz[k]/N) + (dS_y[k]/N)*(dS_xz[k]/N) + (dS_z[k]/N)*(dS_xy[k]/N)) +
                        gll_t3[k]*(dS_x[k]/N)*(dS_y[k]/N)*(dS_z[k]/N)
    end
    return thD_obs
end


function combine_hess(N::Int64,
                    gll_t1::Vector{T},gll_t2::Vector{T},
                    dS_xy::Vector{T},dS_x::Vector{T},
                    dS_y::Vector{T}) where T
    hess_obs = 0.0
    K = length(dS_xy)

    @inbounds @fastmath @simd for k in 1:K
        hess_obs += gll_t1[k]*(dS_xy[k]/N) - gll_t2[k]*(dS_x[k]/N)*(dS_y[k]/N)
    end
    return hess_obs
end


function combine_grad(N::Int64,
                    gll_t1::Vector{T},
                    dS_x::Vector{T}) where T
    grad_obs = 0.0
    K = length(dS_x)

    @inbounds @fastmath @simd for k in 1:K
        grad_obs += gll_t1[k]*dS_x[k]
    end
    return grad_obs/N
end
function calc_derSums_xyz!(n::Int64,s_n::Vector{T},
                    X_mat::Matrix{Float64},Y::Matrix{Float64},
                    Z::Matrix{Float64},
                    μ_ij::Matrix{T},
                    μ_ij_sums_n::T) where T

    dμ_ij_x_sums = 0.0
    dμ_ij_y_sums = 0.0
    dμ_ij_z_sums = 0.0

    dμ_ij_xy_sums = 0.0
    dμ_ij_xz_sums = 0.0
    dμ_ij_yz_sums = 0.0

    dμ_ij_xyz_sums = 0.0

    K = length(s_n)
    for k in 1:K
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[n,k]
        @inbounds @fastmath y = Y[n,k]
        @inbounds @fastmath z = Z[n,k]

        @inbounds @fastmath s_n[k] = u/μ_ij_sums_n
        @inbounds @fastmath dμ_ij_x_sums+= u*x
        @inbounds @fastmath dμ_ij_y_sums+= u*y
        @inbounds @fastmath dμ_ij_z_sums+= u*z
        @inbounds @fastmath dμ_ij_xy_sums+= u*x*y
        @inbounds @fastmath dμ_ij_xz_sums+= u*x*z
        @inbounds @fastmath dμ_ij_yz_sums+= u*y*z
        @inbounds @fastmath dμ_ij_xyz_sums+= u*x*y*z
    end

    return dμ_ij_x_sums, dμ_ij_y_sums, dμ_ij_z_sums,dμ_ij_xy_sums,dμ_ij_xz_sums,dμ_ij_yz_sums,dμ_ij_xyz_sums
end

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

function calc_derSums_x!(n::Int64,s_n::Vector{T},
                    X_mat::Matrix{Float64},
                    μ_ij::Matrix{T},
                    μ_ij_sums_n::T) where T

    dμ_ij_x_sums = 0.0

    K = length(s_n)
    for k in 1:K
        @inbounds @fastmath u = μ_ij[n,k]
        @inbounds @fastmath x = X_mat[n,k]

        @inbounds @fastmath s_n[k] = u/μ_ij_sums_n
        @inbounds @fastmath dμ_ij_x_sums+= u*x
    end

    return dμ_ij_x_sums
end

function calc_prodTerms_xyz!(n::Int64,dS_xyz::Vector{Float64},
                        X_mat::Matrix{Float64},Y::Matrix{Float64},
                        Z::Matrix{Float64},s_n::Vector{Float64},
                        Γ_x::Float64,Γ_y::Float64,Γ_z::Float64,
                        Γ_xy::Float64,Γ_xz::Float64,Γ_yz::Float64,
                        Γ_xyz::Float64)

    K = length(dS_xyz)
    t_last = Γ_x*Γ_yz + Γ_y*Γ_xz + Γ_z*Γ_xy - Γ_xyz - 2*Γ_x*Γ_y*Γ_z
    t_1 = (Γ_x*Γ_y - Γ_xy)
    t_2 = (Γ_x*Γ_z - Γ_xz)
    t_3 = (Γ_y*Γ_z - Γ_yz)
    for k in 1:K
        @inbounds @fastmath dS_xyz[k]+= s_n[k]*(
        (X_mat[n,k] - Γ_x)*(Y[n,k] - Γ_y)*(Z[n,k] - Γ_z) +
        t_1*(Z[n,k] - Γ_z) +
        t_2*(Y[n,k] - Γ_y) +
        t_3*(X_mat[n,k] - Γ_x) + t_last )
    end
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


function hess_calc!(dS_xy::Vector{Float64},s_n::Vector{Float64},
                    μ_ij::Matrix{Float64},
                    X_mat::Matrix{Float64},Y::Matrix{Float64},
                    μ_ij_sums::Vector{Float64})
    dS_xy[:] .= 0.0


    (N,K) = size(μ_ij)

    ### 0 Draw Calculation ###
    for n in 1:N
        μ_ij_sums_n = μ_ij_sums[n]

        dμ_ij_x_sums, dμ_ij_y_sums, dμ_ij_xy_sums = calc_derSums_xy!(n,s_n,X_mat,Y,μ_ij,μ_ij_sums_n)

        @fastmath Γ_x = dμ_ij_x_sums/μ_ij_sums_n
        @fastmath Γ_y = dμ_ij_y_sums/μ_ij_sums_n
        @fastmath Γ_xy = dμ_ij_xy_sums/μ_ij_sums_n
        # s_n = μ_ij[n,:].*δ/μ_ij_sums_n

        calc_prodTerms_xy!(n,dS_xy,
                            X_mat,Y,s_n,
                            Γ_x,Γ_y,Γ_xy)
    end
    return nothing
end


function thD_calc!(dS_xyz::Vector{Float64},s_n::Vector{Float64},
                    μ_ij::Matrix{Float64},
                    X_mat::Matrix{Float64},Y::Matrix{Float64},Z::Matrix{Float64},
                    μ_ij_sums::Vector{Float64})
    dS_xyz[:] .= 0.0


    (N,K) = size(μ_ij)


    for n in 1:N
        μ_ij_sums_n = μ_ij_sums[n]

        dμ_ij_x_sums,dμ_ij_y_sums,dμ_ij_z_sums,dμ_ij_xy_sums,dμ_ij_xz_sums,dμ_ij_yz_sums,dμ_ij_xyz_sums =
            calc_derSums_xyz!(n,s_n,X_mat,Y,Z,μ_ij,μ_ij_sums_n)

        @fastmath Γ_x = dμ_ij_x_sums/μ_ij_sums_n
        @fastmath Γ_y = dμ_ij_y_sums/μ_ij_sums_n
        @fastmath Γ_z = dμ_ij_z_sums/μ_ij_sums_n
        @fastmath Γ_xy = dμ_ij_xy_sums/μ_ij_sums_n
        @fastmath Γ_xz = dμ_ij_xz_sums/μ_ij_sums_n
        @fastmath Γ_yz = dμ_ij_yz_sums/μ_ij_sums_n
        @fastmath Γ_xyz = dμ_ij_xyz_sums/μ_ij_sums_n
        # s_n = s_n./μ_ij_sums_n

        calc_prodTerms_xyz!(n,dS_xyz,
                    X_mat,Y,Z,s_n,Γ_x,Γ_y,Γ_z,Γ_xy,Γ_xz,Γ_yz,Γ_xyz)
    end
    return nothing
end


function grad_calc!(dS_x::Vector{T},
                    s_n::Vector{T},
                    μ_ij::Matrix{T},
                    X_mat::Matrix{Float64},
                    μ_ij_sums::Vector{T}) where T

    dS_x[:] .= 0.0

    (N,K) = size(μ_ij)

    for n in 1:N
        μ_ij_sums_n = μ_ij_sums[n]
        dμ_ij_x_sums = calc_derSums_x!(n,s_n,X_mat,μ_ij,μ_ij_sums_n)

        @fastmath Γ_x = dμ_ij_x_sums/μ_ij_sums_n
        # s_n = s_n./μ_ij_sums_n

        calc_prodTerms_x!(n,dS_x,X_mat,s_n,Γ_x)
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

function returnParameter!(q::Int64,X_mat::Matrix{Float64},
                        Z::Vector{Float64},X_0_t::Matrix{Float64},
                        X_t::Matrix{Float64},draws::Matrix{Float64},
                        F_t::SubArray{Float64,2},
                        γlen::Int64,β0len::Int64,βlen::Int64,σlen::Int64)
    (N,K) = size(X_mat)
    if q<0
        X_mat[:] .= 1.0
    elseif q<=γlen
        X_mat[:] .= Z[q]
    elseif q<=β0len
        for n in 1:N
            @inbounds X_mat[n,:] = X_t[q-γlen,:]
        end
    elseif q<=βlen
        # Characteristic Interactions
        for n in 1:N
            @inbounds X_mat[n,:] = X_t[1,:].*Z[q-β0len]
        end
    elseif q<=σlen
        #Quality Random Effect
        for n in 1:N,k in 1:K
            @inbounds X_mat[n,k] = draws[n,q-(βlen)]*X_0_t[q-(βlen),k]
        end
    else
        #Fixed Effect
        for n in 1:N
            @inbounds X_mat[n,:] = F_t[q-σlen,:]
        end
    end
    return Nothing
end
