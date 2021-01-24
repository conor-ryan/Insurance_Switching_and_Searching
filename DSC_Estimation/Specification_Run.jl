function estimate_specification(df::DataFrame;
                            haltonDim = 1,
                            spec_per = [:hh_year_id],
                            spec_prd = [:product],
                            spec_ch = [:choice],
                            spec_ch_last = [:iplan],
                            spec_prodchr = [:padj,:iplan],
                            spec_prodchr_0= Vector{Symbol}(undef,0),
                            spec_inertchr=Vector{Symbol}(undef,0),
                            spec_demR=Vector{Symbol}(undef,0),
                            spec_prodInt=Vector{Symbol}(undef,0),
                            spec_fixInt=Vector{Symbol}(undef,0),
                            spec_fixEff=Vector{Symbol}(undef,0),
                            spec_wgt=[:constant],
                            nested = false,
                            x_start::Union{Missing,Vector{Float64}} = missing,
                            method="nr",
                            ga_itr=200,
                            ll_start = false,
                            use_active_var=false)

    ## Build Model
    c_data = ChoiceData(df;
        per = spec_per,
        prd = spec_prd,
        ch = spec_ch,
        ch_last = spec_ch_last,
        prodchr = spec_prodchr,
        prodchr_0= spec_prodchr_0,
        inertchr=spec_inertchr,
        demR=spec_demR,
        prodInt=spec_prodInt,
        fixInt=spec_fixInt,
        fixEff=spec_fixEff,
        wgt=spec_wgt,
        check_colin=true)

    if length(spec_prodInt)==0
        spec_prodInt = copy(spec_prodchr)
    end

    spec_Dict = Dict("per" => spec_per,
    "prd" => spec_prd,
    "ch" => spec_ch,
    "ch_last" => spec_ch_last,
    "prodchr" => spec_prodchr,
    "prodchr_0"=> spec_prodchr_0,
    "inertchr"=>spec_inertchr,
    "demR"=>spec_demR,
    "prodInt"=>spec_prodInt,
    "fixInt"=>spec_fixInt,
    "fixEff"=>spec_fixEff,
    "wgt"=>spec_wgt,
    "haltonDim"=>haltonDim)

    m = InsuranceLogit(c_data,haltonDim,nested=nested,use_active_var=use_active_var)
    println("Data Loaded")

    spec_labels = unpack_labels(m,spec_Dict)

    p0 = zeros(m.parLength[:All])
    if ll_start
        println("Compute Starting Point with No Heterogeneity")
        ll_ga = 0
        if length(spec_inertchr)>0
            ll_ga = 500
        end
        ll_res = estimate_specification(df,
                            haltonDim = 1,
                                spec_per = spec_per,
                                spec_prd = spec_prd,
                                spec_ch = spec_ch,
                                spec_ch_last = spec_ch_last,
                                spec_prodchr = spec_prodchr,
                                spec_prodchr_0= Vector{Symbol}(undef,0),
                                spec_inertchr=spec_inertchr,
                                spec_demR=spec_demR,
                                spec_prodInt=spec_prodInt,
                                spec_fixInt=spec_fixInt,
                                spec_fixEff=spec_fixEff,
                                spec_wgt=spec_wgt,
                                method = "ga",ga_itr = ll_ga,use_active_var=use_active_var)
        x_est = ll_res[2]
        ind1 = 1:(length(p0) - m.parLength[:FE] - m.parLength[:σ])
        ind2 = (length(p0) - m.parLength[:FE] + 1):m.parLength[:All]
        p0[ind1] = x_est[ind1]
        p0[ind2] = x_est[ind2 .- m.parLength[:σ]]
        indσ = (maximum(ind1) + 1):(minimum(ind2) - 1)
        p0[indσ].= rand(m.parLength[:σ])/10 .+ 1e-4
    elseif ismissing(x_start)
        Istart = rand(m.parLength[:I])/10 .-.05
        βstart = rand(m.parLength[:β])/10 .-.05
        γstart = rand(m.parLength[:γ]*length(m.data._prodInteract))/10 .-.05
        σstart = rand(m.parLength[:σ])/10 .+ 1e-5
        FEstart = rand(m.parLength[:FE])/10 .-.05

        p0 = vcat(Istart,βstart,γstart,σstart,FEstart)
    else
        p0[:] = x_start[:]
        println(log_likelihood(m,x_start))
    end
    println("Begin Estimation")
    println(log_likelihood(m,p0))
    ## Estimate
    if method == "non_gradient"
        flag,fval,p_est = estimate_ng!(m,p0,method=:LN_NELDERMEAD)
    elseif method =="ga"
        if ga_itr>1
            p_est,fval,flag = gradient_ascent(m,p0,grad_tol = 1e-8,max_itr=ga_itr)
        else
            p_est = p0
        end
        println("Run NR Estimation")
        p_disp = p_est[1:20]
        println("Starting at $p_disp")
        p_est,fval = newton_raphson_ll(m,p_est)
    elseif method=="nr"
        p_est, fval = newton_raphson_ll(m,p0)
    end

    return m, p_est, spec_Dict, spec_labels, fval
end




function res_process(model::InsuranceLogit,p_est::Vector{Float64})
    ##Total Population
    Pop = length(model.data._personIDs)

    ## Create Param Dictionary
    # paramFinal = parDict(model,p_est)
    # println("Gradient-based Variance Estimate")
    # Var_2 = calc_Avar(model,paramFinal)
    # Var_2 = Var_2./Pop


    println("Hessian-based Variance Estimate")
    hess = Matrix{Float64}(undef,length(p_est),length(p_est))
    grad = Vector{Float64}(undef,length(p_est))
    res = log_likelihood!(hess,grad,model,p_est)
    Var_1 = -inv(hess)./Pop
    Var_2 = Var_1


    if any(diag(Var_2.<0)) | any(diag(Var_1.<0))
        println("Some negative variances")
        stdErr_1 = sqrt.(abs.(diag(Var_1)))
        stdErr_2 = sqrt.(abs.(diag(Var_2)))
    else
        stdErr_1 = sqrt.(diag(Var_1))
        stdErr_2 = sqrt.(diag(Var_2))
    end
    t_stat = p_est./stdErr_1

    stars = Vector{String}(undef,length(t_stat))
    for i in 1:length(stars)
        if abs(t_stat[i])>2.326
            stars[i] = "***"
        elseif abs(t_stat[i])>1.654
            stars[i] = "**"
        elseif abs(t_stat[i])>1.282
            stars[i] = "*"
        else
            stars[i] = ""
        end
    end

    return Var_1, stdErr_1, stdErr_2, t_stat, stars
end

function MainSpec(df::DataFrame,filename::String;
                            haltonDim = 1,
                            spec_per = [:hh_year_id],
                            spec_prd = [:product],
                            spec_ch = [:choice],
                            spec_ch_last = [:iplan],
                            spec_prodchr = [:padj,:iplan],
                            spec_prodchr_0= Vector{Symbol}(undef,0),
                            spec_inertchr=Vector{Symbol}(undef,0),
                            spec_demR=Vector{Symbol}(undef,0),
                            spec_prodInt=Vector{Symbol}(undef,0),
                            spec_fixInt=Vector{Symbol}(undef,0),
                            spec_fixEff=Vector{Symbol}(undef,0),
                            spec_wgt=[:constant],
                            nested = false,
                            x_start::Union{Missing,Vector{Float64}} = missing,
                            method="nr",
                            ga_itr = 200,
                            ll_start = false,
                            use_active_var=false)

    println("Estimate Specification")
    spec = estimate_specification(df,
                        haltonDim = haltonDim,
                            spec_per = spec_per,
                            spec_prd = spec_prd,
                            spec_ch = spec_ch,
                            spec_ch_last = spec_ch_last,
                            spec_prodchr = spec_prodchr,
                            spec_prodchr_0= spec_prodchr_0,
                            spec_inertchr=spec_inertchr,
                            spec_demR=spec_demR,
                            spec_prodInt=spec_prodInt,
                            spec_fixInt=spec_fixInt,
                            spec_fixEff=spec_fixEff,
                            spec_wgt=spec_wgt,
                            x_start = x_start,
                            method = method,
                            ga_itr = ga_itr,
                            ll_start = ll_start,
                            use_active_var=use_active_var)

    println("Save Results")
    # Unpack
    m, p_est, spec_Dict, spec_labels, fval = spec

    file = "$(homedir())/Documents/Research/CovCAInertia/Output/Estimation_Results/$filename.jld2"
    @save file p_est spec_Dict fval

    println("Calculate Standard Errors")
    Var, se1, se2,t_stat, stars = res_process(m,p_est)
    out1 = DataFrame(labels=spec_labels,pars=p_est,se=se1,ts=t_stat,sig=stars)
    file1 = "$(homedir())/Documents/Research/CovCAInertia/Output/Estimation_Results/$filename.csv"
    CSV.write(file1,out1)

    return p_est, fval
end

function unpack_labels(m::InsuranceLogit,spec_Dict::Dict{String,Any})
    labels = Vector{String}(undef,m.parLength[:All])

    ind = 1
    for i in 1:length(spec_Dict["inertchr"])
        search_spec = String(spec_Dict["inertchr"][i])
        labels[ind] = "Search::$search_spec"
        ind+=1
    end

    for i in 1:length(spec_Dict["prodchr"])
        labels[ind] = String(spec_Dict["prodchr"][i])
        ind+=1
    end

    for i in 1:length(spec_Dict["demR"])
        dem_spec = String(spec_Dict["demR"][i])
        for j in 1:length(spec_Dict["prodInt"])
            prod_spec = String(spec_Dict["prodInt"][j])
            labels[ind] = "$prod_spec::$dem_spec"
            ind+=1
        end
    end

    for i in 1:length(spec_Dict["prodchr_0"])
        prod_spec = String(spec_Dict["prodchr_0"][i])
        labels[ind] = "Var::$prod_spec"
        ind+=1
    end

    for i in 1:m.parLength[:FE]
        labels[ind] = "FE_$i"
        ind+=1
    end
    return labels
end


function predict_switching(m::InsuranceLogit,p_vec::Vector{Float64},spec::Dict{String,Any};
    fullAtt::Bool=false,
    noCont::Bool=false,
    noHass::Bool=false,
    useActiveVar::Bool=false)

    p_est = copy(p_vec)

    Ilength = m.parLength[:I]

    parBase = parDict(m,p_est)
    if useActiveVar
        parBase.ω_i[:] = m.data.active[:]
    end
    test = mean(parBase.ω_i)
    println("Test Inattention 1 : $test")

    # cont_pars = vcat(Ilength .+ (2:4),(Ilength.+ m.parLength[:β]).+ vcat((2:4), 4 .+ (2:4),8 .+ (2:4),12 .+ (2:4)))
    # cont_pars = vcat(Ilength .+ (3:4),(Ilength.+ m.parLength[:β]).+ vcat((3:4), 4 .+ (3:4),8 .+ (3:4),12 .+ (3:4)))
    # hassle_pars = vcat([Ilength + 2],(Ilength.+ m.parLength[:β]).+ [2, 6,10,14])
    if noCont
        β_ind = inlist(spec["prodchr"],[:inet,:iiss])
        parBase.β[β_ind,:].=0.0
        parBase.β_0[β_ind].= 0.0
    end
    if noHass
        β_ind = findall(spec["prodchr"].==:iplan)
        parBase.β[β_ind,:].=0.0
        parBase.β_0[β_ind].= 0.0
    end



    individual_values!(m,parBase)

    if fullAtt
        parBase.ω_i[:] .= 1.0
    end
    test = mean(parBase.ω_i)
    println("Test Inattention 2 : $test")
    individual_shares(m,parBase)

    pred, data = count_switchers(m,parBase)

    return pred, data
end

function count_switchers(m::InsuranceLogit,par::parDict{Float64})
    # inertPlan = choice_last(m.data)
    # obsPlan = choice(m.data)

    All_Return = 0.0
    All_Stay = 0.0
    All_Stay_Obs = 0.0
    for app in eachperson(m.data)
        i = person(app)[1]
        # println(i)
        idx = m.data._personDict[i]
        stay_prob,inertPlan,obsPlan = calc_switch_prob(app,m,par)
        # for (yr, val) in stay_prob
        #     if isnan(val)
        #         println("Found Nan")
        #     end
        # end
        # years = sort(Int.(keys(app._personYearDict[ind])))
        for (yr,per_idx_yr) in app._personYearDict[i]
            # idx_yr = idx[per_idx_yr]
            returning = sum(inertPlan[per_idx_yr])
            if returning==0.0
                continue
            end
            if (isnan(stay_prob[yr]))
                println("NaN Value for person $i in year $yr")
                return "Break"
            end
            All_Return +=  1.0
            All_Stay += stay_prob[yr]
            retplan = findall(inertPlan[per_idx_yr].>0)
            obs = findall(obsPlan[per_idx_yr].>0)
            if retplan==obs
                All_Stay_Obs += 1.0
            end
        end
    end
    data = All_Stay_Obs/All_Return
    pred = All_Stay/All_Return
    return pred, data
end

function returning_index(m::InsuranceLogit,par::parDict{Float64})
    # inertPlan = choice_last(m.data)
    # obsPlan = choice(m.data)

    ret_idx = Vector{Int}(undef,0)
    ylast_large = choice_last(m.data)[:]
    for app in eachperson(m.data)
        i = person(app)[1]
        # println(i)
        idx = m.data._personDict[i]
        # stay_prob,inertPlan,obsPlan = calc_switch_prob(app,m,par)
        inertPlan = ylast_large[idx]
        for (yr,per_idx_yr) in app._personYearDict[i]
            # idx_yr = idx[per_idx_yr]
            returning = sum(inertPlan[per_idx_yr])
            if returning==0.0
                continue
            end
            ret_idx = vcat(ret_idx,idx[per_idx_yr])
        end
    end
    return ret_idx
end

function returning_peryear_index(m::InsuranceLogit,par::parDict{Float64})
    # inertPlan = choice_last(m.data)
    # obsPlan = choice(m.data)

    ret_idx = Vector{Int}(undef,0)
    ylast_large = choice_last(m.data)[:]
    for app in eachperson(m.data)
        i = person(app)[1]
        # println(i)
        idx = m.data._personDict[i]
        # stay_prob,inertPlan,obsPlan = calc_switch_prob(app,m,par)
        inertPlan = ylast_large[idx]
        years = sort(Int.(keys(m.data._personYearDict[i])))
        yr_ind = 0
        ω_index = m.data._searchDict[i]
        for year in years
            per_idx_yr = m.data._personYearDict[i][year]
            yr_ind+=1
            # idx_yr = idx[per_idx_yr]
            returning = sum(inertPlan[per_idx_yr])
            if returning==0.0
                continue
            end
            ret_idx = vcat(ret_idx,ω_index[yr_ind])
        end
    end
    return ret_idx
end



function activePredict(m::InsuranceLogit,par::parDict{Float64},df::DataFrame,spec,rundate)
    inertPlan = choice_last(m.data)
    obsPlan = choice(m.data)
    active_long = df[:active]
    active_obs = Vector{Float64}(undef,length(par.ω_i))
    active_pred = Vector{Float64}(undef,length(par.ω_i))
    returning = Vector{Float64}(undef,length(par.ω_i))
    stay = zeros(length(par.ω_i))
    stay_pred = zeros(length(par.ω_i))
    personID = Vector{Float64}(undef,length(par.ω_i))
    yearID = Vector{Float64}(undef,length(par.ω_i))

    ind = 0
    for app in eachperson(m.data)
        i = person(app)[1]
        yr_ind = 0
        idx = m.data._personDict[i]
        stay_prob,inertPlan,obsPlan = calc_switch_prob(app,m,par)
        s_idx = app._searchDict[i]
        dict = app._personYearDict[i]
        years = sort(Int.(keys(dict)))
        for yr in years
            per_idx_yr = dict[yr]
            ind+=1
            yr_ind+=1
            idx_yr = idx[per_idx_yr]
            returning[ind] = sum(inertPlan[per_idx_yr])
            active_obs[ind] = active_long[idx_yr[1]]
            active_pred[ind] = par.ω_i[s_idx[yr_ind]]

            retplan = findall(inertPlan[per_idx_yr].>0)
            obs = findall(obsPlan[per_idx_yr].>0)
            if retplan==obs
                stay[ind] = 1.0
            end
            stay_pred[ind] = stay_prob[yr]
            personID[ind] = i
            yearID[ind] = yr
        end
    end
    out1 = DataFrame(Person = personID, Year = yearID,
                        stay_obs=stay,stay_pred=stay_pred,
                        returning=returning,active_obs=active_obs,active_pred=active_pred)
    file1 = "$(homedir())/Documents/Research/CovCAInertia/Output/Estimation_Results/active_$spec$rundate.csv"
    CSV.write(file1,out1)

    return nothing
end



function calc_switch_prob(app::ChoiceData,d::InsuranceLogit,p::parDict{T}) where T

        ind, S_ij, wgt, idxitr, X_t, X_0_t,X_int, Z, F_t, X_last, y_last = unPackChars(app,d)
        wgt = convert(Array{Float64,2},wgt)
        S_ij = convert(Array{Float64,2},S_ij)
        ω_i = p.ω_i[app._searchDict[ind]]
        choice_ind = findall(S_ij[:].==1)

        draws = d.draws

        # Get Utility and derivative of Utility
        μ_ij,s_hat,s_uncond = unPackParChars(p,idxitr)

        # Year Indices
        yr_next = findYearInd(app._personYearDict[ind])

        # Pre-Calculate Shares
        # μ_ij_sums = preCalcμ(μ_ij,app._personYearDict[ind])
        s_hat, s_n, ll_n, μ_ij_sums = calc_shares_mat(μ_ij,ω_i,y_last,choice_ind,yr_next)

        stay_prob = Dict{Int64,Float64}()
        years = sort(Int.(keys(app._personYearDict[ind])))
        last_choice = Vector{Int64}(undef,0)
        for (i,yr) in enumerate(years)
            idx_yr = app._personYearDict[ind][yr]
            if sum(y_last[idx_yr])==0.0
                last_choice = vcat(last_choice,choice_ind[i])
                stay_prob[yr] = 0.0
                continue
            end
            ret_choice = idx_yr[findall(y_last[idx_yr].>0)]
            if (yr-1) in years
                choice_seq = vcat(last_choice,ret_choice)
                denom = mean(prod(s_hat[last_choice,:],dims=1))
                stay_prob[yr] = mean(prod(s_hat[choice_seq,:],dims=1))/denom
                last_choice = vcat(last_choice,choice_ind[i])
            else
                last_choice = Vector{Int64}(undef,0)
                stay_prob[yr] = mean(s_hat[ret_choice,:])
            end
        end

    return stay_prob,y_last,S_ij
end
