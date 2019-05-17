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
                            ga_itr=500)

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
        wgt=spec_wgt)

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

    m = InsuranceLogit(c_data,haltonDim,nested=nested)
    println("Data Loaded")

    spec_labels = unpack_labels(m,spec_Dict)

    if ismissing(x_start)
        Istart = rand(m.parLength[:I])/10 .-.05
        βstart = rand(m.parLength[:β])/10 .-.05
        γstart = rand(m.parLength[:γ]*length(m.data._prodInteract))/10 .-.05
        σstart = rand(m.parLength[:σ])/10 .- .05
        FEstart = rand(m.parLength[:FE])/10 .-.05

        p0 = vcat(Istart,βstart,γstart,σstart,FEstart)
    else
        p0 = x_start
    end
    println("Begin Estimation")

    ## Estimate
    if method == "non_gradient"
        flag,fval,p_est = estimate_ng!(m,p0,method=:LN_NELDERMEAD)
    elseif method =="ga"
        p_est,fval,flag = gradient_ascent(m,p0,grad_tol = 1e-8,max_itr=ga_itr)
        println("Run NR Estimation")
        p_disp = p_est[1:20]
        println("Starting at $p_disp")
        p_est,fval,flag = newton_raphson_ll(m,p_est)
    elseif method=="nr"
        p_est, fval, flag = newton_raphson_ll(m,p0)
    end

    return m, p_est, spec_Dict, spec_labels, fval, flag

end




function res_process(model::InsuranceLogit,p_est::Vector{Float64})
    ##Total Population
    Pop = sum(weight(model.data).*choice(model.data))

    ## Create Param Dictionary
    paramFinal = parDict(model,p_est)
    println("Gradient-based Variance Estimate")
    Var_2 = calc_Avar(model,paramFinal)
    Var_2 = Var_2./Pop


    println("Hessian-based Variance Estimate")
    hess = Matrix{Float64}(undef,length(p_est),length(p_est))
    grad = Vector{Float64}(undef,length(p_est))
    res = log_likelihood!(hess,grad,model,p_est)
    Var_1 = -inv(hess)./Pop


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
                            ga_itr = 500)

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
                            ga_itr = ga_itr)

    println("Save Results")
    # Unpack
    m, p_est, spec_Dict, spec_labels, fval, flag = spec

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
