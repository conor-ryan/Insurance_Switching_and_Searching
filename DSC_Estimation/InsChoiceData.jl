using DataFrames
import Base.getindex, Base.setindex!

abstract type
    ModelData
end


# Make a type to interface with the data for discrete choice models
struct ChoiceData <: ModelData
    # Matrix of the data (transposed, pre-sorted)
    data::Matrix{Float64}
    # Matrix of the product level data (pre-sorted)
    # pdata::DataFrame
    # Matrix of Fixed Effects
    fixedEffects::Matrix{Float64}
    # Index of the data column names
    index
    # Names of rows (columns of input data)
    prodchars   # Product Characteristics
    prodchars_0   # Product Characteristics
    choice      # Binary choice indicator
    demoRaw    # Household Demographics - raw
    wgt     # Number of People in each type
    # unins     # Outside Option Share

    # Precomputed Indices
    _person::Array{Int,1}
    _product::Array{Int,1}
    _panel::Array{Int,1}
    _prodchars::Array{Int,1}
    _prodchars_0::Array{Int,1}
    _inertchars::Array{Int,1}
    _choice::Array{Int,1}
    _choice_last::Array{Int,1}
    _demoRaw::Array{Int,1}
    _wgt::Array{Int,1}

    # Random Coefficient Specification
    _randCoeffs::Array{Int,1}
    _prodInteract::Array{Int,1}

    # ID Lookup Mappings
    _personIDs::Array{Float64,1}
    _personDict::Dict{Int, UnitRange{Int}}
    _personYearDict::Dict{Int, Dict{Int,UnitRange}}
    _searchDict::Dict{Int, UnitRange{Int}}
    _productDict::Dict{Int, Array{Int,1}}

    _rel_fe_Dict::Dict{Real,Array{Int64,1}}
end

function ChoiceData(data_choice::DataFrame;
        per=[:hh_id],
        panel=[:year],
        prd=[:product],
        prodchr=Vector{Symbol}(undef,0),
        prodchr_0=Vector{Symbol}(undef,0),
        inertchr= Vector{Symbol}(undef,0),
        ch=[:choice],
        ch_last = [:iplan],
        demR=Vector{Symbol}(undef,0),
        prodInt=Vector{Symbol}(undef,0),
        # demoRaw=[:F0_Y0_LI1,
        #          :F0_Y1_LI0,:F0_Y1_LI1,
        #          :F1_Y0_LI0,:F1_Y0_LI1,
        #          :F1_Y1_LI0,:F1_Y1_LI1],
        fixEff=Vector{Symbol}(undef,0),
        fixInt=Vector{Symbol}(undef,0),
        wgt=[:constant])

    # Get the size of the data
    n, k = size(data_choice)

    # Convert everything to an array once for performance
    i = convert(Matrix{Float64},data_choice[per])
    yr = convert(Matrix{Float64},data_choice[panel])
    j = convert(Matrix{Float64},data_choice[prd])
    X = convert(Matrix{Float64},data_choice[prodchr])
    # X_0 = convert(Array{Float64},data_choice[prodchars_0])
    In = convert(Matrix{Float64},data_choice[inertchr])
    y = convert(Matrix{Float64},data_choice[ch])
    y_last = convert(Matrix{Float64},data_choice[ch_last])
    Z = convert(Matrix{Float64},data_choice[demR])
    w = convert(Matrix{Float64},data_choice[wgt])


    println("Create Fixed Effects")
    F, feNames = build_FE(data_choice,fixEff,hasConstant=true)
    F = permutedims(F,(2,1))

    F_int, feNames_int = build_FE(data_choice,fixInt,hasConstant=true)
    F_int = permutedims(F_int,(2,1))



    index = Dict{Symbol, Int}()
    dmat = Matrix{Float64}(undef,n,0)


    # Create a data matrix, only including person id
    println("Put Together Data non FE data together")
    k = 0
    for (d, var) in zip([i,yr,j,X, Z, y,In,y_last,w], [per,panel,prd,prodchr,
        demR, ch,inertchr,ch_last, wgt])
        for l=1:size(d,2)
            k+=1
            dmat = hcat(dmat, d[:,l])
            index[var[l]] = k
        end
    end

    println("Build Interaction Terms")
    if length(X)==0
        # X_int = ones(n,1)
        X_int = permutedims(F,(2,1))
    else
        # X_int = hcat(X,ones(n))
        X_int = hcat(X,permutedims(F,(2,1)))
    end

    L = size(X_int,2)
    M = size(F_int,1)
    println("$n by $L*$M")
    interactions = Matrix{Float64}(undef,L*M,n)

    for l in 1:L, m in 1:M
        ind = (l-1)*M+m
        # int_vec = Matrix{Float64}(undef,1,n)
        for i in 1:n
            @inbounds @fastmath interactions[ind,i] = F_int[m,i]*X_int[i,l]
        end
        # keep_int = sum(y[findall(abs.(int_vec[:]).>0)])>1.0
        # if keep_int
        #     interactions = vcat(interactions,int_vec)
        # end
    end

    ### Drop 0 valued Interaction Terms
    # M = size(interactions,1)
    # keep_int = BitArray{1}(undef,M)
    # for m in 1:M
    #     keep_int[m] = sum(y[findall(abs.(interactions[m,:]).>0)])>1.0
    # end
    # interactions = interactions[keep_int,:]

    # positive_values = sum(abs.(interactions),dims=2)[:]
    # interactions = interactions[positive_values.>0,:]
    println("Concatenate Interactions")
    F = vcat(F,interactions)

    ### Drop 0 valued Fixed Effects
    zero_ind = findall(sum(F,dims=2)[:].==0)
    F_ind = 1:size(F,1)
    for k in zero_ind
        F_ind = F_ind[F_ind.!=k]
    end
    F = F[F_ind,:]
    # F_all = interactions



    #Transpose data to store as rows
    println("Transpose")
    dmat = permutedims(dmat,(2,1))
    i = permutedims(i,(2,1))
    j = permutedims(j,(2,1))

    # Precompute the row indices
    _person = getDictArray(index,per)
    _panel = getDictArray(index,panel)
    _product = getDictArray(index,prd)
    _prodchars = getDictArray(index, prodchr)
    _prodchars_0 = getDictArray(index, prodchr_0)
    _inertchars = getDictArray(index, inertchr)
    _choice = getDictArray(index, ch)
    _choice_last =  getDictArray(index, ch_last)
    _demoRaw = getDictArray(index, demR)
    _wgt = getDictArray(index, wgt)

    println("Check Collinearity")
    if length(_inertchars)>0
        all_ind = vcat(_inertchars)
        all_data = dmat[all_ind,:]
        X = all_data*all_data'
        smallest_ev = minimum(abs.(eigvals(X)))
        println("Smallest Data Eigenvalue (search): $smallest_ev")
    else
        println("Smallest Data Eigenvalue (search):No Search Variable")
    end

    smallest_ev = 0
    while abs(smallest_ev).<1e-10
        all_ind = vcat(_choice,_prodchars,_demoRaw)
        all_data = vcat(dmat[all_ind,:],F)
        X = all_data*all_data'
        all_vals = abs.(eigvals(X))
        num_zero_vals = length(findall(all_vals.<1e-10))
        smallest_ev = minimum(all_vals)
        println("Smallest Data Eigenvalue (choice): $smallest_ev")
        println("Number of Zero Eigenvalues (choice): $num_zero_vals")
        pdim = size(dmat[all_ind,:],1)
        if abs(smallest_ev)<1e-10
            v = eigvecs(X)
            zero_vals = findall(all_vals.<1e-10)
            drop_list = Vector{Int}(undef,0)
            for ind_v in zero_vals
                ind_colin = maximum(findall(abs.(v[:,ind_v]).>1e-10))
                ind_colin_fe = ind_colin - pdim
                if ind_colin<pdim
                    println("Collinearity in Product Characteristics")
                    break
                elseif ind_colin_fe in drop_list
                    continue
                else
                    drop_list = vcat(drop_list,[ind_colin_fe])
                end
            end
            F_ind = 1:size(F,1)
            for k in drop_list
                F_ind = F_ind[F_ind.!=k]
            end
            F = F[F_ind,:]
            println("Dropping: $drop_list")
        end
    end

    #### Check Constant
    all_ind = vcat(_choice,_prodchars,_demoRaw)
    all_data = vcat(dmat[all_ind,:],F,ones(1,n))
    X = all_data*all_data'
    all_vals = abs.(eigvals(X))
    num_zero_vals = length(findall(all_vals.<1e-10))
    smallest_ev = minimum(all_vals)
    println("Constant Check (choice): $smallest_ev")


    ## Rand Coefficient Index
    _randCoeffs = Array{Int,1}(undef,length(prodchr_0))
    for (i,var) in enumerate(prodchr_0)
        _randCoeffs[i] = findall(var.==prodchr)[1]
    end

    ## Product Interaction Index
    if length(prodInt)==0
        _prodInteract = 1:length(prodchr)
    else
        _prodInteract = Array{Int,1}(undef,length(prodInt))
        for (i,var) in enumerate(prodInt)
            _prodInteract[i] = findall(var.==prodchr)[1]
        end
    end
    println("Product Interaction Terms: $_prodInteract")

    # Get Person ID Dictionary Mapping for Easy Subsets
    println("Person ID Mapping")
    _personDict = Dict{Real, UnitRange{Int}}()
    _searchDict = Dict{Real, UnitRange{Int}}()
    _personYearDict = Dict{Int, Dict{Int,UnitRange{Int}}}()
    allids = dmat[_person,:][1,:]
    allyears = dmat[_panel,:][1,:]
    uniqids = sort(unique(allids))
    last_ind = 1

    for id in uniqids
        idx1 = searchsortedfirst(allids,id)
        idxJ = searchsortedlast(allids,id)
        _personDict[id] = idx1:idxJ
        per_years = allyears[idx1:idxJ]
        uniq_years = sort(unique(per_years))
        _panelDict = Dict{Int, UnitRange{Int}}()
        for yr in uniq_years
            idx1 = searchsortedfirst(per_years,yr)
            idxJ = searchsortedlast(per_years,yr)
            _panelDict[yr] = idx1:idxJ
        end
        _personYearDict[id] = _panelDict
        new_ind = last_ind + length(uniq_years)
        _searchDict[id] = (last_ind+1):new_ind
        last_ind = copy(new_ind)
    end


    #Create Product Dictionary
    println("Product Dictionary")
    prod_vec = j[1,:]
    _productDict = build_ProdDict(prod_vec)
    # allprods = sort(unique(j))
    # _productDict = Dict{Real, Array{Int}}()
    # for id in allprods
    #     _productDict[id] = findin(j,id)
    # end

    # Relevant Parameters Per Person
    println("Find Relevant Parameters")
    rel_fe_Dict = Dict{Real,Array{Int64,1}}()
    for (id,idxitr) in _personDict
        F_t = view(F,:,idxitr)
        any_positive = maximum(abs.(F_t),dims=2)[:,1]
        pars_relevant = findall(any_positive .>0)
        rel_fe_Dict[id] = pars_relevant
    end

    # Make the data object
    m = ChoiceData(dmat, F, index,
            prodchr,prodchr_0,ch, demR,wgt,
             _person,_panel,_product, _prodchars,_prodchars_0,_inertchars,
            _choice,_choice_last, _demoRaw, _wgt,
             _randCoeffs,_prodInteract,
             uniqids,_personDict,_personYearDict,_searchDict,
             _productDict,
            rel_fe_Dict)
    return m
end


function build_ProdDict(j::Array{T,N}) where {T,N}
    allprods = unique(j)
    sort!(allprods)
    _productDict = Dict{Real, Array{Int64,1}}()
    for id in allprods
        _productDict[id] = findall(j.==id)
    end
    return _productDict
end

# function build_FE(data_choice::DataFrame,fe_list::Vector{T}) where T
#     # Create Fixed Effects
#     n, k = size(data_choice)
#     L = 0
#
#     # No Fixed effects for empty lists
#     if typeof(fe_list)!=Vector{Symbol}
#         println("No Fixed Effects")
#         F = Matrix{Float64}(undef,n,L)
#         feNames = Vector{Symbol}(undef,0)
#         return F,feNames
#     end
#
#     for fe in fe_list
#         fac_variables = data_choice[fe]
#         factor_list = sort(unique(fac_variables))
#         # if fe==:constant
#         #     num_effects=1
#         # elseif (!(:constant in fe_list)) & (fe==fe_list[1])
#         #     num_effects = length(factor_list)
#         #     # if fe==:Market
#         #     #     num_effects = length(factor_list) - 3
#         #     # end
#         # else
#         num_effects = length(factor_list)-1
#         # end
#         L+=num_effects
#     end
#
#     F = zeros(n,L)
#     feNames = Vector{Symbol}(undef,0)
#     ind = 1
#     for fe in fe_list
#         if fe==:constant
#             F[:,ind] = 1
#             ind+=1
#             continue
#         end
#         fac_variables = data_choice[fe]
#         factor_list = sort(unique(fac_variables))
#         filter!(!ismissing,factor_list)
#         # if (!(:constant in fe_list)) & (fe==fe_list[1])
#         #     st_ind = 1
#         # else
#         st_ind = 2
#         # end
#
#         for fac in factor_list[st_ind:length(factor_list)]
#             # fac_data = zeros(n)
#             # fac_data[fac_variables.==fac] = 1.0
#             # if fac in ["ND_4","MD_4","IA_7"]
#             #     continue
#             # end
#             fac_index = fac_variables.==fac
#             fac_index[ismissing.(fac_index)] .= false
#             fac_index = Bool.(fac_index)
#             F[fac_index,ind] .= 1
#             ind+= 1
#
#             feNames = vcat(feNames,Symbol(fac))
#         end
#     end
#     return F, feNames
# end


function build_FE(data_choice::DataFrame,fe_list::Vector{T};bigFirm=false,hasConstant=false) where T
    # Create Fixed Effects
    n, k = size(data_choice)
    L = 0

    # No Fixed effects for empty lists
    if typeof(fe_list)!=Vector{Symbol}
        println("No Fixed Effects")
        F = Matrix{Float64}(undef,n,L)
        feNames = Vector{Symbol}(undef,0)
        return F,feNames
    end

    for fe in fe_list
        fac_variables = data_choice[fe]
        factor_list = sort(unique(fac_variables[map(!,ismissing.(fac_variables))]))
        if fe==:constant
            num_effects=1
        elseif (!((:constant in fe_list) | hasConstant)) & (fe==fe_list[1])
            num_effects = length(factor_list)
        else
            num_effects = length(factor_list)-1
        end
        L+=num_effects
    end

    F = zeros(n,L)
    feNames = Vector{Symbol}(undef,0)
    ind = 1
    for fe in fe_list
        if fe==:constant
            F[:,ind] = 1
            ind+=1
            continue
        end
        fac_variables = data_choice[fe]
        factor_list = sort(unique(fac_variables[map(!,ismissing.(fac_variables))]))
        if (!((:constant in fe_list) | hasConstant)) & (fe==fe_list[1])
            st_ind = 1
        else
            println("Skip")
            st_ind = 2
        end

        for fac in factor_list[st_ind:length(factor_list)]
            F[map(!,ismissing.(fac_variables)) .& (fac_variables.==fac),ind] .= 1
            ind+= 1

            feNames = vcat(feNames,Symbol(fac))
        end
    end
    return F, feNames
end

# Defining Indexing Methods on ChoiceData
Symbols = Union{Symbol, Vector{Symbol}}
getindex(m::ChoiceData, idx) = m.data[idx,:]
getindex(m::ChoiceData, idx, cols) = m.data[idx, cols]
getindex(m::ChoiceData, idx::Array{Int,1}) = m.data[idx,:]
getindex(m::ChoiceData, idx::Array{Int,1}, cols::Array{Int,1}) = m.data[idx, cols]
getindex(m::ChoiceData, idx::Symbols) = m.data[getindex.(m.index, idx),:]
getindex(m::ChoiceData, idx::Symbols, cols) = m.data[getindex.(m.index, idx),cols]

# Define other retrieval methods on ChoiceData
person(m::ChoiceData)      = m[m._person]
product(m::ChoiceData)      = m[m._product]
@inbounds prodchars(m::ChoiceData)   = m[m._prodchars]
@inbounds prodchars0(m::ChoiceData)   = m[m._prodchars_0]
inertchars(m::ChoiceData)   = m[m._inertchars]
choice(m::ChoiceData)      = m[m._choice]
choice_last(m::ChoiceData)      = m[m._choice_last]
@inbounds demoRaw(m::ChoiceData)     = m[m._demoRaw]
weight(m::ChoiceData)      = m[m._wgt]


fixedEffects(m::ChoiceData)= m.fixedEffects
fixedEffects(m::ChoiceData,idx)= view(m.fixedEffects,:,idx)

# fixedInt(m::ChoiceData)= m.fixedInteractions
# fixedInt(m::ChoiceData,idx)= view(m.fixedInteractions,:,idx)

########################################################################
#################### Iterating over People ############################
########################################################################

# Quickly Generate Subsets on People
function subset(d::T, idx) where T<:ModelData

    data = d.data[:,idx]
    fixedEf = d.fixedEffects
    #fixedEf = view(d.fixedEffects,:,idx)
#    people = data[d._person,:]

    # Don't subset any other fields for now...
    return T(data,
    fixedEf,
    # d.fixedInteractions,
    # Index of the column names
    d.index,
    # Names of rows (columns of input data)
    d.prodchars,   # Product Characteristics
    d.prodchars_0,   # Product Characteristics
    d.choice,      # Binary choice indicator
    d.demoRaw,    # Household Demographics - raw
    d.wgt,     # Demographic Fixed Effects
    # Precomputed Indices
    d._person,
    d._panel,
    d._product,
    d._prodchars,
    d._prodchars_0,
    d._inertchars,
    d._choice,
    d._choice_last,
    d._demoRaw,
    d._wgt,
    d._randCoeffs,
    d._prodInteract,
    d._personIDs,
    d._personDict,
    d._personYearDict,
    d._searchDict,
    d._productDict,
    d._rel_fe_Dict)
end

########## People Iterator ###############
# Define an Iterator Type
mutable struct PersonIterator
    data
    id
end

# Construct an iterator to loop over people
function eachperson(m::ChoiceData)
    #ids = sort(unique(person(m)))
    ids = m._personIDs
    return PersonIterator(m, ids)
end

# Iterate Code in 0.6
# start(itr::PersonIterator) = 1
# function next(itr::PersonIterator, state)
#     # Get the current market
#     id = itr.id[state]
#
#     # Find which indices to use
#     idx = itr.data._personDict[id]
#
#     # Subset the data to just look at the current market
#     submod = subset(itr.data, idx)
#
#     return submod, state + 1
# end
# done(itr::PersonIterator, state) = state > length(itr.id)

function Base.iterate(iter::PersonIterator, state=1)

    if state> length(iter.id)
        return nothing
    end

    # Get the current market
    id = iter.id[state]

    # Find which indices to use
    idx = iter.data._personDict[id]

    # Subset the data to just look at the current market
    submod = subset(iter.data, idx)

    return (submod, state + 1)
end


###########################################################
### Model Object ########


abstract type LogitModel end

mutable struct InsuranceLogit <: LogitModel
    # Dictionary of Parameters and implied lengths
    parLength::Dict{Symbol, Int64}
    # ChoiceData struct
    data::ChoiceData

    #Store Halton Draws
    draws::Array{Float64,2}
end


# Construct the model instance
function InsuranceLogit(c_data::ChoiceData,haltonDim::Int;
    nested=false)

    # Get Parameter Lengths
    println(size(demoRaw(c_data)))
    println(c_data.demoRaw)
    println(c_data._demoRaw)

    println(size(prodchars(c_data)))
    println(c_data.prodchars)
    println(c_data._prodchars)

    γlen = size(demoRaw(c_data),1)
    β0len = size(prodchars0(c_data),1)
    βlen = size(prodchars(c_data),1)
    intLen = length(c_data._prodInteract)
    Ilen = size(inertchars(c_data),1)
    flen = size(fixedEffects(c_data),1)

    println(γlen)
    println(βlen)

    if haltonDim<=1 & !nested
        σlen = 0
    elseif (haltonDim>1) & (!nested)
        # σlen = (size(prodchars(c_data),1)-1)
        σlen = β0len
    elseif (haltonDim<=1) & nested
        σlen =1
    else
        error("Nesting Parameter not right")
        return
    end

    #total = 1 + γlen + β0len + γlen + flen + σlen
    # total = γlen + βlen + γlen + flen + σlen
    # total = βlen + γlen + Ilen + flen + σlen
    total = Ilen + βlen + intLen*γlen + σlen + flen
    parLength = Dict(:γ=>γlen,:β=>βlen,:I=>Ilen,:FE=>flen,
    :σ => σlen,:βγ=>intLen,:All=>total)

    println(parLength)

    # Initialize Halton Draws
    # These are the same across all individuals
    #draws = permutedims(MVHaltonNormal(haltonDim,2),(2,1))
    if haltonDim>1
        draws = MVHaltonNormal(haltonDim,σlen;scrambled=false)
    else
        dim_min = max(βlen,1)
        draws = zeros(1,dim_min)
    end


    d = InsuranceLogit(parLength,
                        c_data,
                        draws)
    return d
end
