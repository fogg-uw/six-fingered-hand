#= Estimate (with simulations) and calculate (exactly) the expected CFs
   of all quartets, for all simulated networks within a "job", where
   1 job = 1 scenario = 1 combination of parameter choices

usage: navigate to the directory "jobi" with i = index of scenario
       where the input files should exist.
input files: ./SiPhyNetwork_output/sim*.tree
output file: ./quartets.csv

usage: julia ../2_extract_quartet_subnetworks.jl seed ngt rho delete1

options:
ngt = number of gene trees for PCS to simulate per quartet network.
delete1 = 1,0 is whether to delete 1 (one) leaf from each network before doing anything else with it (like examining quartets).
=#

using PhyloNetworks
using Random
using DataFrames
using Combinatorics
using CSV

include("find_anomalies.jl")

seed    = parse(Int64, ARGS[1])
ngt     = parse(Int64, ARGS[2])
rho     = parse(Float64, ARGS[3])
delete1 = parse(Int8,  ARGS[4])

if !(delete1 == 0 || delete1 == 1)
	print("delete1 not 0 or 1, interpreting as 0")
	delete1 = 0
end

jobid = replace(basename(pwd()), r"job(\d+)$" => s"\1") # string, but could be parsed as Int if needed
inputdir = "SiPhyNetwork_output/"
netfile_regex = r"sim(\d+)\.tree"
files =  filter!(x -> occursin(netfile_regex, x), readdir(inputdir))
N = length(files)

outfile = "quartets.csv"

## first: find anomalies by simulating gene trees
# @info "job $jobid: starting simulations to estimate CFs, loop over $N networks"

Random.seed!(seed)
dfts = repeat([DataFrame()], N) # 1 data frame per network, later concatenated into a single data frame

function analyzeTreeFile(treefile::String)

	m = match(netfile_regex, treefile)
	sim_num = parse(Int16, m.captures[1])
	tree = readTopology(joinpath(inputdir, treefile))
	deleteaboveLSA!(tree)
	# randomly sample 1 tip to prune, if requested by user
	if delete1==1
		taxa = tipLabels(tree)
		pruneit = taxa[rand(1:length(taxa))]
		deleteleaf!(tree, pruneit, simplify=false, nofuse=false)
	end

	taxa = sort(tipLabels(tree))
	numTaxa = length(taxa)
	numTaxa < 4 && error("fewer than 4 taxa, no quartets exist. simulated networks should have 4+ taxa!")

	# loop over quartets
	quartets = collect(combinations(1:numTaxa,4))
	nquartets = length(quartets)
	dft = DataFrame(
		sim_num = repeat([sim_num], nquartets),
		quartet_num  = repeat([""], nquartets),
		num3blob_col = repeat([-1], nquartets),
		num4blob_col = repeat([-1], nquartets),
		nsplit       = repeat([-1], nquartets),
		qCF_n        = repeat([-1], nquartets),
		split1       = repeat([Float64(-1)], nquartets),
		split2       = repeat([Float64(-1)], nquartets),
		split3       = repeat([Float64(-1)], nquartets),
		is32blob     = zeros(Bool, nquartets), # is there a 3_2 blob?
		flag_class   = zeros(Bool, nquartets), # is the class be perhaps underestimated?
		ngenes = zeros(Int, nquartets),
		is_anomalous = Vector{Union{Missing, Symbol}}(missing, nquartets),
		s1_from      = Vector{Symbol}(undef, nquartets),
		s2_from      = Vector{Symbol}(undef, nquartets),
		s3_from      = Vector{Symbol}(undef, nquartets),
	)
	Threads.@threads for j = 1:nquartets
		dft[j,2] = string(quartets[j])
		dft[j,3:end] = analyzeQuartet(quartets[j], taxa, tree; seed=seed)[1,:] # each quartet returns a DataFrame with one row and no sim_num
	end
	return(dft)
end

function analyzeQuartet(quartet, taxa, tree; seed=nothing)

	dfq = DataFrame(
		num3blob_col = -1,
		num4blob_col = -1,
		nsplit       = -1,
		qCF_n        = -1,
		split1       = Float64(-1),
		split2       = Float64(-1),
		split3       = Float64(-1),
		is32blob   = false,
		flag_class = false,
		ngenes = 0,
		is_anomalous = Vector{Union{Missing, Symbol}}(missing,1),
		s1_from = :undef,
		s2_from = :undef,
		s3_from = :undef,
		)

	quartet_taxa = taxa[quartet]
	notquartet_taxa = setdiff(taxa, quartet_taxa)
	quartettree = deepcopy(tree)
	for pruneit in notquartet_taxa
		deleteleaf!(quartettree, pruneit, simplify=false, nofuse=false)
	end
	deleteaboveLSA!(quartettree)
	
	quartet_blob_degree = blob_degree(quartettree)

	quartet_num3blob = 0
	quartet_num4blob = 0
	
	for degree in quartet_blob_degree[2]
		if degree == 3
			quartet_num3blob += 1
		elseif degree == 4
			quartet_num4blob += 1
		end
	end

	dfq[1,"num3blob_col"] = quartet_num3blob
	dfq[1,"num4blob_col"] = quartet_num4blob

	ns, qCF, _, tmpdf, is32blob, isanomalous, flag, o = quartettype_qCF(quartettree, ngt, rho;
	    verbose=false, blob_degrees=quartet_blob_degree, seed=seed)
	# using default attempt_resolve_ambiguity: false

	dfq[1,"nsplit"] = ns
	dfq[1,"qCF_n"] = ngt
	dfq[1,"split1"] = qCF[1]
	dfq[1,"split2"] = qCF[2]
	dfq[1,"split3"] = qCF[3]
	dfq[1,:is32blob]   = is32blob
	dfq[1,:flag_class] = flag
	dfq[1,:ngenes] = Int(tmpdf[1,:ngenes])
	dfq[1,:is_anomalous] = isanomalous
	dfq[1,:s1_from] = o[1]
	dfq[1,:s2_from] = o[2]
	dfq[1,:s3_from] = o[3]

	return(dfq)
end

for i = 1:N
	dfts[i] = analyzeTreeFile(files[i])
end
bigdf = vcat(dfts...)
# CSV.write(outfile, bigdf)

## second: find anomalies by calculating exact quartet CFs
# @info "job $jobid: starting calculation of CFs, loop over $N networks"

using QuartetNetworkGoodnessFit
const PN = PhyloNetworks

function quartetT_as_df(quartets::Vector{PN.QuartetT{T}}, sim_num) where T <: AbstractVector
  V = eltype(T)
  V <: Real || error("expected real data values")
  fourtax(q) = string(q.taxonnumber) # string(indexorder[q.taxonnumber])
  df = DataFrame(sim_num=Int[], quartet_num=String[],
            CF12_34=V[],CF13_24=V[],CF14_23=V[])
  for q in quartets
    push!(df, (sim_num, fourtax(q), q.data...) )
  end
  return df
end

dfts = repeat([DataFrame()], N) # 1 data frame per network, later concatenated

# loop over simulated networks
Threads.@threads for i in 1:N
  netfile = files[i]
  m = match(netfile_regex, netfile)
  sim_num = parse(Int16, m.captures[1])
  net = readTopology(joinpath(inputdir, netfile))
  q,t = network_expectedCF(net, showprogressbar=false, inheritancecorrelation=rho) # loops over 4-taxon sets
  t == sort(tipLabels(net)) || error("job $jobid, $netfile: sorted taxa and t don't match. taxa=$(taxa) and t=$t")
  dfts[i] = quartetT_as_df(q,sim_num)
end
# combine data frames across all networks
bigdf_exp = vcat(dfts...)
# CSV.write(outfile, bigdf_exp)
ntotal = nrow(bigdf)
nrow(bigdf_exp) == ntotal || error("simulated and exact CF tables with different # of rows")

# use simulation df to find which of the 3 CFs are split1,2,3
df_both = outerjoin(bigdf, bigdf_exp; on=[:sim_num,:quartet_num])
nrow(df_both) == ntotal || error("combined df with $(nrow(df_both)) rows instead of $ntotal")
function findsplit123(o123, cfe)
  names = [:CF12_34,:CF13_24,:CF14_23]
  nt = (; zip(names, cfe)...)
  return collect(nt[o123])
end
transform!(df_both,
  [:s1_from,:s2_from,:s3_from, :CF12_34,:CF13_24,:CF14_23] =>
  ByRow( (o1,o2,o3, c1,c2,c3) -> findsplit123((o1,o2,o3), (c1,c2,c3)) ) =>
  [:s1_exact,:s2_exact,:s3_exact])

CSV.write(outfile, df_both) # overwrites any pre-existing file of same name
