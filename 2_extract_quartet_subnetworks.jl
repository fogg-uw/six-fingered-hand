# usage: julia 2_extract_quartet_subnetworks.jl [ngt] [delete1]
# ngt is the number of gene trees for PCS to simulate per quartet network.
# delete1 = 1,0 is whether to delete 1 (one) leaf from each network before doing anything else with it (like examining quartets).

using PhyloNetworks
using Random
using DataFrames
using Combinatorics
using CSV

include("find_blobs_and_degree.jl")
include("find_3blobqCF.jl")

ngt = parse(Int64, ARGS[1])
delete1 = parse(Int8, ARGS[2])

if !(delete1 == 0 || delete1 == 1)
	print("delete1 not 0 or 1, interpreting as 0")
	delete1 = 0
end

seed = 9 # nine rings for men
Random.seed!(seed)

#need to loop over trees

inputdir = "SiPhyNetwork_output/"
files = readdir(inputdir)
files = filter(x -> occursin(r"sim\d+\.tree", x), files)
N = length(files)
dfts = repeat([DataFrame()], N) # N data frames, one for each tree, later compressed into one data frame

# analog to "zeros" function for strings

blankstrings = function(N=1)
	empty = String[]
	for i in 1:N
		push!(empty, "")
	end
	return(empty)
end

# "global" scalar variables (exist outside of trees; mostly for debugging)

readTopologySuccess = zeros(Int8, N)
readTopologyFailures = blankstrings(N)

function analyzeTreeFile(treefile::String, treenum::Int64)

	#print(file * '\n')

	m = match(r"sim(\d+)\.tree", treefile)
	sim_num = parse(Int16, m.captures[1])

	tree = PhyloNetworks.Network

	try
		tree = readTopology(joinpath(inputdir, treefile))
	catch
		print("couldn't read topology")
		readTopologyFailures[treenum] = treefile
		return(NA)
	end
	readTopologySuccess[treenum] = 1

	deleteaboveLSA!(tree)

	# randomly sample 1 tip to prune, if requested by user
	if delete1==1
		taxa = tipLabels(tree)
		pruneit = taxa[rand(1:length(taxa))]
		deleteleaf!(tree, pruneit, simplify=false, nofuse=false)
	end

	taxa = tipLabels(tree)
	numTaxa = length(taxa)

	# loop over quartets

	if numTaxa < 4
		print("fewer than 4 taxa, no quartets exist\n")
		return(NA)
	end

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
		split3       = repeat([Float64(-1)], nquartets)
	)

	for j in nquartets
		dft[j,2] = string(quartets[j])
		dft[j,3:9] = analyzeQuartet(quartets[j], taxa, tree)[1,1:7] # each quartet returns a DataFrame with one row and no sim_num
	end

	return(dft)

end

function analyzeQuartet(quartet, taxa, tree)

	dfq = DataFrame(
		num3blob_col = -1,
		num4blob_col = -1,
		nsplit       = -1,
		qCF_n        = -1,
		split1       = Float64(-1),
		split2       = Float64(-1),
		split3       = Float64(-1)
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

	ns, qCF, hwc, df = quartettype_qCF(quartettree, ngt; seed=321, verbose=false)

	dfq[1,"nsplit"] = ns
	dfq[1,"qCF_n"] = ngt
	dfq[1,"split1"] = qCF[1]
	dfq[1,"split2"] = qCF[2]
	dfq[1,"split3"] = qCF[3]

	return(dfq)

end

Threads.@threads for i = 1:N
	dfts[i] = analyzeTreeFile(files[i], i)
end

global bigdf = dfts[1]
for i in 2:N
	global bigdf = vcat(bigdf, dfts[i])
end

try
	rm("quartets.csv")
catch
	nothing
end
CSV.write("quartets.csv", bigdf)




