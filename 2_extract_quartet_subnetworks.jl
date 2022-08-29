# usage: julia 2_extract_quartet_subnetworks.jl [ngt]
# ngt is the number of gene trees for PCS to simulate per quartet network.

print("2_extract_quartet_subnetworks.jl" * "\n")
#print(pwd())
using PhyloNetworks # master version
#cd("/media/john/Phylo/research/2022-05-18 six-fingered hand/vary_params/six-fingered-hand")
include("find_blobs_and_degree.jl")
include("find_3blobqCF.jl")
ngt = parse(Int64, ARGS[1])
inputdir = "SiPhyNetwork_output/"
#need to loop over trees
files = readdir(inputdir)
files = filter(x -> occursin(r"sim\d+\.tree", x), files)

#scalar variables (updated as we loop through trees; mostly for debugging)

global numquartet = 0
global num3blob = 0
global num4blob = 0
global readTopologySuccess = 0
global readTopologyFail = 0

#vector variables (record new entries as we loop through trees/networks; export at end)

global sim_num = Int16[]
global quartet_num = String[]
global num3blob_col = Int16[]
global num4blob_col = Int16[]
global qCF_n = Int16[]
global split1 = Float16[]
global split2 = Float16[]
global split3 = Float16[]

# this last is for debugging.

global readTopologyFailures = String[]

for file in files
	#print(file * '\n')

	try
		global tree = readTopology(joinpath(inputdir, file))
	catch
		print("couldn't read topology")
		global readTopologyFail += 1
		push!(readTopologyFailures, file)
		continue
	end
	global readTopologySuccess += 1

	taxa = tipLabels(tree)
	numTaxa = length(taxa)

	# loop over quartets
	using Combinatorics
	if numTaxa < 4
		print("fewer than 4 taxa, no quartets exist\n")
		continue
	end
	
	quartets = collect(combinations(1:numTaxa,4))
	for quartet in quartets
		m = match(r"sim(\d+)\.tree", file)
		push!(sim_num, parse(Int16, m.captures[1]))
		global numquartet += 1
		push!(quartet_num, string(quartet))

		quartet_taxa = taxa[quartet]
		notquartet_taxa = setdiff(taxa, quartet_taxa)
		quartettree = deepcopy(tree)
		for pruneit in notquartet_taxa
			deleteleaf!(quartettree, pruneit, simplify=false, nofuse=false)
		end
		deleteaboveLSA!(quartettree)
		
		quartet_blob_degree = blob_degree(quartettree)

		global quartet_num4blob = 0
		global quartet_num3blob = 0
		for degree in quartet_blob_degree[2]
			if degree == 3
				quartet_num3blob += 1
				global num3blob += 1
			end
			if degree == 4
				quartet_num4blob += 1
				global num4blob += 1
			end
		end

		push!(num3blob_col, quartet_num3blob)
		push!(num4blob_col, quartet_num4blob)
		
		ngenes = -1
		qCF = [-1,-1,-1]

		if quartet_num4blob > 0
			#print("   4. if there is a 4-degree blob: all good. B-quartet, no pathology." * "\n")
			nothing
		elseif quartet_num4blob == 0 && quartet_num3blob == 0
			#print("   5. if there is no 4-degree blob and no 3-degree blob: probably all good (?)." * "\n")
			nothing
		elseif quartet_num3blob > 0
			#print("   6. if there is a 3-degree blob (therefore no 4-degree blob), then:")
			ngenes=ngt
			try
				ns, qCF, hwc, df = quartettype_qCF(quartettree, ngenes; seed=321, verbose=false)
			catch y
				print("  encountered error in quartettype_qCF" * "\n")
				print(string(y) * "\n")
			end
		end

		push!(qCF_n, ngenes)
		push!(split1, qCF[1])
		push!(split2, qCF[2])
		push!(split3, qCF[3])
	end
end

# remove SiPhyNetwork_output
files = readdir(inputdir)
for file in files
	rm(joinpath(inputdir, file))
end
rm(inputdir, recursive=true)

using DataFrames
df = DataFrame(
	sim_num=sim_num,
	quartet_num=quartet_num,
	num3blob_col=num3blob_col,
	num4blob_col=num4blob_col,
	qCF_n=qCF_n,
	split1=split1,
	split2=split2,
	split3=split3
	)
using CSV
try
	rm("quartets.csv")
catch
	nothing
end
CSV.write("quartets.csv", df)




