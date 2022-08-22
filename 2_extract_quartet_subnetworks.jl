print("2_extract_quartet_subnetworks.jl" * "\n")
using PhyloNetworks # master version
using Random
using QuartetNetworkGoodnessFit
Random.seed!(1718) # current time
cd("/media/john/Phylo/research/2022-05-18 six-fingered hand/six-fingered-hand")
include("3_find_blobs_and_degree.jl")
include("find_3blobqCF.jl")
inputdir = "SiPhyNetwork_output/"
#need to loop over trees
files = readdir(inputdir)
files = filter(x -> occursin(r"sim\d+\.tree", x), files)

global numquartet = 0
global num2cycle = 0
global num2blob = 0
global num3blob = 0
global num4blob = 0
global num5blob = 0
global num6blob = 0

global sim_num = Int16[]
global quartet_num = String[]
global num2cycle_col = Int16[]
global num2blob_col = Int16[]
global num3blob_col = Int16[]
global num4blob_col = Int16[]
global num5blob_col = Int16[]
global num6blob_col = Int16[]
global qCF_n = Int16[]
global split1 = Float16[]
global split2 = Float16[]
global split3 = Float16[]

global readTopologySuccess = 0
global readTopologyFail = 0

global readTopologyFailures = String[]

for file in files
	print(file * '\n')

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
		#print(m.captures[1] * '\n')
		#print(string(parse(Int16, m.captures[1])) * '\n')
		push!(sim_num, parse(Int16, m.captures[1]))
		global numquartet += 1
		push!(quartet_num, string(quartet))
		#print(" " * string(quartet) * "\n")

		quartet_taxa = taxa[quartet]
		notquartet_taxa = setdiff(taxa, quartet_taxa)
		quartettree = deepcopy(tree)
		for pruneit in notquartet_taxa
			deleteleaf!(quartettree, pruneit, simplify=false, nofuse=false)
		end
		deleteaboveLSA!(quartettree)
		
		quartet_blob_degree = blob_degree(quartettree)

		global quartet_num2cycle = 0
		global quartet_num2blob = 0
		global quartet_num3blob = 0
		global quartet_num4blob = 0
		global quartet_num5blob = 0
		global quartet_num6blob = 0

		for blob in biconnectedComponents(quartettree)
			if length(blob)==2
				quartet_num2cycle += 1
				global num2cycle += 1
			end
		end
		
		for degree in quartet_blob_degree[2]
			if degree == 2
				quartet_num2blob += 1
				global num2blob += 1
			elseif degree == 3
				quartet_num3blob += 1
				global num3blob += 1
			elseif degree == 4
				quartet_num4blob += 1
				global num4blob += 1
			elseif degree == 5
				quartet_num5blob += 1
				global num5blob += 1
			elseif degree == 6
				quartet_num6blob += 1
				global num6blob += 1
			end
		end

		push!(num2cycle_col, quartet_num2cycle)
		push!(num2blob_col, quartet_num2blob)
		push!(num3blob_col, quartet_num3blob)
		push!(num4blob_col, quartet_num4blob)
		push!(num5blob_col, quartet_num5blob)
		push!(num6blob_col, quartet_num6blob)
		
		#print(" quartet " * string(quartet) * "\n")
		#print("  quartet_num4blob=" * string(quartet_num4blob) * "\n")
		#print("  quartet_num3blob=" * string(quartet_num3blob) * "\n")
		
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
			ngenes=3200
			outputdir = "hl_output";
			ispath(outputdir) || mkdir(outputdir);
			gt = joinpath(outputdir, "d3blob"); # this file will be created then deleted
			try
				ns, qCF, hwc, df = quartettype_qCF(quartettree, gt, ngenes; seed=321, verbose=false)
				#print("  " * string(qCF) * "\n")
			catch y
				#print("  encountered error in quartettype_qCF" * "\n")
				#print(string(y) * "\n")
			end
		end

		push!(qCF_n, ngenes)
		push!(split1, qCF[1])
		push!(split2, qCF[2])
		push!(split3, qCF[3])

		if quartet == quartets[1]
			global subnetworks = [quartettree]
		else
			push!(subnetworks, quartettree)
		end
	end
	#print("numquartet " * string(numquartet) * "\n")
	#print("num4blob=" * string(num4blob) * "\n")
	#print("num3blob=" * string(num3blob) * "\n")
	if file == files[1]
		global sim_subnetworks = [subnetworks]
	else
		push!(sim_subnetworks, subnetworks)
	end
end

using DataFrames
df = DataFrame(
	sim_num=sim_num,
	quartet_num=quartet_num,
	num2cycle_col=num2cycle_col,
	num2blob_col=num2blob_col,
	num3blob_col=num3blob_col,
	num4blob_col=num4blob_col,
	num5blob_col=num5blob_col,
	num6blob_col=num6blob_col,
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




