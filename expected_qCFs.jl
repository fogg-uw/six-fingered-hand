# obsolete! folded into step 2. kept for code to re-order taxa, lines ~ 76-79

#= Calculate the expected CFs of all quartets, for all simulated networks
   across all "jobs", where 1 job = 1 combination of parameter choices

   assumes: independent inheritance (rho=0) --so far at least

usage: navigate to the directory parent to all "job$i" with i = index of scenario
where the input files should exist.
input files: job$i/SiPhyNetwork_output/sim*.tree
output file: job$i/quartets_expCF.csv

usage: julia ../expected_qCFs.jl i_start i_end

options
i_start, i_end = start and end IDs of the job folders to treat.

rows in "quartets.csv" contain this:
sim_num,quartet_num,num3blob_col,num4blob_col,nsplit,qCF_n,split1,split2,split3,is32blob,flag_class,ngenes,is_anomalous
1,"[1, 2, 3, 4]",0,0,1,300,0.53,0.25666666666666665,0.21333333333333335,false,false,300,good
where:
[1, 2, 3, 4] correspond to taxon indices in tipLabels(net) -- unsorted
split1 = major split, if only 1 split displayed in network
split1, split2 = splits displayed in network, if only 2 of them
=#

using PhyloNetworks
using QuartetNetworkGoodnessFit
using CSV, DataFrames
using Combinatorics
const PN = PhyloNetworks

function quartetT_as_df(quartets::Vector{PN.QuartetT{T}}, indexorder, sim_num) where T <: AbstractVector
  V = eltype(T)
  V <: Real || error("expected real data values")
  fourtax(q) = string(indexorder[q.taxonnumber])
  df = DataFrame(sim_num=Int[], quartet_num=String[],
            f12_34=V[],f13_24=V[],f14_23=V[])
  for q in quartets
    push!(df, (sim_num, fourtax(q), q.data...) )
  end
  return df
end
  
# read job_{start,end} arguments and check that "job" folders exits
jobid_regex = r"^job(\d+)$"
jobid_seen = sort(parse.(Int, replace.(filter(x -> occursin(jobid_regex,x), readdir()), jobid_regex => s"\1")))
isempty(jobid_seen) && error("no folder named 'jobxxx'")
jobid_start = (length(ARGS) > 0 ? parse(Int, ARGS[1]) : jobid_seen[1])
if jobid_start < jobid_seen[1]
    @error "starting job $jobid_start too small: changing to $(jobid_seen[1])"
    jobid_start = jobid_seen[1]
end
jobid_end   = (length(ARGS) > 1 ? parse(Int, ARGS[2]) : jobid_seen[end])
if jobid_end > jobid_seen[end]
    @error "ending job $jobid_end too large: changing to $(jobid_seen[end])"
    jobid_end = jobid_seen[end]
end

netfile_regex = r"sim(\d+)\.tree"

for jobid in jobid_start:jobid_end
  jobdir = "job$jobid"
  inputdir = joinpath(jobdir, "SiPhyNetwork_output")
  files = readdir(inputdir)
  files = filter(x -> occursin(netfile_regex, x), files)
  outfile = joinpath(jobdir, "quartets_expCF.csv")
  @info "starting job $jobid, found $(length(files)) simulated networks"
  dfts = DataFrame[] # there will be one for each network, later concatenated
  # loop over simulated networks
  Threads.@threads for netfile in files
    m = match(netfile_regex, netfile)
    sim_num = parse(Int16, m.captures[1])
    net = readTopology(joinpath(inputdir, netfile))
    q,t = network_expectedCF(net) # loops over 4-taxon sets
    taxa = tipLabels(net) # compare to t, then get o such that taxa[o[i]] = t[i]
    o = sortperm(taxa) # for taxon number that match between quartets.csv and quartets_expCF.csv
    taxa[o] == t || error("job $jobid, $netfile: sorted taxa and t don't match. taxa=$(taxa) and t=$t")
    push!(dfts, quartetT_as_df(q,o,sim_num))
  end # of loop over networks
  # combine data frames across all networks
  bigdf = vcat(dfts...)
  CSV.write(outfile, bigdf)
end # of loop over jobs
