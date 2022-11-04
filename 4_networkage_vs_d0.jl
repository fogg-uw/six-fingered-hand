#=
code to look at the age of the simulated networks, and compare that
to the distance threshold for d0 above which reticulation fails.

All simulated networks are ultrametric, because extinct taxa were pruned
using SiPhyNetwork's option `complete = FALSE`, and no rate variation was
modelled in any way. Ultrametricity makes it easier to get the network's age.

run from main folder here with `julia --project` on franklin01
=#

using PhyloNetworks
using Statistics
using CSV, DataFrames

# where all intermediate files were placed: on franklin01 at
netdir = "/nobackup2/ane/simulation-pathologyquartets/results20221103_10params/"

# extract network ages
networkage = Vector{Vector{Float64}}() # one vector for each job
njobs = sum(occursin.(r"^job\d+$", readdir(netdir)))
for jobid in 1:njobs
  agevec = Vector{Float64}()
  jobdir = joinpath(netdir, "job$jobid", "SiPhyNetwork_output")
  nsims = sum(occursin.(r"^sim\d+.tree$", readdir(jobdir)))
  for simid in 1:nsims
    netfile = joinpath(jobdir, "sim$(lpad(simid, 3, '0')).tree")
    net = readTopology(netfile)
    deleteaboveLSA!(net)
    push!(agevec, maximum(PhyloNetworks.getHeights(net)))
  end
  push!(networkage, agevec)
end

# extract d0 associated with each job
df = CSV.read(joinpath(netdir,"results_withtime.csv"), DataFrame)
d0relative = round.(df.d_0 .* df.lambda, digits=5)
unique(d0relative) # 3 values only: good.
# d0: threshold in coalescent units
# d0relative: threshold in average # of speciations
df_age = DataFrame(
    :job    => 1:njobs,
    :ntaxa  => df.ntaxa,
    :lambda => df.lambda,
    :d0_relative => d0relative,
    :d0    => df.d_0,
    :mu     => df.mu,
    :networkage_median => median.(networkage),
    :networkage_mean   => mean.(networkage))
CSV.write("networkage.csv", df_age)

# look at summaries
combine(groupby(df_age, [:lambda], sort=true),
    :networkage_median => mean,
    :networkage_mean => mean)
#=
 lambda   networkage_median_mean  networkage_mean_mean 
───────────────────────────────────────────────────────
    0.1               13.6532               21.079
    0.3                4.46732               7.83383
    1.0                1.36025               2.18843
    3.0                0.452177              0.717225
=#
tmp = combine(groupby(df_age, [:lambda, :d0], sort=true),
    :networkage_median => mean => :age_median,
    :networkage_mean   => mean => :age_mean )
tmp.ratio_d0_2age = tmp.d0 ./ (2 .* tmp.age_median)
tmp.ratio_2age_d0 = (2 .* tmp.age_median) ./ tmp.d0
tmp
#=
lambda   d0          age_median  age_mean   ratio_d0_2age  ratio_2age_d0 
─────────────────────────────────────────────────────────────────────────
    0.1   2.0         14.1384    18.7829        0.0707291       14.1384
    0.1   6.0         13.4934    21.3735        0.222331         4.49779
    0.1  12.0         13.3277    23.0807        0.450189         2.22129
    0.3   0.666667     4.70761    6.20979       0.0708073       14.1228
    0.3   2.0          4.42093    7.77029       0.226197         4.42093
    0.3   4.0          4.2734     9.52141       0.468011         2.1367
    1.0   0.2          1.38956    1.89177       0.0719654       13.8956
    1.0   0.6          1.36869    2.22927       0.219188         4.5623
    1.0   1.2          1.3225     2.44424       0.453685         2.20417
    3.0   0.0666667    0.461922   0.645232      0.0721622       13.8577
    3.0   0.2          0.459902   0.716559      0.217438         4.59902
    3.0   0.4          0.434706   0.789883      0.460082         2.17353
=#
