using CSV, DataFrames
using AlgebraOfGraphics, CairoMakie
set_aog_theme!()

#= need result files stored in /nobackup2/
done on franklin01:
cd /nobackup2/ane/simulation-pathologyquartets/results_base/
=#

## check agreement btw qCFs exact and estimated with simulation
#  using job117 only, but all 800 networks

jobid = 117
jobdir = joinpath("results_base", "job$jobid")
netdir = joinpath(jobdir, "SiPhyNetwork_output")
qfile = joinpath(jobdir, "quartets.csv")

##

qdat = DataFrame(CSV.File(qfile))
select!(qdat, r"^[sqni]")
filter!(:nsplit => n -> n<3, qdat) # CFs not estimated if 3 splits
@show nrow(qdat) # 52585 four-taxon sets

tck = [0,1/3,0.5,1]; tcklab = ["0","1/3","0.5","1"]

plt1 = data(qdat) * mapping(
  :split1,
  :s1_exact => "exact CF";
  row = :ngenes => nonnumeric,
) * visual(; alpha=0.2, markersize=2);
plt2 = data(qdat) *
  mapping(:split2, :s2_exact, row = :ngenes => nonnumeric) *
  visual(; alpha=0.2, color=:purple, markersize=2);
plt3 = data(qdat) *
  mapping(:split3, :s3_exact, row = :ngenes => nonnumeric) *
  visual(; alpha=0.2, color=:green, markersize=2);
plt = plt1 + plt2 + plt3;
fig = draw(plt;
 axis=(width=200, height=200, aspect=1,
       xticks=(tck,tcklab), yticks=(tck,tcklab)),
)
save("qCF_exactvssimulated.pdf", fig)

## job117, network 607: why no 3_2 blob yet CF1 < 1/3 ?

using PhyloNetworks
using RCall, PhyloPlots
include("find_anomalies.jl")

## read offending network and 4-taxon subnet

simid = "607"
net = readTopology(joinpath(netdir, "sim$simid.tree")) # 8 tips, h=5
taxa = sort(tipLabels(net))
q1i = [1, 3, 5, 7]; q1 = taxa[q1i] # l10 l8 t15 t5
q2i = [3, 4, 5, 7]; q2 = taxa[q2i] #  l8 l9 t15 t5

n1 = deepcopy(net)
for t in taxa  t in q1 || deleteleaf!(n1, t); end
# (((t5:0.361,#H21:0.0::0.496):0.51,#H12:0.0::0.296):0.126,(((((((l8:0.282,l10:0.282):0.064,#H24:0.0::0.219):0.015)#H21:0.056::0.504)#H18:0.0::0.619,(t15:0.346)#H24:0.071::0.781):0.012,#H18:0.012::0.381):0.442)#H12:0.126::0.704);

n2 = deepcopy(net)
for t in taxa  t in q2 || deleteleaf!(n2, t); end

## plot full network and quarnets

R"par"(mar=[0,0,0,0]); R"layout"([1 2 3]);
plot(net, useedgelength=true, showedgenumber=true);
plot(n1, useedgelength=false);
plot(n2, useedgelength=false);

bcc, bdegree = blob_degree(n1)
ns, qCF, hwc, df, is32blob, isanomalous, flag, o = quartettype_qCF(n1, 10, 0.0; verbose=true, blob_degrees=(bcc, bdegree))
is32blob # false!
blobexit = PhyloNetworks.biconnectedcomponent_exitnodes(n1, bcc, false) # correct

#= conclusions:
- how to best define 3_2 blobs?
- if: "3-blob with a sister pair below hybrid node",
  then there are anomalous networks with a 3-blob that's not a 3_2 blob.
=#
