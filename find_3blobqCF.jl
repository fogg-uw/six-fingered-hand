#= code to determine if:

1. a 4-taxon set is a pathological T-quartet, in the sense that
   it has a cut edge ab|cd but the associated quartet CF
   is the smallest: CF(ab|cd) < CF(ac|bd) = CF(ad|bc)

2. a 4-taxon set with a circular ordering is a pathological B-quartet,
   in the sense that it displays 2 splits, say ab|cd and ac|bc,
   but the third split ad|bc has a quartet CF that is not the smallest:
   CF(ad|bc) > either CF(ab|cd) or CF(ac|bd).

=#

using Random
using PhyloNetworks # ] add PhyloNetworks#master on 2020-05-20
using PhyloPlots
using QuartetNetworkGoodnessFit # to install & use HybridLambda
#hybridlambda = QuartetNetworkGoodnessFit.hybridlambda # path to hybrid-lambda simulator, on local machine
hybridlambda = "/media/john/Phylo/research/2022-05-18 six-fingered hand/hybrid-Lambda-v0.6.2-beta-exec" # path to hybrid-lambda simulator, on local machine
using HypothesisTests
pcspath = "/media/john/Phylo/research/2022-06-15 PhyloCoalSimulations/PhyloCoalSimulations.jl/src/PhyloCoalSimulations.jl" # PhyloCoalSimulations, on local machine
include(pcspath) 

"""
    quartettype_qCF(net, genetreefile, nsim=200; seed=nothing, verbose=true)

Calculate the quartet type of quartet concordance factors (CFs) of a
4-taxon network. The quartet type is described by a matrix with 1 row
per split in the displayed tree, and 1 column per taxon (so 4 columns).
If the network has a non-trivial cut-edge (and therefore no 4-degree blob),
then this matrix has a single row, and describes the split of the cut-edge.

Quartet CFs are estimated by simulating gene trees under the coalescent model
using HybridLambda. Warning: HybridLambda seems to assume that the network is
ultrametric.

output:
- number of splits in the displayed trees.
  If there is only 1 split, then the network contains a non-trivial cut-edge
  (based on proof in "tree of blobs" Allman et al. 2022 preprint).
  If there are 2 or more splits, then the semidirected network contains a 4-blob.
  If there are 3 splits, then the network is not outerplanar.
- estimated quartet CFs, in the following order:
  split1, split2 (if relevant), split3
  where split1 is the first split in the list of splits in the displayed trees,
  split2 is the second split, if there is more than 1, and
  split3 is the last split not already listed.
- matrix of splits: of size #splits × 4, with the 4 taxa listed in the same
  order as in the data frame
- data frame with a single row, containing the 4 taxon labels and
  the 3 estimated quartet CFs.

# examples

```julia
julia> include("find_3blobqCF.jl")

julia> net_3blob = readTopology("((B:0.6,((A:0.4,C:0.4):0.1)#H1:0.1::0.51):1.0,(#H1:0.1::0.49,O:0.6):1.0);");

julia> ngenes = 10_000 # number of genes to be simulated

julia> outputdir = "hl_output";

julia> ispath(outputdir) || mkdir(outputdir);

julia> gt = joinpath(outputdir, "d3blob"); # this file will be created then deleted

julia> ns, qCF, hwc, df = quartettype_qCF(net_3blob, gt, ngenes; seed=321, verbose=true);
Default Kingman coalescent on all branches.
Default population size of 10000 on all branches. 
Random seed: 321
Produced gene tree files: 
hl_output/d3blob_coal_unit
Reading in trees, looking at 1 quartets in each...
0+--------------------------------------------------+100%
  **************************************************

julia> ns # only 1 split: it's AC|BD based on hwc below
1
  
julia> df
1×8 DataFrame
 Row │ t1      t2      t3      t4      CF12_34  CF13_24  CF14_23  ngenes  
     │ String  String  String  String  Float64  Float64  Float64  Float64 
─────┼────────────────────────────────────────────────────────────────────
   1 │ A       B       C       O        0.3521   0.2904   0.3575  10000.0

julia> hwc # hardwired clusters, 4 taxa in same order as t1-t4 in data frame above
1×4 BitMatrix:
 1  0  1  0

julia> qCF # order corresponding to 1st split in hwc matrix
(major = 0.2904, alt1 = 0.3521, alt2 = 0.3575)

julia> n_major = Int(round(qCF[:split1] * ngenes, digits=8));

julia> n_alt1  = Int(round(qCF[:split2] * ngenes, digits=8));

julia> n_alt2  = Int(round(qCF[:split3] * ngenes, digits=8));

julia> bt = BinomialTest(n_alt1, n_alt1 + n_alt2); # equal alternative CFs?

julia> pvalue(bt) # the data are consistent with equal alternative CFs.
0.5292396850392722

julia> bt = BinomialTest(n_major, ngenes, 1/3); # is CF_major < 1/3?

julia> pvalue(bt, tail=:left) # reject: so this is a pathological network
3.492803208145045e-20
```

Here is another example where the network has a 4-degree blob, and is outerplanar.
This case is not pathological because the first 2 quartet CFs, which correspond
to the 2 splits in the displayed trees, indeed have the 2 largest CFs,
so they tell us the correct circular ordering.

```julia
julia> net_4blob = readTopology("((((T:0.5)#H1:0.6::0.51,C:1.1):0.7,(#H1:0.0::0.49,E:0.5):1.3):1.0,O:2.8);");

julia> gt = joinpath(outputdir, "d4blob");

julia> ns, qCF, hwc, df = quartettype_qCF(net_4blob, gt, ngenes; seed=321, verbose=true);

julia> ns # 2 splits in the displayed trees: CT and ET based on hwc and df below
2

julia> hwc
2×4 BitMatrix:
 1  0  0  1
 0  1  0  1

julia> df
1×8 DataFrame
 Row │ t1      t2      t3      t4      CF12_34  CF13_24  CF14_23  ngenes  
     │ String  String  String  String  Float64  Float64  Float64  Float64 
─────┼────────────────────────────────────────────────────────────────────
   1 │ C       E       O       T        0.1276   0.4831   0.3893  10000.0

julia> qCF
(split1 = 0.3893, split2 = 0.4831, split3 = 0.1276)

julia> n_cf1    = Int(round(qCF[:split1] * ngenes, digits=8));

julia> n_cf2    = Int(round(qCF[:split2] * ngenes, digits=8));

julia> n_minor  = Int(round(qCF[:split3] * ngenes, digits=8));

julia> bt = BinomialTest(n_cf1, n_cf1 + n_minor); # is CF1 > CF_minor?

julia> pvalue(bt, tail=:left) # the data are consistent with CF1 >= CF_minor
1.0

julia> pvalue(BinomialTest(n_cf2, n_cf2 + n_minor), tail=:left) # also with CF2 >= CF_minor
1.0
```
"""
function quartettype_qCF(net::HybridNetwork, genetreefile::AbstractString,
        nsim=200; seed=nothing, verbose=true)

    """
    # ulrametrize to the degree of precision needed by hybridlambda

    for edge in net.edge
      edge.length = floor(edge.length, digits=10)
    end
    """

    """
    # label the quartet tips q1, q2, q3, q4 to avoid sorting problems later

    for leaf in net.leaf

    end
    """

    taxonlist = sort(tipLabels(net))
    length(taxonlist) == 4 || error("there aren't 4 tips: $taxonlist")
    # collect splits appearing in any displayed tree
    dtree = displayedTrees(net, 0.0)
    mat = BitMatrix(undef, (0,4)) # initialize: 1 row per split
    for tree in dtree
        mat = vcat(mat, Bool.(hardwiredClusters(tree, taxonlist)[1:end,2:(end-1)]))
    end
    # keep non-trivial splits that partition the 4 taxa in 2 vs 2
    mat = mat[ sum(mat, dims=2)[:,1] .== 2 ,:]
    # remove duplicate rows, from nodes with the same hardwired clusters
    tokeep = Bool.(ones(size(mat,1)))
    for i in 1:size(mat,1)
        for j in 1:(i-1) # consider unrooted splits
            if mat[i,:] == mat[j,:] || mat[i,:] == .!mat[j,:]
                tokeep[i] = false
                break
            end
        end
    end
    mat = mat[tokeep,:]
    nsplits = size(mat,1)
    isnothing(seed) || Random.seed!(seed)

    #gtcu = genetreefile * "_coal_unit"  # name of output file created by hybrid-Lambda
    #netHL = hybridlambdaformat(net) # string format for the network, as required by hybrid-lambda
    
    # hack-y fix to hybridlamdaformat(); need to change hybrid node notation from (e.g.) #H22:1.0 to H22#:1.0
    # cecile is working on a (less hack-y?) fix to the hybridlambdaformat() function.
    #re = r"(\#)(H\d+)(\:-?\d+)" # e.g #H22:1.0 or #H22:-0.1
    #su = s"\2\1\3" # move # in position 1 to position 2
    #netHL = replace(netHL, re => su)
    
    #netHL = "'$netHL'" # quoted
    #hlcommand = `$hybridlambda -spcu $netHL -num $nsim -seed $seed -o $genetreefile`
    #hlout = ( verbose ? stdout : devnull )
    #run(pipeline(hlcommand; stderr = hlout));
    #treelist = readMultiTopology(gtcu)
    #rm(gtcu)

    treelist = PhyloCoalSimulations.simulatecoalescent(net, nsim, 1)

    obsCF, t = countquartetsintrees(treelist; showprogressbar=verbose)

    # hybridlambda adds suffix "_1" to each taxon name: remove it
    #map!(name -> replace(name, r"_1$" => s""), t, t)

    taxonlist == t || @error("different order of taxa used by countquartetsintrees (but maybe john's code will help)")
    splittype(s) = (s==[1,1,0,0] || s == [0,0,1,1] ? 0 :
                    (s== [1,0,1,0] || s == [0,1,0,1] ? 1 : 2))
    splitsymbol(st) = (st==0 ? :CF12_34 : (st==1 ? :CF13_24 : :CF14_23))

    if taxonlist != t

      # find how countquartetsintrees permuted the taxa
      perm_t = Int16[]
      for taxon in taxonlist
        push!(perm_t, findall(x->x==taxon, t)[1])
      end

      # taxon order (and therefore split order) has been "corrupted" by countquartetsintrees
      corrupted_CF =
        [obsCF[1].data[1],
         obsCF[1].data[2],
         obsCF[1].data[3]]

      # we must "reform" it by un-permuting the data
      perm_splits = BitMatrix(undef, (0,4))
      pat1 = [1,1,0,0]
      pat2 = [1,0,1,0]
      pat3 = [1,0,0,1]
      perm_splits = vcat(perm_splits, transpose([pat1[perm_t[1]], pat1[perm_t[2]], pat1[perm_t[3]], pat1[perm_t[4]]]))
      perm_splits = vcat(perm_splits, transpose([pat2[perm_t[1]], pat2[perm_t[2]], pat2[perm_t[3]], pat2[perm_t[4]]]))
      perm_splits = vcat(perm_splits, transpose([pat3[perm_t[1]], pat3[perm_t[2]], pat3[perm_t[3]], pat3[perm_t[4]]]))

      perm_splits_2 =
      [splittype(perm_splits[1,:]) + 1,
       splittype(perm_splits[2,:]) + 1,
       splittype(perm_splits[3,:]) + 1]

      obsCF[1].data[perm_splits_2[1]] = corrupted_CF[1]
      obsCF[1].data[perm_splits_2[2]] = corrupted_CF[2]
      obsCF[1].data[perm_splits_2[3]] = corrupted_CF[3]

      t = taxonlist

    end

    df = writeTableCF(obsCF, t)

    st1 = splittype(mat[1,:])
    st2 = (nsplits > 1 ? splittype(mat[2,:]) : mod(st1+1,3))
    st3 = setdiff(0:2, [st1,st2])[1]
    o = splitsymbol.([st1,st2,st3])
    return nsplits, (split1=df[1,o[1]], split2=df[1,o[2]], split3=df[1,o[3]]), mat, df
end

"""
    net6: example network

example with 6 taxa and 1 reticulation, from which we can extract various
4-taxon network cases, with a 2-degree or 3-degree or 4-degree blob.

```julia
net6 = readTopology("((((D:0.2,C:0.2):2.4,((A:0.4,B:0.4):1.1)#H1:1.1::0.51):2.0,(#H1:0.0::0.49,E:1.5):3.1):1.0,O:5.6);");
```

how to exact a 4-taxon network with a 3-blob, and change it edge lengths
to make it pathological --but still make it ultrametric:

```julia
net_3blob = deepcopy(net6);
for tip in ["E","D"] deleteleaf!(net_3blob, tip); end
net_3blob.edge[1].length = 0.6 # external to C
net_3blob.edge[9].length = -0.4 # external to 0
net_3blob.edge[4].length = 0.1 # edge number 6 below the hybrid node
net_3blob.edge[5].length = 0.1 # edge number 7: major hybrid edge
net_3blob.edge[7].length = 0.1 # edge number 9: minor hybrid edge
rootonedge!(net_3blob,8)
plot(net_3blob, :R, useEdgeLength=true, showEdgeNumber=true, showEdgeLength=true);
QuartetNetworkGoodnessFit.ultrametrize!(net_3blob, true)
# output is true and no warning: it is now ultrametric
writeTopology(net_3blob, round=true)

net_4blob = deepcopy(net6);
for tip in ["B","D"] deleteleaf!(net_4blob, tip); end
"""


"""
Failed attempt to find the number of blobs of degrees 2,3,4 in a 4-taxon network,
by calculating the total number "n" of blobs (including trivial blobs, which are
single cut edges), and by calculating the number of non-trivial blobs
(blobs with 3 or more nodes).
If we denote ni = number of i-debree blobs for i=2,3,4, then we must have:

number of non-trivial blobs = n2+n3+n4

and also, with n = total number of blobs (including trivial ones):

n = n4 + 2n2 + 4 if there's a 4-degree blob
n = n3 + 2n2 + 5 if there's 1 or 2 3-degree blob

in all cases: n = n3 + 2n2 + 5.

problem: using the number of trivial / non-trivial blobs is not enough.
example: n4=1, n2=1 is indistinguishable from n3=2, n2=0.
both have 2 non-trivials blobs, and 7 blobs total.

bi = PhyloNetworks.blobInfo(net, false)[2]
# false: to include trivial blobs. [2]: to get the lists of major hybrid edges
nblobs = sum(.!isempty.(bi)) # trivial blobs have an empty list of hybrid edges
# nblobs = number of non-trivial blobs

n3p2n2 = length(bi) - 5 # n-5
n2,n3 = divrem(n3p2n2,2) # correct if n3 = 0 or 1, incorrect if n3=2
if n3 == 1
    n4 = 0
    n2+n3 == nblobs || error("blob counting has a problem")
else # n3=0 truly or n3 should be 2
    # then we cannot tell for sure.
end

bcc = PhyloNetworks.biconnectedComponents(net_4blob, true)
"""

nothing # to avoid screen output when we include this file
