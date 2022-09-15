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
using PhyloNetworks
using PhyloCoalSimulations
using DataFrames

"""
    quartettype_qCF(net, nsim=200; seed=nothing, verbose=true,
                      threshold_h=nothing)

Calculate the quartet type of quartet concordance factors (CFs) of a
4-taxon network. The quartet type is described by a matrix with 1 row
per split in the displayed tree, and 1 column per taxon (so 4 columns).
If the network has a non-trivial cut-edge (and therefore no 4-degree blob),
then this matrix has a single row, and describes the split of the cut-edge.

Quartet CFs are estimated by simulating gene trees under the coalescent model
using [PhyloCoalSimulations](https://cecileane.github.io/PhyloCoalSimulations.jl/dev).

output:
- class: number of splits in the displayed trees.
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
- flag: true if minor hybrid edges near the root were discarded to speed up the
  class calculation. Some displayed splits may be missing, and the true class
  may be greater than the output class (first item).

# warning

The function can take a long time if the network has many reticulations.
An approximation can be used: set argument `threshold_h` to some integer
(e.g. 15). Then, if h>15, minor hybrid edges are deleted (starting from those
closest to the root) until h<=15.
If 3 splits are found, then the full network must also have 3 splits
(there cannot be more!). If 1 or 2 splits are found only, then the full network
does have this/these split(s), but may have more. In that case, you may want
to re-run the function with the default `threshold_h=nothing` to be 100% sure
about the splits in the displayed trees.

# examples

```julia
julia> include("find_3blobqCF.jl")

julia> net_3blob = readTopology("((B:0.6,((A:0.4,C:0.4):0.1)#H1:0.1::0.51):1.0,(#H1:0.1::0.49,O:0.6):1.0);");

julia> ngenes = 10_000; # number of genes to be simulated

julia> ns, qCF, hwc, df, flag = quartettype_qCF(net_3blob, ngenes; seed=321, verbose=true);
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
   1 │ A       B       C       O        0.3459   0.2954   0.3587  10000.0

julia> hwc # hardwired clusters, 4 taxa in same order as t1-t4 in data frame above
1×4 BitMatrix:
 1  0  1  0

julia> qCF # order corresponding to 1st split in hwc matrix
(split1 = 0.2954, split2 = 0.3587, split3 = 0.3459)

julia> n_major = Int(round(qCF[:split1] * ngenes, digits=8));

julia> n_alt1  = Int(round(qCF[:split2] * ngenes, digits=8));

julia> n_alt2  = Int(round(qCF[:split3] * ngenes, digits=8));

julia> using HypothesisTests

julia> bt = BinomialTest(n_alt1, n_alt1 + n_alt2); # equal alternative CFs?

julia> pvalue(bt) # the data are consistent with equal alternative CFs.
0.1302795750013259

julia> bt = BinomialTest(n_major, ngenes, 1/3); # is CF_major < 1/3?

julia> pvalue(bt, tail=:left) # reject: so this is a pathological network
2.377765063361621e-16
```

Here is another example where the network has a 4-degree blob, and is outerplanar.
This case is not pathological because the first 2 quartet CFs, which correspond
to the 2 splits in the displayed trees, indeed have the 2 largest CFs,
so they tell us the correct circular ordering.

```julia
julia> net_4blob = readTopology("((((T:0.5)#H1:0.6::0.51,C:1.1):0.7,(#H1:0.0::0.49,E:0.5):1.3):1.0,O:2.8);");

julia> ns, qCF, hwc, df, flag = quartettype_qCF(net_4blob, ngenes; seed=321, verbose=false);

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
   1 │ C       E       O       T        0.1276   0.4882   0.3842  10000.0

julia> qCF
(split1 = 0.3842, split2 = 0.4882, split3 = 0.1276)

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
function quartettype_qCF(net::HybridNetwork, 
        nsim=200; seed=nothing, verbose=true,
        threshold_h=nothing, threshold_gamma=0.10)

    taxonlist = sort(tipLabels(net))
    length(taxonlist) == 4 || error("there aren't 4 tips: $taxonlist")
    # collect splits appearing in any displayed tree
    # their splits only depend on their unrooted topologies, so simplify the network first
    net2 = deepcopy(net) # keep 'net' intact for simulation below
    PhyloNetworks.removedegree2nodes!(net2) # unroots the network, if degree-2 root
    shrink3cycles!(net2, true) # true to unroot the network, e.g. if root is part of a 2- or 3-cycle
    h_true = net2.numHybrids
    h_used = h_true
    if !isnothing(threshold_h) && h_true > threshold_h
      # then delete reticulations with "small" γ: not the best choice though
      # deleteHybridThreshold!(net2, threshold_gamma, false, true)
      @info "h=$h_true, will delete early minor hybrid edges"
      while h_used > threshold_h
        hnum = deleteearlyminorhybrid!(net2)
        PhyloNetworks.removedegree2nodes!(net2)
        shrink3cycles!(net2, true)
        h_used = net2.numHybrids
        # @info "hybrid edge $hnum deleted. h=$h_used after shrinking small cycles."
      end
    end
    flag_class = (h_true != h_used)
    # below: easy, but too slow and too memory-hungry for what we need.
    # dtree = displayedTrees(net2, 0.0)
    # instead: custom function to extract unrooted displayed tree topologies
    dtree = displayed_unrootedtreetopologies!(net2)
    mat = BitMatrix(undef, (0,4)) # initialize: 1 row per split
    for tree in dtree
        split = treesplit(tree, taxonlist)
        oldsplit = any( x-> isequal_split(split, x), eachrow(mat))
        if !oldsplit
          mat = [mat; split']
        end
        size(mat,1) >= 3 && break
    end
    nsplits = size(mat,1)

    # if 3 splits, no qCFs can be anomalous, so need to simulate gene trees:
    #              and even if flag_class is true, class 3 is in fact correct
    nsplits == 3 && return nsplits, (split1=-1, split2=-1, split3=-1), mat, DataFrame(), false

    # otherwise: simulate gene trees!
    isnothing(seed) || Random.seed!(seed)
    treelist = simulatecoalescent(net, nsim, 1)
    obsCF, t = countquartetsintrees(treelist; showprogressbar=verbose)

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
    return nsplits, (split1=df[1,o[1]], split2=df[1,o[2]], split3=df[1,o[3]]), mat, df, flag_class
end

isequal_split(x, y) = x == y || x == .!y

"""
    treesplit(tree::HybridNetwork, taxa)

BitVector representing the split in `tree`. Assumption *not* checked:
`tree` is a tree, it has 4 taxa, it is resolved (not a star polytomy).
"""
function treesplit(tree::HybridNetwork, taxa)
  mat = Bool.(hardwiredClusters(tree, taxa)[1:end,2:(end-1)])
  # find first non-trivial split that partitions the 4 taxa in 2 vs 2
  tokeep = findfirst( x -> sum(x) == 2, eachrow(mat))
  isnothing(tokeep) && error("none of the clades have 2 taxa")
  return(mat[tokeep,:])
end

"""
    deleteearlyminorhybrid(net)

Delete the minor hybrid parent edge of the earliest hybrid node.
Assumptions *not* checked: `net` has edge lengths, and is time-consistent.
Output: number of the hybrid edge that was deleted (which no longer exists).
"""
function deleteearlyminorhybrid!(net::HybridNetwork)
  heights = PhyloNetworks.getHeights(net) # runs preorder! by default, good Here
  hyb_idx = findfirst(n -> n.hybrid, net.nodes_changed)
  isnothing(hyb_idx) && error("no hybrid: can't delete an early hybrid edge")
  hyb_edge = PhyloNetworks.getMinorParentEdge(net.nodes_changed[hyb_idx])
  hyb_num = hyb_edge.number
  PhyloNetworks.deletehybridedge!(net, hyb_edge, false, true)
  return hyb_num
end

"""
    displayed_unrootedtreetopologies!(net)

Vector of all unrooted tree topologies displayed in the network.
Some may be repeated. Their edge lengths and root should be ignored.
Compared to the output of `displayedTrees`, some displayed trees with the same
unrooted topologies, but differ in their root or edge lengths could be
represented only once in the output of `displayed_unrootedtreetopologies!`.
"""
function displayed_unrootedtreetopologies!(net)
  trees = HybridNetwork[]
  displayed_unrootedtreetopologies!(trees, net)
end
function displayed_unrootedtreetopologies!(trees, net)
  if net.numHybrids==0
    # warning: no update of edges' containRoot (true) or edges' and nodes' inCycle (-1)
    push!(trees, net)
  else
    netmin = PhyloNetworks.displayedNetworks!(net, net.hybrid[1], false, true, false, false)
    PhyloNetworks.removedegree2nodes!(netmin) # in case degree-2 root
    PhyloNetworks.removedegree2nodes!(net)
    deleteaboveLSA!(netmin)
    deleteaboveLSA!(net)
    shrink3cycles!(netmin, true)
    shrink3cycles!(net, true)
    displayed_unrootedtreetopologies!(trees, net)
    displayed_unrootedtreetopologies!(trees, netmin)
  end
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
