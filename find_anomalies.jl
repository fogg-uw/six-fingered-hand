#= code used by `2_extract_quartet_subnetworks.jl` to determine:

1. the number of splits displayed in 4-taxon networks.
   class-1 networks have 1 split.
   class-2 networks have 2 splits (and must have a 4-blob).

2. if a class-1 4-taxon network is anomalous, in the sense that
   it has a cut edge ab|cd but the associated quartet CF
   is the smallest: CF(ab|cd) < CF(ac|bd) = CF(ad|bc)

2. if a class-2 4-taxon set is anomalous,
   in the sense that it displays 2 splits, say ab|cd and ac|bc,
   but the third split ad|bc has a quartet CF that is not the smallest:
   CF(ad|bc) > either CF(ab|cd) or CF(ac|bd).

Tools: find blobs, find the degree of each blob, find all splits displayed
by a 4-taxon network (without blindly extracting all displayed trees),
simulate gene trees to estimate quartet concordance factors.
=#

using Random
using PhyloNetworks
using PhyloCoalSimulations
using DataFrames
using HypothesisTests

"""
    blob_degree(net)

Find the blobs and their degrees (number of outgoing cut edges), by counting
A = number of edges in the blob and
B = number of edges adjacent to some node in the blob.
Then degree = B - A.

Output: (blobs, blob_degrees)
"""
function blob_degree(net::HybridNetwork)
    blobs = biconnectedComponents(net, true) # true: to ignore trivial blobs
    blob_degrees = Int16[]
    for blob in blobs
        # list of nodes adjacent to the blob
        nodes = PhyloNetworks.ANode[]
        for edge in blob
            append!(nodes, edge.node)
        end
        nodes = unique(nodes)
        # list edges incident to nodes adjacent to blobs: either inside or touching the blob
        bigblob = PhyloNetworks.Edge[]
        for node in nodes
            append!(bigblob, node.edge)
        end
        bigblob = unique(bigblob)
        # find cut edges adjacent to the blob
        blob_neighbors = setdiff(bigblob, blob)
        degree = length(blob_neighbors)
        push!(blob_degrees, degree)
    end
    return blobs, blob_degrees
end

isequal_split(x, y) = x == y || x == .!y
splittype(s) = ( s== [1,1,0,0] || s == [0,0,1,1] ? 0 :
                (s== [1,0,1,0] || s == [0,1,0,1] ? 1 :
                                                   2  ))
splitsymbol(st) = (st==0 ? :CF12_34 : (st==1 ? :CF13_24 : :CF14_23))

"""
    has32cycle!(net, bcc=missing, blobdegree=missing, blobexit=missing)
    has32cycle(blobdegree, exitnodes, blobedges)

Whether network `net` or one of its blob with known `exitnodes` contains
a 3_2 cycle as a subnetwork.
Assumptions *not* checked: `net` should have 4 taxa, no 4-blob, and be unrooted.

Strategy, to avoid extracting all subnetworks:
1. a lowest hybrid node has 2 descendants => has 32 cycle
2. exit nodes have at most 2 descendants total => does NOT have 32 cycle
3. 1 single hybrid in the blob with 1 descendant => NOT a 32 cycle
4. else: remove a lowest hybrid node and determine recursively.
   when doing so: the network is deepcopied to remove the minor hybrid parent,
   and is modified to remove the major hybrid parent edge.

The first method calls the second, returns true or false.
The second method returns true in case 1, false in cases 2 and 3, missing in case 4.
"""
function has32cycle!(net::HybridNetwork)
  bcc, blobdegree = blob_degree(net)
  blobexit = PhyloNetworks.biconnectedcomponent_exitnodes(net, bcc, false)
  res_byblob = [has32cycle(blobdegree[i],blobexit[i],bcc[i]) for i in eachindex(bcc)]
  res = any(res_byblob)   # at least 1 of the blobs contains a 3_2 cycle
  !ismissing(res) && res && return true
  res = all(.!res_byblob) # none of the blobs contain a 3_2 cycle
  !ismissing(res) && res && return false
  # by now: all blobs have either false or missing, and at least one blob has missing
  # missing: when at most 2 descendants, and only 1 per hybrid exit.
  # unrooting net ensures that there's at most 1 blob with missing info
  missingbi = findall(ismissing, res_byblob)
  length(missingbi) == 1 || error("more than 1 blob has missing 32cycle info: res_byblob=$res_byblob, blobexit=$blobexit, net=$(writeTopology(net))")
  bi = missingbi[1]
  ni = findfirst(en -> en.hybrid, blobexit[bi])
  isnothing(ni) && error("blob $bi without an exit hybrid node: blobexit=$blobexit, net=$(writeTopology(net))")
  hybnum = blobexit[bi][ni].number # number of exit hybrid node
  hybind = findfirst(n -> n.number == hybnum, net.node) # index: nodes will be deepcopied
  isnothing(hybind) && error("hybnode not found in net: hybnum=$hybnum, index $ni in blobexit[bi]=$(blobexit[bi]), net=$(writeTopology(net))")
  # delete minor parent: return true if subnetwork has 3_2 cycle
  minornet = deepcopy(net)
  hp = PhyloNetworks.getMinorParentEdge(minornet.node[hybind])
  PhyloNetworks.deletehybridedge!(minornet, hp, false, true)
  has32cycle!(minornet) && return true
  # delete major parent: return true if subnetwork has 3_2 cycle
  hp = PhyloNetworks.getMajorParentEdge(net.node[hybind])
  PhyloNetworks.deletehybridedge!(net, hp, false, true)
  return has32cycle!(net)
end

function has32cycle(bdegree::Integer, exitnodes::AbstractVector, blobedges::AbstractVector)
  bdegree < 3 && return false
  totaldes = 0
  for en in exitnodes
    # find exit edge: *the* child edge of exit node 'en' that is not in the blob
    exitedge = nothing # child edge not in the blob
    for ee in en.edge
      PhyloNetworks.getParent(ee) === en || continue
      ee in blobedges && continue
      exitedge = ee
      break
    end
    ndes = length(PhyloNetworks.descendants(exitedge))
    en.hybrid && ndes == 2 && return true
    totaldes += ndes
  end
  # by now, all exit hybrid have 1 descendant only
  totaldes < 3 && return false
  # if the blob is a cycle: not a 3_2 cycle
  sum(e.hybrid for e in blobedges) <= 2 && return false
  return missing
end

"""
    quartettype_qCF(net, nsim=200, inheritancecorrelation=0.0;
                      seed=nothing, verbose=true,
                      threshold_h = Inf,
                      blob_degrees = nothing,
                      attempt_resolve_ambiguity = false)

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
- symbol `:anomalous` `:ambiguous` or `:good` from comparing the confidence interval
  for the quartet CF of interest (major or minor) with the anomalous threshold
  (1/3 or either 1st or 2st split CF). If :ambiguous with `nsim` gene trees,
  and if `attempt_resolve_ambiguity` is true,
  the simulation is repeated with `100 * nsim` gene trees to attempt to remove
  the ambiguity (this number of genes is reflected in the data frame).
  `missing` if 3 splits.
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
julia> include("find_anomalies.jl")

julia> net_3blob = readTopology("((B:0.6,((A:0.4,C:0.4):0.1)#H1:0.1::0.51):1.0,(#H1:0.1::0.49,O:0.6):1.0);");

julia> ngenes = 10_000; # number of genes to be simulated

julia> ns, qCF, hwc, df, is32blob, isanomalous, flag, o = quartettype_qCF(net_3blob, ngenes; seed=321, verbose=true);
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

julia> isanomalous
:anomalous

julia> o # order of split1, split2, split3
3-element Vector{Symbol}:
 :CF13_24
 :CF14_23
 :CF12_34

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

julia> confint(bt, level=0.95) # fully below 1/3: anomalous
(0.28646883314426824, 0.30445046101054496)

julia> test_split1_onethird(qCF, ngenes)
:anomalous
```

Here is another example where the network has a 4-degree blob, and is outerplanar.
This case is not pathological because the first 2 quartet CFs, which correspond
to the 2 splits in the displayed trees, indeed have the 2 largest CFs,
so they tell us the correct circular ordering.

```julia
julia> net_4blob = readTopology("((((T:0.5)#H1:0.6::0.51,C:1.1):0.7,(#H1:0.0::0.49,E:0.5):1.3):1.0,O:2.8);");

julia> ns, qCF, hwc, df, is32blob, isanomalous, flag, o = quartettype_qCF(net_4blob, ngenes; seed=321, verbose=false);

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

julia> o
3-element Vector{Symbol}:
 :CF14_23
 :CF13_24
 :CF12_34

julia> qCF
(split1 = 0.3842, split2 = 0.4882, split3 = 0.1276)

julia> isanomalous
:good

julia> n_cf1    = Int(round(qCF[:split1] * ngenes, digits=8));

julia> n_cf2    = Int(round(qCF[:split2] * ngenes, digits=8));

julia> n_minor  = Int(round(qCF[:split3] * ngenes, digits=8));

julia> bt = BinomialTest(n_cf1, n_cf1 + n_minor); # is CF1 > CF_minor?

julia> pvalue(bt, tail=:left) # the data are consistent with CF1 >= CF_minor
1.0

julia> confint(bt, level=0.90) # fully above 0.5: all good
(0.7405349763394173, 0.7606208489695341)

julia> test_split3_split12(qCF, ngenes)
:good

julia> pvalue(BinomialTest(n_cf2, n_cf2 + n_minor), tail=:left) # also with CF2 >= CF_minor
1.0
```
"""
function quartettype_qCF(net::HybridNetwork, 
        nsim=200, inheritancecorrelation=0.0; seed=nothing, verbose=true,
        threshold_h=Inf, blob_degrees=nothing,
        attempt_resolve_ambiguity=false)

    taxonlist = sort(tipLabels(net))
    length(taxonlist) == 4 || error("there aren't 4 tips: $taxonlist")
    if isnothing(blob_degrees)
      blob_degrees = blob_degree(net)
    end
    bcc     = blob_degrees[1]
    bdegree = blob_degrees[2]
    is32blob = false
    flag_class = false
    if all(bdegree .< 4)
    # if the network doesn't have any 4-blob, then it's of class 1:
    # extract hardwired clusters, and determine if there's a 3_2 blob
      netcp = deepcopy(net)
      PhyloNetworks.fuseedgesat!(netcp.root, netcp)
      is32blob = has32cycle!(netcp) #, bcc, bdegree)
      # get the one split from hardwired clusters
      nsplits = 1
      mat = Bool.(hardwiredClusters(net, taxonlist)[1:end,2:(end-1)])
      # find splits that partition the 4 taxa in 2 vs 2
      tokeep = findall( x -> sum(x) == 2, eachrow(mat))
      # if there are several, check they all code for the same split
      length(tokeep)>0 || error("didn't find any 2+2 split in net without 4-blob")
      for i in 2:length(tokeep)
        isequal_split(mat[tokeep[1],:], mat[tokeep[i],:]) ||
          error("several distinct 2+2 splits in net without 4-blob")
      end
      mat = mat[[tokeep[1]],:]
    else
    # if the network has a 4-blob: collect splits appearing in any displayed tree
    # their splits only depend on their unrooted topologies, so simplify the network first
      net2 = deepcopy(net) # keep 'net' intact for simulation below
      unroot_shrink!(net2)
      h_true = net2.numHybrids
      h_used = h_true
      if h_true > threshold_h
        @info "h=$h_true, will delete early minor hybrid edges"
        while h_used > threshold_h
          hnum = deleteearlyminorhybrid!(net2)
          unroot_shrink!(net2)
          h_used = net2.numHybrids
          # @info "hybrid edge $hnum deleted. h=$h_used after shrinking small cycles."
        end
      end
      flag_class = (h_true != h_used)
      # custom function to extract displayed splits
      dsplit = displayed_splits!(net2, taxonlist)
      mat = BitMatrix(undef, (0,4)) # initialize: 1 row per split
      for split in dsplit
      	mat = [mat; split']
      end
      nsplits = size(mat,1)
      if nsplits == 1
        flag_class || @error "found only 1 split on network with a 4-blob, yet no underestimation. ???"
        @error "found only 1 split on 4-blob: flag is true and there's an underestimation for sure"
      end
    end

    # if 3 splits, no qCFs can be anomalous, no need to simulate gene trees:
    #              and even if flag_class is true, class 3 is in fact correct
    nsplits == 3 &&
      return (nsplits, (split1=-1, split2=-1, split3=-1), mat,
              DataFrame(:ngenes => 0),
              false, missing, false, [:CF12_34,:CF13_24,:CF14_23])

    # otherwise: simulate gene trees!
    st1 = splittype(mat[1,:])
    st2 = (nsplits > 1 ? splittype(mat[2,:]) : mod(st1+1,3))
    st3 = setdiff(0:2, [st1,st2])[1]
    o = splitsymbol.([st1,st2,st3])

    df = estimate_qCFs(net, taxonlist, nsim, inheritancecorrelation, seed, verbose)
    qCF = (split1=df[1,o[1]], split2=df[1,o[2]], split3=df[1,o[3]])
    isanomalous = test_anomaly(qCF, nsim, nsplits)
    if attempt_resolve_ambiguity && isanomalous == :ambiguous
      # then estimate with 100 times more genes
      nsim = 100 * nsim
      df = estimate_qCFs(net, taxonlist, nsim, inheritancecorrelation, seed, verbose)
      qCF = (split1=df[1,o[1]], split2=df[1,o[2]], split3=df[1,o[3]])
      isanomalous = test_anomaly(qCF, nsim, nsplits)
    end
    return nsplits, qCF, mat, df, is32blob, isanomalous, flag_class, o
end

"""
    estimate_qCFs(net, taxonlist, nsim, inheritancecorrelation, seed, verbose)

Estimate quartet CFs using simulations from PhyloCoalSimulations.
`net`: should be a 4-taxon network.
`nsim`: number of simulated gene trees.
`taxonlist`: to control that the order of taxa is correctly interpreted.

Output: data frame with 1 row and 8 columns: 4 taxon names, 3 qCFs, and number of simulated genes.
"""
function estimate_qCFs(net, taxonlist, nsim, rho, seed, verbose)
  isnothing(seed) || Random.seed!(seed)
  treelist = simulatecoalescent(net, nsim, 1; inheritancecorrelation=rho)
  obsCF, t = countquartetsintrees(treelist; showprogressbar=verbose)

  taxonlist == t || @error("different order of taxa used by countquartetsintrees (but maybe john's code will help)")

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
  return df
end

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
  mat = mat[tokeep,:]
#  if size(mat,2) != 4
#    mat = transpose(mat) # make sure it is row vector, not column vector
#  end 
  return(mat)
end

"""
    unroot_shrink!(net)

1. remove what's above the LSA
2. unroot the network (for a degree-3 root), then
3. shrink small cycles (of degree 2 or 3).
"""
function unroot_shrink!(net::HybridNetwork)
  deleteaboveLSA!(net)
  rootedges = net.node[net.root].edge
  if length(rootedges) == 2
    if all(e.hybrid for e in rootedges) # should never happen
      # fuseedgesat! (and removedegree2nodes!) would truly error in this case
      @error "after deleting things above the LSA, the root is still connected to 2 hybrid edges..."
    else
      PhyloNetworks.fuseedgesat!(net.root, net)
    end
  end
  shrink3cycles!(net, true) # true: to unroot the network
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
    unroot_shrink!(netmin)
    unroot_shrink!(net)
    displayed_unrootedtreetopologies!(trees, net)
    displayed_unrootedtreetopologies!(trees, netmin)
  end
end

"""
    displayed_splits!(splits, net, taxonlist)

Like displayed_unrootedtreetopologies!, but uses "tricks" that may speed the
calculation when the network is a quartet network.  Specifically, if there is
no 4-blob, then we know that the quartet network has only one split, and it
can be identified quickly.

Returns a vector of binary splits
(e.g. [[1 1 0 0], [1 0 1 0]] for a quartet network displaying the splits
AB|CD and AC|BD but not AD|BC.)
"""
function displayed_splits!(net, taxonlist)
  splits = BitVector[] # could be improved? could preallocate all 12 bits needed
  #splits = BitMatrix(undef, (3,4))
  splits = displayed_splits!(splits, net, taxonlist)
  return(splits)
end
function displayed_splits!(splits, net, taxonlist)
  # if fewer than 3 splits, then we have work to do.  o/w skip to the end
  if length(splits) < 3
    # splits are invariant to root and 3-cycles, so "flatten" them out
    unroot_shrink!(net)
    # if it's a tree or there are no 4-blobs...
    if net.numHybrids==0 || all(blob_degree(net)[2] .< 4)
      # ...then it's class 1.  grab the split!
      newsplit = treesplit(net,taxonlist)
      isempty(splits) && push!(splits, newsplit)
      any(isequal_split(newsplit, x) for x in splits) || push!(splits, newsplit)
    # otherwise (if there's a 4-blob)...
    else
      # ...then break the network apart at the lowest hybrid...
      preorder!(net)
      i = findlast(nn -> nn.hybrid, net.nodes_changed)
      isnothing(i) && error("can't find any hybrid node")
      netmin = PhyloNetworks.displayedNetworks!(net, net.nodes_changed[i], false, true, false, false)
      # ... and then recursion.
      displayed_splits!(splits, net,    taxonlist)
      displayed_splits!(splits, netmin, taxonlist)
    end
  end
  return(splits)
end

"""
    test_anomaly(qCF, ngenes, nsplits)
    test_split1_onethird(qCF, ngenes)
    test_split3_split12( qCF, ngenes)

Binomial confidence intervals, for the alternative hypotheses:
- for test_split1_onethird: qCF1 < 1/3
- for test_split3_split12: qCF3 > qCF1 or qCF3 > qCF2 (using a Bonferroni correction)

Output:
`:anomalous` if there is sufficient evidence for the alternative hypothesis,
`:good` if there is sufficient evidence for the null hypothesis,
`:ambiguous` if the confidence interval includes the threshold (1/3 or 1/5 resp.)
"""
function test_anomaly(qCF, ngenes, nsplits)
  if nsplits == 1
    return test_split1_onethird(qCF, ngenes)
  elseif nsplits == 2
    return test_split3_split12(qCF, ngenes)
  else
    @error "nsplits should be 1 or 2"
    return missing
  end
end
function test_split1_onethird(qCF, ngenes)
  n_cf1 = Int(round(qCF[:split1] * ngenes, digits=8))
  # p = pvalue(BinomialTest(n_cf1, ngenes, 1/3), tail=:left) # CF_major >= 1/3 plausible?
  cf1_lo, cf1_hi = confint(BinomialTest(n_cf1, ngenes, 1/3), level=0.95)
  if 1/3 < cf1_lo
    return :good
  elseif 1/3 > cf1_hi
    return :anomalous
  end
  return :ambiguous
end
function test_split3_split12(qCF, ngenes)
  qCF12 = min(qCF[:split1], qCF[:split2]) # split that will give the smallest p-value
  n_cf12 = Int(round(qCF12        * ngenes, digits=8))
  n_cf3  = Int(round(qCF[:split3] * ngenes, digits=8))
  # p = pvalue(BinomialTest(n_cf12, n_cf12 + n_cf3), tail=:left) # CF1 (and CF2) >= CF_minor plausible?
  # 2*p # for Bonferroni correction
  # level: 1-0.05/2 = 0.975 for the Bonferroni correction
  cf1_lo, cf1_hi = confint(BinomialTest(n_cf12, n_cf12 + n_cf3), level=0.975)
  if 0.5 < cf1_lo
    return :good
  elseif 0.5 > cf1_hi
    return :anomalous
  end
  return :ambiguous
end

nothing # to avoid screen output when we include this file
