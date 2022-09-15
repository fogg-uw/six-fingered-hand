using PhyloNetworks

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

