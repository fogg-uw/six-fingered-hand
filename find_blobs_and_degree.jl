using PhyloNetworks

# given a network "net",
# find the degree (number of outgoing cut edges) of each blob in net 
# by counting A = number of edges in the blob and
# B = number of edges adjacent to some node in the blob.
# degree = B - A.

function blob_degree(net::HybridNetwork)

    blobs = PhyloNetworks.biconnectedComponents(net, true) #true: ignore trivial blobs.  docs suggest shoud be ignoreTrivial=true?

    blob_degrees = Int16[]

    for blob in blobs

        nodes = PhyloNetworks.ANode[]
        for edge in blob
            nodes = [nodes; edge.node]
        end
        nodes = unique(nodes)

        bigblob = PhyloNetworks.Edge[]
        for node in nodes
            bigblob = [bigblob; node.edge]
        end
        bigblob = unique(bigblob)

        blob_neighbors = setdiff(bigblob, blob)
        degree = length(blob_neighbors)
        blob_degrees = [blob_degrees; degree]

    end

    return blobs, blob_degrees

end

