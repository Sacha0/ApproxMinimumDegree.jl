module NWLLPriorityQueues
const DEBUG = false

### Node weight tracking data structure ###
#
# We track the weight (minimum degree, minimum fill, approximate minimum degree, or
# what have you) of each live supervariable (quotient graph node) in the following
# data structure. Specifically, we need be able to
# 1 --- Efficiently identify and remove a minimum-weight node from the data structure
#       in each gross iteration of the algorithm.
# 2 --- For each such gross iteration of the algorithm, we need efficiently update
#       the node weights of (typically several) nodes after eliminating the chosen
#       minimum-weight node from the quotient graph.
# 3 --- For each such gross iteration of the algorithm, we need efficiently remove
#       nodes (often several) that are merged into other nodes (superveriables)
#       due to indistinguishability.
# Additionally, once we've constructed the data structure, we never insert new nodes
# (outside of updating, where updating involves insertion).
#
# This seems to be the traditional data structure in the local preorderings literature

immutable NWLLPriorityQueue{T <: Signed}
    # Node weight linked-list priority queue
    order::T
    minweight::Vector{T}
    weights::Vector{T}
    fwdlinks::Vector{T}
    bwdlinks::Vector{T}
    headlinks::Vector{T}
    function NWLLPriorityQueue(order::T)
        minweight = T[order+1]
        fwdlinks = zeros(T, order) # FUTURE necessary to zero-initialize? no
        bwdlinks = zeros(T, order) # FUTURE necessary to zero-initialize? no
        headlinks = zeros(T, order) # FUTURE necessary to zero-initialize? yes
        weights = Vector{T}(order)
        new(order, minweight, weights, fwdlinks, bwdlinks, headlinks)
    end
end

@inline getnodeweight(pq::NWLLPriorityQueue, node) = pq.weights[node]
@inline setnodeweight!(pq::NWLLPriorityQueue, node, weight) = pq.weights[node] = weight

function enqueue!{T}(pq::NWLLPriorityQueue{T}, node::T, weight::T)
    fwdlink = pq.headlinks[weight+1]
    pq.headlinks[weight+1] = node
    if fwdlink != 0
        pq.bwdlinks[fwdlink] = node
    elseif pq.minweight[1] > weight
        pq.minweight[1] = weight
    end
    pq.fwdlinks[node] = fwdlink
    pq.bwdlinks[node] = 0
    pq.weights[node] = weight
end

function dequeue!{T}(pq::NWLLPriorityQueue{T}, node::T)
    weight = pq.weights[node]
    fwdlink = pq.fwdlinks[node]
    bwdlink = pq.bwdlinks[node]
    fwdlink != 0 && (pq.bwdlinks[fwdlink] = bwdlink)
    bwdlink != 0 ? (pq.fwdlinks[bwdlink] = fwdlink) : (pq.headlinks[weight+1] = fwdlink)
    if bwdlink == 0 && fwdlink == 0 && weight == pq.minweight[1]
        # FUTURE Get rid of unnnecessary branch in preceding comparison
        while true
            pq.minweight[1] += 1
            pq.minweight[1] == pq.order && break
            pq.headlinks[pq.minweight[1]+1] != 0 && break
        end
    end
    return node, weight
end
function dequeue!{T}(pq::NWLLPriorityQueue{T})
    node = pq.headlinks[pq.minweight[1]+1]
    dequeue!(pq, node)
end

@inline function update!{T}(pq::NWLLPriorityQueue{T}, node::T, weight::T)
    dequeue!(pq, node)
    enqueue!(pq, node, weight)
end
# FUTURE Write specialized update routine.

end