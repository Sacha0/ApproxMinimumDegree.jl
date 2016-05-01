module ApproxMinimumDegree
const DEBUG = false

include("NodeWeightPQs/NWLLPriorityQueues.jl")
import .NWLLPriorityQueues
include("QuotientGraphs/QuotientGraphs.jl")
import .QuotientGraphs

function amd{Ti,Tv}(C::SparseMatrixCSC{Tv,Ti})
    matorder = C.n

    ### DATA STRUCTURE SETUP ###
    
    # Quotient Graph (elimination and supervariable quotient)
    #
    # As with most local preordering codes, this code maintains a quotient graph model
    # of the elimination graph rather than an explicit representation of the elimination
    # graph. Here by quotient graph we mean quotient in both elimination and
    # supervariable senses.
    #
    # The quotient graph data structure we've abstracted into the module QuotientGraphs.
    # In this method we primarily manipulate it via a relatively opaque, high-level
    # interface. For more on the QuotientGraph data structure, see Quotientgraphs.
    #
    # Build the initial quotient graph from input matrix C. The method below currently
    # assumes that C is symmetric and conatins all entries, not merely an upper or lower
    # triangle.
    qg = QuotientGraphs.QuotientGraph(C)
    # FUTURE Allow input matrix C to be asymmetric or contain only the upper or lower
    # triangle of a symmetric matrix, and build the initial quotient graph for (C + C^T)
    
    # Node weight tracking data structure
    #
    # We track the weight (minimum degree, minimum fill, approximate minimum degree, or
    # what have you) of each live supervariable (quotient graph node) in a specialized
    # priority-queue-like data structure. Here we primarily manipulate this data structure
    # via a relatively opaque, high-level interface. See the NWLLPriorityQueues module
    # for more information on the data structure.
    #
    # Initialize the node weight priority queue and then populate it
    nodeweightpq = NWLLPriorityQueues.NWLLPriorityQueue{Ti}(matorder)
    for k in 1:matorder
        NWLLPriorityQueues.enqueue!(nodeweightpq, k, qg.setstorage[qg.setpointers[k]])
    end
    # FUTURE Initializing the priority queue in one shot from the quotient graph may be
    # more efficient than building it incrementally as we do just above.
     
    # Allocate various workspaces
    #
    # Generic workspaces for various operations (for example by the QuotientGraph)
    workpoolA = Vector{Ti}(matorder)
    workpoolB = Vector{Ti}(matorder)
    #
    # Workspaces used for specific steps in approximate external degree calculation
    # W = Vector{Ti}(n)
    # We now attach W to workpoolA as operations on either are temporally segregated.
    W = workpoolA
    # lL = Vector{Ti}(n)
    # We now attached lL to nodeweightpq.weights as entries go into lL only after the
    # corresponding weight no longer matters.
    lL = nodeweightpq.weights        
    
    ### MAIN ELIMINATION LOOP ###
    elimstepn = 1
    while elimstepn <= matorder
        # Identify a minimum-weight supervariable / quotient graph node by querying the
        # the node weight priority queue. This minimum-weight supervariable / node,
        # commonly the "pivot" p in the literature, we eliminate from the quotient
        # graph in this iteration of the main loop and then add to the ordering.
        p, pw = NWLLPriorityQueues.dequeue!(nodeweightpq)
        DEBUG && println("Beginning elimination step $elimstepn. Pivot $p of weight $pw.")
                
        # Transform the quotient graph, eliminating supervariable p. Specifically, we
        # transform supervariable / quotient graph node p into quotient graph element p.
        QuotientGraphs.elementize!(qg, nodeweightpq, p, workpoolA, workpoolB)
        
        # Update the remaining supervariables' / quotient graph nodes' weights.
        # Specifically, we need udpate the approximate external degree of each supervariable
        # i in [the set of supervariables to which now-element p links, L_p].
        #
        # First we compute the number of simple variables that now-element p links to, i.e.
        # the simple-variable count of set L_p. We calculate this value once and cache it
        # as we potentially use it many times and it does not change during its lifetime
        # (as L_p does not change the set of _simple_ variables it contains, though the
        # set of _super_variables it contains may change, till another element absorbs
        # element p and we no longer need this information). For efficiency we've fused
        # the loop over L_p calculating this quantity into similar loops over L_p for
        # Alg. 2 described below. The old, unfused code read:
        # lLp = 0
        # for kLp in QuotientGraphs.rangeL(qg, p)
        #     iLp = qg.setstorage[kLp]
        #     lLp += QuotientGraphs.numsimplevars(qg, iLp)
        # end
        # lL[p] = lLp
        #
        # Second, for each supervariable i in L_p, and for each element e in [the set of
        # elements to which supervariable i links, E_i], we must calculate
        # length(setdiff(L_e, L_p)). We accomplish this via Algorithm 2 from the seminal
        # AMD paper by Amestoy, Davis, and Duff. Notes:
        #
        # Algorithm 2 calculates the relevant quantity for all live elements,
        # but it does so implicitly. That is, it explicitly calculate the relevant quantity
        # for all (e in E[i] for i in L[p]). Calculating the relevant quantity
        # for live elements outside that set is trivial, as for element e outside that set
        # L[e] and L[p] are necessarily disjoint, not being connected to any of the same
        # live supervariables. Which...
        #
        # .. is fine, as updating the approximate external degree of supervariables i
        # in L[p] only requires the relevant quantity for (e in E[i] for i in
        # L[p]), i.e. only where the calculation is nontrivial. We restrict our calculation
        # to only the required subset.
        #
        # We begin Algorithm 2 with the implicit marking step, setting each W[e] = -1
        for kLp in QuotientGraphs.rangeL(qg, p)
            iLp = qg.setstorage[kLp]
            for kEi in QuotientGraphs.rangeE(qg, iLp)
                eEi = qg.setstorage[kEi]
                W[eEi] = -1
            end
        end
        # FUTURE We should be able to eliminate this marking step by introducing slight
        # additional complexity into the latter part of Alg. 2 and using some assumptions
        # on the contents of the workspace in which we place W.
        #
        # Then we perform Algorithm 2 proper, calculating the set of values simultaneously.
        #
        # NOTE We only access lL[e = p] in this loop to set W[e = p], which hypothetically
        # becomes simplevarcount(setdiff(L_{e = p}, L_p)), but of course we know that value
        # trivially (zero) and, though we might use it below, we don't actually need it.
        # This allowed us to fuse the lL[p] calculation discussed above with the following
        # loop. This way we avoid evaluating QuotientGraphs.numsimplevars(qg, i) twice
        # for each i in L_p, we may geet better cache behavior once we fuse lL with
        # qg.nodestatus instead of nodeweightpq.weights, and we are free to get rid of
        # the above marking loop if we find a nice way to do so.
        lLp = zero(Ti)
        for kLp in QuotientGraphs.rangeL(qg, p)
            i = qg.setstorage[kLp]
            lenIi = QuotientGraphs.numsimplevars(qg, i)
            lLp += lenIi
            for kEi in QuotientGraphs.rangeE(qg, i)
                e = qg.setstorage[kEi]
                We = W[e]
                if We < 0
                    W[e] = lL[e] - lenIi
                else
                    W[e] = We - lenIi
                end
            end
        end
        W[p] = 0
        lL[p] = lLp
        #
        # Third, with the quantities calculated above on hand, we calculate and update
        # the approximate external degrees of the supervariables i in L_p.
        for kLp in QuotientGraphs.rangeL(qg, p)
            iLp = qg.setstorage[kLp]
            # Three terms contribute to the approximate external degree for supervariable
            # i. We calculate each in turn, then the approximate external degree proper.
            #
            # Term 1, the simple-variable count of setdiff(L_p, I_i). We know that
            # supervariable i is an entry of set L_p, and hence in terms of simple
            # variables all entries in I_i are entries of L_p. Hence, having calculated
            # the simple-variable count of L_p above, we can calculate the simple-var count
            # of setdiff(L_p, I_i) trivially as simplevarcount(L_p) - simplevarcount(I_i):
            lengthLpLessIi = lLp - QuotientGraphs.numsimplevars(qg, iLp)
            #
            # Term 2, the simple-variable count of setdiff(A_i, I_i). If we are
            # maintaining A_i properly, we know that supervariable i should not be an
            # entry of set A_i; that is, when we build our QuotientGraph, we exclude
            # i from A_i, and if we maintain A_i properly, we should not reintroduce i
            # into A_i, and when merging supervariables etcetera we carefully strip
            # entries in I_i from A_i. So setdiff(A_i, I_i) should just be A_i in our
            # case. Hence the simple-variable count of that setdiff is just the
            # simple-variable count of A_i, which we calculate readily.
            #
            # Term 3, sum_{e in E_i} simplevarcount(setdiff(L_e, L_p)). We precalculated
            # each term in the sum above via Algorithm 2. Hence here we need simply sum
            # the appropriate precalculated terms.
            #
            # Given that A_i and E_i are stored contiguously in qg.setstorage and in that
            # order, we calculate terms 2 and 3 fused as follows for efficiency.
            #
            rangeAi, rangeEi = QuotientGraphs.rangesAandE(qg, iLp)
            # FUTURE This range-pair-retrieval operation is slow. Fix.
            lengthAiLessIi = zero(Ti)
            for kAi in rangeAi
                jAi = qg.setstorage[kAi]
                lengthAiLessIi += QuotientGraphs.numsimplevars(qg, jAi)
            end
            sumEi_lengthLeLessLp = zero(Ti)
            for kEi in rangeEi
                eEi = qg.setstorage[kEi]
                # eEi == p && continue # W[p] == 0 from above
                sumEi_lengthLeLessLp += W[eEi]
            end
            # We combine the three terms calculated above into the approximate
            # external degree proper.
            oldapproxextdeg = NWLLPriorityQueues.getnodeweight(nodeweightpq, iLp)
            bound1_maxextdegree = matorder - elimstepn
            bound2_crudemaxfill = oldapproxextdeg + lengthLpLessIi
            bound3_tightextdegest = lengthLpLessIi + lengthAiLessIi + sumEi_lengthLeLessLp
            newapproxextdeg = min(bound1_maxextdegree, bound2_crudemaxfill, bound3_tightextdegest)
            #
            # Finally, we update this supervariable's / quotien graph node's
            # approximate external degree / weight in the node weight priority queue.
            NWLLPriorityQueues.update!(nodeweightpq, iLp, newapproxextdeg)
        end        
        elimstepn += QuotientGraphs.numsimplevars(qg, p)
    end
    
    postorder!(qg, workpoolA, workpoolB)    
    return workpoolB
end

"""
    postorder!{T<:Signed}(qg::QuotientGraph{T}, perm::Vector{T}, iperm::Vector{T})
    
Given a `QuotientGraph` on which we completed elimination/ordering, returns an assembly-tree
postordering as a permutation vector in `perm` and as an inverse permutation vector in
`iperm`.
"""
function postorder!{T<:Signed}(
        qg::QuotientGraphs.QuotientGraph,
        perm::Vector{T},
        iperm::Vector{T})
    supervarcards = qg.nodestatus
    rootwardlinks = qg.setpointers
    leafwardlinks = perm
    # On entry, rootwardlinks and supervarcards together define the assembly tree generated
    # during the ordering process. Specifically, for index i into supervarcards,
    # superacards[i] indicates one of two possible things:
    #
    # (1) supervarcards[i] == 0 indicates that supervariable i was amalgmated into another
    # supervariable prior to (mass) elimination during the ordering process.
    #
    # (2) supervarcards[i] (== k) > 0 indicates that supervariable i was representative
    # when it was eliminated, and hence became an element. k is the number of simple
    # variables that supervariable i represented on (mass) elimination.
    #
    # For index i into rootwardlinks, rootwardlinks[i] indicates one of three things:
    #
    # (1) rootwardlinks[i] >= 0 indicates that as an element i was never absorbed into
    # another element, and hence constitutes a root of the assembly tree.
    #
    # (2) If supervarcards[i] == 0 and hence as above i was a supervariable amalgamated
    # into another supervariable prior to (mass) elimination during the ordering process,
    # then rootwardslinks[i] (== -j) < 0 indicates that supervariable i was amalgamated
    # into supervariable j.
    #
    # (3) If supervarcards[i] > 0 and hence as above supervariable i was representative
    # when it was eliminated and so became an element, then rootwardlinks[i] (== -j) < 0
    # indicates that element i was absorbed into elemnt j during the ordering process.
    #
    # In this fashion rootwardlinks and supervarcards together define the assembly tree
    # generated during the ordering process by way of child-to-parent links.
    #
    # Our objective is to generate a postordering of this assembly tree as a permutation
    # vector `perm` and an inverse permutation vector `iperm`. We do this via a
    # depth-first-search of the assembly tree (which can be a forest, and we handle that
    # case as well). To perform an efficient DFS, we need parent-to-child links in addition
    # to the child-to-parent links we have. Hence this method involves two stages: In
    # the first stage, we build such parent-to-child links in leafwardlinks. In the second
    # stage, we perform the depth-first-search of the assembly tree.
    #
    # Moreover, we need not only parent-to-child links, but we need an ordering among the
    # children of a given parent. Specifically, we need that those children corresponding to
    # assembly operations be 'elder' to those corresponding to supervariable amalgamations.
    # Hence, for nodes with younger siblings, we replace the link to the node's parent in
    # rootwardlinks to a link to the node's next-younger sibling. Having this information
    # allows us to efficiently walk the tree ordering elder siblings (assembly operations)
    # before younger siblings (supervariable amalgamations). To effect this, we build
    # leafwardlinks in two stages, first sweeping over amalgamated supervariables, and then
    # over elements (assembly operations), ensuring the ordering we need.
    #
    # For more information on this process, see the paper and technical report on MA27
    # from Duff and Reid.
    #
    # Additionally, as a side effect of stage one we build a linked list of assembly tree
    # roots. nextroot points to the first root. rootwardlinks[nextroot] then points to the
    # second root if it exists, rootwardslinks[rootwardlinks[nextroot]] then points
    # to the third root if it exists, and so on. A zero in the chain indicates its end.

    # Stage one / initialization.
    nextroot = 0
    fill!(leafwardlinks, 0)
    # Stage one / substage one, building parent-to-child and child-to-ysibling links
    # for amalgamated supervariables.
    for child in 1:length(supervarcards)
        supervarcards[child] == 0 || continue # Skip all but amalgamated supervariables
        nextroot = parentchildlink!(child, rootwardlinks, leafwardlinks, nextroot)
    end
    # Stage one / substage two, building parent-to-child and child-to-ysibling links
    # for absorbed elements.
    for child in 1:length(supervarcards)
        supervarcards[child] != 0 || continue # Skip all but elements
        nextroot = parentchildlink!(child, rootwardlinks, leafwardlinks, nextroot)
    end
    
    # Stage two, depth-first-search through the assembly tree to generate the postorder.
    korder = 1
    while nextroot != 0
        # We begin the DFS from the assembly-tree root identified in nextroot. Presently
        # rootwardlinks[nextroot] either contains 0 indicating nextroot is the last root
        # in the rootchain, or contains a positive integer identifying the next root in
        # the root chain. During the DFS, we need rootwardlinks[root] = 0 such that we
        # can terminate DFS once we've returned to the root. Hence we need set
        # rootwardlinks[nextroot] = 0, but we need the present value after we've completed
        # this DFS. Hence, retrieve rootwardlinks[nextroot], store it in nextroot after
        # seeding the present DFS, and then set rootwardlinks[nextroot] = 0.
        presnode = nextroot
        nextroot = rootwardlinks[presnode]
        rootwardlinks[presnode] = 0
        while true
            child = leafwardlinks[presnode]
            if child != 0
                # If the present node has a child, step to the child and then burn the
                # bridge so that we cannot walk down this path again.
                leafwardlinks[presnode] = 0
                presnode = child
                continue
            else
                # If the present node has no child, add this node to both representations
                # of the postordering.
                leafwardlinks[presnode] = korder # permutation representation
                iperm[korder] = presnode # inverse permutation represeentation
                korder += 1
                # Retrieve the present node's parent or younger sibling.
                ysiborpar = rootwardlinks[presnode]
                # If this node has no parent or younger sibling, we've returned to the root,
                # in which case terminate the search.
                ysiborpar >= 0 && break
                # Otherwise, step to this node's parent or sibling and continue the search.
                presnode = -ysiborpar
            end
        end        
    end
end
    
@inline function parentchildlink!(child, rootwardlinks, leafwardlinks, nextroot)
    parent = rootwardlinks[child]
    if parent >= 0
        # rootwardlinks[child] >= 0 indicates child was not an amalgamated supervariable or
        # absorbed element, but rather an ultimately unabsorbed element and hence assembly
        # tree root node. Add this node to the chain of roots.
        rootwardlinks[child] = nextroot
        return child # nextroot = child
    else # parent < 0
        # If leafwardlinks[parent] already contains a link (i.e. != 0 and so parent must
        # already have at least one child attached), then child must have a younger sibling.
        # In this case, we replace child's child-to-parent link in rootwardlinks with a
        # child-to-ysibling link, and then place a parent-to-child link to this child in
        # leafwardlinks over the younger sibling's parent-to-child link.
        ysibling = leafwardlinks[-parent]
        leafwardlinks[-parent] = child
        ysibling != 0 && (rootwardlinks[child] = -ysibling)
    end
    return nextroot
end

end
