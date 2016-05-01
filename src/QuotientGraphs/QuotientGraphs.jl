module QuotientGraphs
const DEBUG = false

import ..NWLLPriorityQueues

"""
    QuotientGraph{T<:Signed}
    
An efficient quotient graph data structure for local preordering algorithms.

`matorder` is the number of nodes in the quotient graph on initialization.

`nodestatus` contains supervariable/element status information: Where `i` is a supervariable
that has neither been amalgamated into another supervariable nor eliminated, `nodestatus[i]`
provides supervariable `i`'s simple-variable count. When such a supervariable `i` is
amalgamated into another supervariable, `nodestatus[i]` changes to `0`, indicating
supervariable `i`'s end in amalgamation. When such a supervariable `i` is instead eliminated
and hence becomes an element, `nodestatus[i]` by default does not change, `nodestatus[i] > 0`
indicating supervariable `i`'s end in elimination. As in such case `nodestatus[i]` merely
need be positive, `nodestatus[i]` can be used to store other positive values if desired; for
example, using `nodestatus[i]` to store un-absorbed element `i`'s simple-variable count `(>0)`
may be useful in some local preordering algorithms.

`setstorage` and `setpointers` together define and contain the linked-supervariable and
linked-element sets for unamalgamated, uneliminated supervariables (respectively `A_i` and
`E_i` in the literature) and the linked-supervariable sets for unabsorbed elements (`L_p`
in the literature). Specifically, for unamalgamated, uneliminated supervariable `i`,
`setpointers[i]` provides the `head` in `setstorage` at which sets `A_i` and `E_i` are
described and stored. `setstorage[headi]` provides the length of set `A_i` and
`setstorage[headi+1]` the length of set `E_i`, followed contiguously by sorted-vector
representations of first `A_i` and then `E_i`. For unabsorbed element `i`, `setpointers[i]`
provides the `head` in `setstorage` at which set `L_p` is described and stored. If
`headi >= stackpointer[1]`, `setstorage[headi]` provides the length of set `L_p` followed
contiguously by a sorted-vector representation of `L_p`; if instead `headi < stackpointer[1]`,
`setstorage[headi]` provides the length of set `L_p` followed contiguously by a zero and
then a sorted-vector representation of `L_p`.

`setpointers` additionally contains information about amalgamated supervariables and
absorbed elements. Where supervariable `i` was amalgamated into supervariable `j`,
`setpointers[i] = -j` indicates this. Where element `i` was absorbed into element `j`,
`setpointers[i] = -j` indicates this.

`setstorage` contains three regions: one preceding `stackpointers[1]`, another from
`stackpointers[1]` to `stackpointers[2]` (the 'L-stack' or 'stack'), and a third from
`stackpointers[2]` on; defining these regions is `stackpointers`'s sole purpose. The first
region stores sets `A_i` and `E_i` for unamalgamated, uneliminated supervariables and also
sets `L_p` for elements for which we were able to store `L_p` 'in-place', i.e. in the storage
for the supervariable that became that element. The second, 'L-stack' region contains sets
`L_p` for elements for which we were not able to store `L_p` 'in-place' in the sense above.
The third region represents open storage into which the L-stack expands as necessary.

The present implementation of `QuotientGraph{T}` allows signed integer `T` only. Though
relaxing this restriction is possible, using the sign bit for marking significantly
simplifies the logic, reduces storage requirements, and improves performance.

This data structure was inspired primarily by the publications: [1] I.S. Duff and J.K. Reid,
"The multifrontal solution of indefinite sparse symmetric linear equations," ACM Transactions
on Mathematical Software (TOMS) 9.3 (1983): 302-325; [2] I.S. Duff and J.K. Reid, "MA27 -
A set of Fortran subroutines for solving sparse symmetric sets of linear equations," UKAEA
Atomic Energy Research Establishment, 1982; and [3] P.R. Amestoy, T.A. Davis, and I.S. Duff,
"An approximate minimum degree ordering algorithm," SIAM Journal on Matrix Analysis and
Applications 17.4 (1996): 886-905.
"""
immutable QuotientGraph{T<:Signed}
    matorder::T
    nodestatus::Vector{T}
    setstorage::Vector{T}
    setpointers::Vector{Int64}
    stackpointers::Vector{Int64}
end

"`numsimplevars(qg, i)` returns the number of simple variables supervariable `i` represents."
@inline numsimplevars(qg, i) = qg.nodestatus[i]

# Functions simplifying iteration over sets A_i, E_i, and L_i associated with nodes in the
# quotient graph follow below. A more opaque, graceful mechanism for iterating over these
# sets would be great. Unfortunately, last I tried building iterators for this task (0.4),
# their performance was significantly worse than the simple if ungraceful facility that
# presently exists revolving around these functions.
"`rangeE(qg, i)` returns the range of set `E_i` in `qg.setstorage`."
@inline function rangeE(qg, i)
    head = qg.setpointers[i]
    lengthA = qg.setstorage[head]
    lengthE = qg.setstorage[head+1]
    firstkE = head + 2 + lengthA
    lastkE = head + 2 + lengthA + lengthE - 1
    return firstkE:lastkE
end
"`rangeA(qg, i)` returns the range of set `A_i` in `qg.setstorage`)"
@inline function rangeA(qg, i)
    head = qg.setpointers[i]
    lengthA = qg.setstorage[head]
    firstkA = head + 2
    lastkA = head + 2 + lengthA - 1
    return firstkA:lastkA
end
"`rangeAE(qg, i)` returns the range containing both sets `A_i` and `E_i` in `qg.setstorage`"
@inline function rangeAE(qg, i)
    head = qg.setpointers[i]
    lengthA = qg.setstorage[head]
    lengthE = qg.setstorage[head+1]
    firstkAE = head + 2
    lastkAE = head + 2 + lengthA + lengthE - 1
    return firstkAE:lastkAE
end
"`rangesAandE(qg, i)` returns the ranges for each of sets `A_i` and `E_i` in `qg.setstorage`"
@inline function rangesAandE(qg, i)
    head = qg.setpointers[i]
    lengthA = qg.setstorage[head]
    lengthE = qg.setstorage[head+1]
    firstkA = head + 2
    firstkE = head + 2 + lengthA
    lastkA = head + 2 + lengthA - 1
    lastkE = head + 2 + lengthA + lengthE - 1
    return firstkA:lastkA, firstkE:lastkE
end
"`rangeL(qg, i)` returns the range of set `L_i` in `qg.setstorage`"
@inline function rangeL(qg, e)
    head = qg.setpointers[e]
    lengthL = qg.setstorage[head]
    firstkL = head < qg.stackpointers[1] ? head+2 : head+1
    # FUTURE An ifelse() perform better than the ternary op
    lastkL = firstkL + lengthL -1
    return firstkL:lastkL
end


"""
    QuotientGraph{Tv,Ti<:Signed}(mat::SparseMatrixCSC{Tv,Ti})

Builds a `QuotientGraph{Ti}` representing `SparseMatrixCSC{Tv,Ti}` `mat`. This method
assumes that `mat` is symmetric and contains all entries rather than simple an upper or
lower triangle. `mat` need not contain numerical values; only `mat`'s column pointers and
row values need exist.
"""
function QuotientGraph{Tv,Ti<:Signed}(mat::SparseMatrixCSC{Tv,Ti})
    matorder = mat.n
    matcolptrs = mat.colptr
    matrowvals = mat.rowval
    matentryn = matcolptrs[matorder+1] - 1
    
    nodestatus = ones(Ti, matorder)
    setpointers = Vector{Ti}(matorder)
    setstorage = Vector{Ti}(matentryn + 2*matorder)
    
    # Iterate through the columns of the given matrix. From column j, generate initial sets
    # A_j and E_j for the corresponding supervariable j. Specifically, store a pointer into
    # qg.setstorage identifying the beginning of the definitions and contents of sets A_j
    # and E_j, i.e. setpointers[jcol] = headj. Store the number of off-diagonal entries
    # in the column, i.e. the length of set A_j, in setstorage[headj], and the length of
    # E_j (zero at the outset) in setstorage[headj+1]. Store the off-diagonal entries in
    # column j, i.e. the contents of A_j, contiguously in setstorage following the
    # preceding lengths.
    kwrite = 1
    for jcol in 1:matorder
        headj = kwrite
        setpointers[jcol] = headj
        coljfirstk = matcolptrs[jcol]
        coljlastk = matcolptrs[jcol+1] - 1
        coljlength = coljlastk - coljfirstk + 1
        setstorage[kwrite] = coljlength - 1; kwrite += 1
        setstorage[kwrite] = zero(Ti); kwrite += 1
        diagonalpresent = false
        for krow in coljfirstk:coljlastk
            irow = matrowvals[krow]
            if irow == jcol
                diagonalpresent = true
            else
                setstorage[kwrite] = irow
                kwrite += 1
            end
        end
        # Prior to reading column j, we assume that it contains a diagonal entry, on that
        # basis estimate the length of A_j as (coljlength - 1), and write that estimate
        # in the appropriate location. Most of the time this assumption (and so estimate)
        # is correct, in which case writing the estimate at the outset keeps memory access
        # in-order. If that assumption fails (and so the estimate is incorrect), we need
        # touch up the estimated length we wrote at the outset:
        diagonalpresent || (setstorage[headj] = coljlength)
    end
    
    QuotientGraph{Ti}(matorder, nodestatus, setstorage, setpointers, Ti[kwrite, kwrite])
end

"""
    elementize!{T}(qg::QuotientGraph{T}, nwpq, p, workspaceA, workspaceB)

Transform `QuotientGraph` `qg` such that 'pivot' supervariable `p` becomes an element.
`nwpq` must be the priority queue in which node weights are maintained during the
elmiination process; though hypothetically `nwpq` need only support `dequeue!(nwpq, j)`
where `j` is a supervariable in the node weight priority queue, at the moment only
`NWLLPriorityQueue` structures are supported. This method drops supervariables that were
amalgamated into other supervariables from the node weight priority queue, and also drops
those supervariables which require node weight updates. `workspaceA` and `workspaceB` must
be `Vector{T}`s accommodating the largest set involved in the transformation (`L_p`);
`QuotientGraph.matorder` is a safe length for these workspaces.
"""
function elementize!{T}(qg::QuotientGraph{T}, nodeweightpq, p, workspaceA, workspaceB)
    # Transmuting supervariable p into element p involves three gross stages.
    #
    # The first stage is building [the set of supervariables to which element p links, L_p],
    # from [the set of supervariables to which supervariable p linked, A_p], and, for each
    # element e in [the set of elements to which supervariable p linked, E_p], [the set of
    # supervariables to which element e linked, L_e]. Specifically, L_p = union(A_p, 
    # {L_e | e in E_p}) \ {p}. In this process, we also absorb all elements involved in
    # construction of the set L_p into element p, and clean up sets A_p and E_p.
    #
    # The second stage is, for each supervariable i in [the set of supervariables to which
    # element p links, L_p]:(1) updating [the set of supervariables to which supervariable i
    # links directly, A_i], by stripping from it the set of supervariables to which
    # supervariable i now links indirectly through element p and also p itself if p was in
    # A_i (altogether A_i -> setdiff(A_i, union(L_p, {p}))); and (2) updating [the set of
    # elements to which supervariable i links, E_i], by stripping from it the set of 
    # elements absorbed into element p, E_p], and injecting a link to element p (altogether
    # E_i -> union(setdiff(E_i, E_p), {p}).
    #
    # The third stage is indistinguishable supervariable detection and amalgamation. This
    # stage is in flux; see the documentation on this stage below.
    #
    # In some cases much of the work descirbed above is unnecessary a priori. Hence several
    # different execution paths exist below. More potential fast paths exist, but for now
    # we only treat some of those which are likely to occur and allow us to build and
    # store L_p in-place somehow.
    setstorage = qg.setstorage
    setpointers = qg.setpointers
    stackpointers = qg.stackpointers
    headp = setpointers[p]
    lengthAp = setstorage[headp]
    lengthEp = setstorage[headp+1]
    
    DEBUG && println("elementize! p $p lEp $lengthEp lAp $lengthAp contents $(setstorage[rangeAE(qg, p)])")    
    
    if lengthEp == 0
        # Fast paths where E_p is empty.
        if lengthAp > 0
            # Fast path where E_p is empty but A_p is nonempty.
            #
            # Concerning stage one:
            # As E_p is empty, L_p is simply A_p. Hence we leave A_p intact, in-place as L_p
            # and stage one becomes trivial: No work constructing L_p, no work storing L_p,
            # no work absorbing elements, no work cleaning up A_p and E_p.
            #
            # Concerning stage two:
            # In general, for each supervariable i in L_p we need strip E_p from E_i. But as
            # E_p is empty this becomes trivial. Additionally, given that E_p is empty we
            # know that supervariables p and i must have been directly linked. Hence we
            # know that p is in A_i. Work for stage two reduces to injecting p into E_p,
            # and stripping union(L_p, {p}) from A_i, which we can do in one pass.
            lengthLp = lengthAp
            firstkLp = headp + 2
            lastkLp = headp + 2 + lengthAp - 1
            for kLp in firstkLp:lastkLp
                iLp = setstorage[kLp]
                headi = setpointers[iLp]
                lengthAi = setstorage[headi]
                lengthEi = setstorage[headi+1]
                setstorage[headi+1] = lengthEi + 1
                # We know a priori that lengthEi must be one greater on exit than on entry.
                # See below for reasoning. Hence we set lengthEi to the expected value
                # immediately following first access in the interest of linearizing memory
                # access as much as possible.
                #
                # We need strip union(L_p, {p}) from A_i and then inject p into E_i.
                # We know that E_p is empty and A_p is nonempty. That A_i necessarily
                # contains p follows. Beyond that we have no additional information about
                # A_i and E_i in general.
                #
                # FUTURE Potential fast paths exist. Perhaps try these at some point.
                #
                # Where lengthEi == 0 and lengthAi >= 1, we can simplify injection of p
                # into E_i to two statements: (1) assign p to the first location following
                # A_i; and (2) set lengthEi = 1.
                #
                # Where lengthEi == 0 and lengthAi == 1, we know that A_i's sole entry must
                # be p. Hence we can simplify stripping of union(L_p, {p}) from A_i and
                # injecting p into E_i: (1) set lengthAi = 0; and (2) set lengthEi = 1.
                #
                # Do other cases exist admitting simple fast paths?
                #
                # General code:
                firstkAi = headi + 2
                firstkEi = headi + 2 + lengthAi
                # Strip union(L_p, {p}) from A_i.
                lengthAi = strip!(setstorage, firstkAi, lengthAi, firstkLp, lengthLp, p)
                # Given p necessarily exists in A_i, this operation must leave a gap between
                # the now-stripped A_i and the un-revised E_i. Hence we need rewrite E_i
                # after A_i and while injecting p.
                shiftinsert!(setstorage, firstkEi, lengthEi, firstkAi+lengthAi, p)
                # We do not require the return value of shiftinsert! (lengthEi) because we
                # know that lengthEi must be one greater than it was on entry. For the same
                # reason we touch up setstorage[headi+1] (stored lengthEi) above rather than
                # below.
                setstorage[headi] = lengthAi
                # setstorage[headi+1] = lengthEi
            end
        # elseif lengthAp == 1
            # FUTURE Fast path where E_p is empty and A_p contains a single entry.
        # elseif lengthAp == 0
        else
            # Fast path where both E_p and A_p are empty.
            #
            # If both E_p and A_p are empty, then L_p must be empty and p's elimination has
            # no impact on any other supervariable or element in the quotient graph. Hence,
            # in this case we can simply remove A_p and E_p from the quotient graph, which
            # we do by marking p as follows:
            setpointers[p] = 0
            # To be clear, we should be able to skip indistinguishable supervariable
            # detection and node weight updating in this case, as we have changed nothing
            # in our graph apart from removing this isolated node. Return'ing immediately
            # skips indistinguishable supervariable detection, but I am not certain whether
            # it skips everything that can be skipped after returning from elementize!.
            # In fact, I do not think this code has been exercised yet; hence chances
            # are it fails.
            # FUTURE Under what circumstances precisely can this case occur?
            # FUTURE Do we skip all skip-able work later in this elimination step?
            return
        end
    # elseif lengthEp == 1
    #     if lengthAp == 0
    elseif lengthEp == 1 && lengthAp == 0
        # Fast path where A_p is empty and E_p contains a single entry.
        #
        # I.e., supervariable p links solely to a single element in the quotient graph, say
        # q. The safe, consistent way of handlign this case is to create a new element p
        # and absorb element q into element p. We take this approach below, but we may be
        # abel to do better (see external notes).
        #
        # Concerning stage one:
        # As A_p is empty and E_p = {q}, L_p = setdiff(L_q, {p}). Hence we form L_p by 
        # stripping p from L_q, pointing L_p to L_q's storage, and indicating absorption of
        # q into p.
        q = setstorage[headp+2]
        headq = setpointers[q]
        setpointers[q] = -p # Absorb immediately after first access
        setpointers[p] = headq
        lengthLp = setstorage[headq] - 1
        setstorage[headq] = lengthLp
        firstkLp = headq < stackpointers[1] ? headq+2 : headq+1
        # FUTURE An ifelse() may perform better than this ternary
        strip!(setstorage, firstkLp, lengthLp+1, p)
        setstorage[headq] = lengthLp
        #
        # Concerning stage two:
        # We need only replace q with p in each E_i for i in L_p.
        lastkLp = firstkLp + lengthLp - 1
        for kLp in firstkLp:lastkLp
            iLp = setstorage[kLp]
            headi = setpointers[iLp]
            lengthAi = setstorage[headi]
            lengthEi = setstorage[headi+1]
            # q necessarily exists in E_i. Beyond that we know nothing about E_i in general
            # that we can exploit. FUTURE But if E_i contains a single entry, we know it
            # must be q. Hence to strip q from E_i and insert p, we need only overwrite
            # E_i's single entry q with p.
            #
            # General code
            firstkEi = headi + 2 + lengthAi
            stripinsert!(setstorage, firstkEi, lengthEi, p, q)
        end
    # elseif lengthEp == 1 && lengthAp == 1
        # FUTURE Fast path where E_p and A_p both contain one entry.
    else
        # Hereafter we have lengthEp >= 2 or (lengthEp == 1 && lengthAp >= 1), though as
        # described just above in the future the latter conditition may become
        # lengthEp == 1 && lengthAp >= 2.
        #
        # Concerning stage one:
        # We need build L_p = union(A_p, L_e for e in E_p). We evaluate the union via
        # a series of incremental merges presently; see FUTURE below regarding more
        # efficient approaches. If A_p is nonempty, we first construct workset =
        # union(A_p, L_e1) and then construct workset = union(workset, L_e) in series
        # for each e in setdiff(E_p, e_1) yielding union(A_p, L_e for e in E_p). If instead
        # A_p is empty, we first construct workset = union(L_e1, L_e2) and then proceed
        # in series with the remaining sets L_e as above.
        #
        # FUTURE Performing the union via a multi-way merge should be more efficient.
        firstkEp = headp + 2 + lengthAp
        lastkEp = headp + 2 + lengthAp + lengthEp - 1
        e1 = setstorage[firstkEp]
        headL1 = setpointers[e1]
        setpointers[e1] = -p # Absorb element e1 into element p
        lengthL1 = setstorage[headL1]
        firstkL1 = headL1 < stackpointers[1] ? headL1+2 : headL1+1
        if lengthAp == 0
            e2 = setstorage[firstkEp+1]
            headL2 = setpointers[e2]
            setpointers[e2] = -p # Absorb element e2 into element p
            lengthAporL2 = setstorage[headL2]
            firstkAporL2 = headL2 < stackpointers[1] ? headL2+2 : headL2+1
            nextkEp = firstkEp + 2
        else
            lengthAporL2 = lengthAp
            firstkAporL2 = headp + 2
            nextkEp = firstkEp + 1
        end
        lengthLp = union!(workspaceA,
            setstorage, firstkL1, lengthL1,
            setstorage, firstkAporL2, lengthAporL2)
        srcworkspace = workspaceA
        destworkspace = workspaceB
        for kEp in nextkEp:lastkEp
            e = setstorage[kEp]
            headLe = setpointers[e]
            setpointers[e] = -p # Absorb element e into element p
            lengthLe = setstorage[headLe]
            firstkLe = headLe < stackpointers[1] ? headLe+2 : headLe+1
            lengthLp = union!(destworkspace,
                srcworkspace, 1, lengthLp,
                setstorage, firstkLe, lengthLe)
            srcworkspace, destworkspace = destworkspace, srcworkspace
        end
        lengthLp -= 1    
        # srcworkspace should now contain union(L_p, {p}) from entry 1 through entry
        # (lengthLp + 1). To complete stage one, we need store L_p permanently. The existing
        # code defers storage till after stage two for two reasons: Storing L_p may require
        # garbage collection. We (1) may be able to collect a bit more garbage after
        # completing stage two; and (2) need E_p during stage two, and garbage collection
        # would destroy E_p. Chances are the first reason is weak: The amount of additional
        # garbage we would collect is likely marginal, hence the performane impact is
        # unlikely worth the code complexity this structure introduces. The second reason
        # is sound, though we could play some game like temporarily storing E_p in
        # workspace or on the L-stack. In any case, defer storage till after stage two now.
        # FUTURE
        #
        # Concerning stage two:
        # We iterate through L_p by iterating through union(L_p, {p}) in srcworkspace,
        # entries 1 through (lengthLp + 1), and skipping p when we encounter it.
        for kLp in 1:(lengthLp+1)
            iLp = srcworkspace[kLp]
            iLp == p && continue
            headi = setpointers[iLp]
            lengthAi = setstorage[headi]
            lengthEi = setstorage[headi+1]
            # In general we need strip union(L_p, {p}) from A_i, strip E_p from E_i, and
            # inject p into E_i. What do we know about the sets involved?
            #
            # We know that lengthEp >= 2 or (lengthEp == 1 && lengthAp >= 1), though as
            # described just above in the future the latter condition may become
            # lengthEp == 1 && lengthAp >= 2.
            #
            # We do not know anything about E_i and A_i in general it seems, other than that
            # they cannot be simultaneously empty, as otherwise i would not have ended up in
            # L_p. If A_p is empty, p must've linked to it through an element, in which
            # case p does not appear in A_i. Similarly, if A_i is empty, then p must've
            # linked to i through an element. If E_i is empty, p must've linked to i
            # directly.
            # FUTURE This and the different union construction paths suggest that at
            # some point separating the lengthAp == 0 case out may be worthwhile.
            #
            # Can we do less work in various cases?
            #
            if lengthAi == 0
                # If A_i is empty, then E_i is nonempty and must contain an entry in E_p. We
                # need only do the following: Strip E_p from E_i, and inject p. FUTURE If
                # E_p contains only a single entry, we may be able to reduce work by
                # calling a specialized method.
                firstkEi = headi + 2 + lengthAi
                lengthEi = stripinsert!(setstorage, firstkEi, lengthEi, firstkEp, lengthEp, p)
                setstorage[headi+1] = lengthEi
            elseif lengthEi == 0
                # If E_i is empty, then A_i is nonempty and must contain p itself. Here we
                # need only do the following: Strip union(L_p, {p}) from  A_i and inject
                # p into the empty E_i.
                firstkAi = headi + 2
                lengthAi = linearsetdiff!(setstorage, firstkAi, lengthAi, srcworkspace, 1, lengthLp+1)
                # FUTURE Determine safe assumptions and weaken linearsetdiff!
                firstkEi = headi + 2 + lengthAi
                setstorage[firstkEi] = p
                setstorage[headi] = lengthAi
                setstorage[headi+1] = one(T)
            else # lengthAi != 0 && lengthEi != 0
                # If neither E_i nor A_i are empty, we have the general case.    
                firstkAi = headi + 2
                firstkEi = headi + 2 + lengthAi
                # Strip union(L_p, {p}) from A_i
                olengthAi = lengthAi
                lengthAi = linearsetdiff!(setstorage, firstkAi, lengthAi, srcworkspace, 1, lengthLp+1)
                # FUTURE Determine safe assumptions and weaken linearsetdiff!
                if lengthAi == olengthAi
                    # We stripped nothing from A_i, hence we must modify E_i in place
                    lengthEi = stripinsert!(setstorage, firstkEi, lengthEi, firstkEp, lengthEp, p)
                else
                    # We stripped _something_ from A_i, hence we shift and modify E_i
                    # Shift E_i and simultaneously strip E_p and inject p
                    lengthEi = shiftstripinsert!(setstorage, firstkEi, lengthEi, firstkAi+lengthAi, firstkEp, lengthEp, p)
                end
                # Update E_i and A_i lengths
                setstorage[headi] = lengthAi
                setstorage[headi+1] = lengthEi
            end
        end
        #
        # Now we store L_p.
        #
        # FUTURE We have enough information to figure out whether storing L_p in-place is
        # possible, in a strong sense: We can find the last set stored in line preceding
        # the sets for p and also that following. So perhaps we can play some games to
        # avoid garbage collection. But garbage collection may not be a significant cost,
        # and scanning for the preceding and following sets may be a significant cost,
        # and is also complex, so for now we simply place all L-sets in the L-stack
        # at the end of setstorage and garbage collect more regularly.
        #
        # If there is adequate space for L_p on the L-stack, store L_p there. If not,
        # garbage collect setstorage and store L_p in the freed space.
        #
        # We need space for lengthLp+1 elements on the stack given we store the length of
        # L_p immediately preceding L_p itself.
        openstacklength = length(setstorage) - stackpointers[2] + 1
        # FUTURE Weaken this length
        lengthLp+1 > openstacklength && garbagecollect!(qg)
        # Write lengthLp at stackpointers[2] to begin L_p
        setstorage[stackpointers[2]] = lengthLp
        # Write L_p proper immediately thereafter
        copyexcept!(setstorage, stackpointers[2]+1, srcworkspace, 1, lengthLp+1, p)
        # Link L_p to former stack write index
        setpointers[p] = stackpointers[2]
        # Advance stack write index
        stackpointers[2] += lengthLp + 1    
    end
    
    # Detect and merge indistinguishable supervariables
    #
    # Here, among the supervariables impacted by the quotient graph transformations 
    # we've performed above, we check for newfound indistinguishability in the quotient
    # graph. If we find indistinguishable supervariables, we merge them into a single
    # supervariable, reducing subsequent work in manipulating the quotient graph.
    #
    # Checking for indistinguishability properly is expensive. Hence, here we check
    # indistinguishabiltiy in two stages. In the first stage we compare hash values
    # among the supervariables we are checking, where hash collisions necessarily occur
    # where supervariables are indistinguishable, and rarely otherwise. In the second
    # stage, if any collisions occurred, we check for indistinguishability proper
    # among the supervariables that experienced hash collisions.
    #
    isvnumpairs = 0 # number of hash-supervariable pairs on the ISVD detection stack
    # for iLp in QuotientGraphs.SetL(qg, p)
    for kLp in QuotientGraphs.rangeL(qg, p)
        iLp = qg.setstorage[kLp]
        # Calculate the hash value for this supervariable
        sumAEi = zero(T)
        for kAEi in QuotientGraphs.rangeAE(qg, iLp)
            jAEi = qg.setstorage[kAEi]
            sumAEi += jAEi
        end
        # Insert this supervariable-hash pair into the indistinguishable supervariable
        # detection stack for processing below.
        isvinsert!(kLp, iLp, sumAEi, workspaceA, isvnumpairs)
        isvnumpairs += 1
        # FUTURE Calculate the hashes while manipulating the sets in earlier parts of
        # quotient graph transformation and insert them then.
    end
    isvdetectmerge!(qg, nodeweightpq, workspaceA, isvnumpairs)
  
    return # kill return type instability warning
end

"""
    isvinsert!(kLp, iLp, hash, pairstack, numpairs)

Push (hash, iLp) onto the hash-supervariable pair stack for later indistinguishable
supervariable detection.
"""
function isvinsert!(kLp, iLp, hash, pairstack, numpairs)
    pairstack[2*(numpairs+1)-1] = hash
    pairstack[2*(numpairs+1)] = iLp
end

function radixsort_pairstack!() error("Pair-stack radixsort not implemented!") end
function quicksort_pairstack!() error("Pair-stack quicksort not implemented!") end
"Insertion sort the hash-supervariable pair stack for indistinguishable supervariable detection."
function insertionsort_pairstack!(pairstack, numpairs)
    # Manually pulling out common subexpressions, replacing multiplications with
    # additions and subtractions, etcetera, all either do not improve performance
    # or curiously worsen it. Avoid spending time on microoptimization.
    #
    # Initial implementation from Wikipedia pseudocode. If-continue-break performance
    # trick later from the Julia standard library's implementation.
    @inbounds for pairk in 2:numpairs
        hashk = pairstack[2*pairk-1]
        svark = pairstack[2*pairk]
        pairm = pairk - 1
        while pairm >= 1
            if pairstack[2*pairm-1] > hashk
                pairstack[2*(pairm+1)-1] = pairstack[2*pairm-1]
                pairstack[2*(pairm+1)] = pairstack[2*pairm]
                pairm -= 1
                continue
            end
            break
        end
        pairstack[2*(pairm+1)-1] = hashk
        pairstack[2*(pairm+1)] = svark
    end
    return
end
"""
    isvdetectmerge!(qg, nodeweightpq, pairstack, numpairs)
    
Detect indistinguishable supervariables among those on the hash-supervariable pair stack,
merging indistinguishable supervariables in the `QuotientGraph` `qg` and dropping the
amalgamated supervariable from the node weight priority queue `nodeweightpq`.
"""
function isvdetectmerge!(qg, nodeweightpq, pairstack, numpairs)
    numpairs < 2 && return
    insertionsort_pairstack!(pairstack, numpairs)

    basek = 1
    while basek < numpairs
        # Select base hash-supervariable pair for comparisons
        basehash = pairstack[2*basek-1]
        basesv = pairstack[2*basek]
        basehash == 0 && (basek += 1; continue)
        # Select first point hash-supervariable pair for comparisons
        pointk = basek + 1
        pointhash = pairstack[2*pointk-1]
        pointsv = pairstack[2*pointk]
        # Iterate through remaining hash-supervariable pairs as point
        while pointhash <= basehash
            # If point and base hashes are equal, then check for indistinguishability and
            # advance point. If they are not equal, point hash must be less than base hash,
            # which can only obtain if point hash is zero, in which case point was
            # previously merged and we simply advance point.
            if pointhash == basehash
                # Check base and point for indistinguishability
                if QuotientGraphs.indistinguishable(qg, basesv, pointsv)
                    # Retrieve point supervariable's cardinality prior to removing
                    # it from the quotient graph, given we need it below
                    pointsvcard = QuotientGraphs.numsimplevars(qg, pointsv)
                    # Merge the point supervariable into the base supervariable
                    QuotientGraphs.merge!(qg, basesv, pointsv)
                    # Remove point supervariable from the node weight priority queue
                    NWLLPriorityQueues.dequeue!(nodeweightpq, pointsv)
                    # Update base supervariable's approx. external degree (aed) estimate
                    oldaedest = NWLLPriorityQueues.getnodeweight(nodeweightpq, basesv)
                    newaedest = oldaedest - pointsvcard
                    NWLLPriorityQueues.update!(nodeweightpq, basesv, newaedest)
                    # Remove point supervariable from the stack
                    pairstack[2*pointk-1] = 0
                end
            end
            pointk == numpairs && break
            pointk += 1
            pointhash = pairstack[2*pointk-1]
            pointsv = pairstack[2*pointk]
        end
        basek += 1
    end
    return
end

"""
    indistinguishable{T}(qg::QuotientGraph{T}, i::T, j::T)

Check whether supervariables `i` and `j` are indistinguishable in the `QuotientGraph` `qg`.
Returns `true` if so, `false` otherwise.
"""
function indistinguishable{T}(qg::QuotientGraph{T}, i::T, j::T)
    DEBUG && println("  Checking indistinguishability of i $i and j $j")
    # Grossly we need ascertain E_i == E_j and union(A_i, i) == union(A_j, j).
    #
    # To ascertain E_i == E_j, we first check the necessary condition
    # length(E_i) == length(E_j). If that inexpensive check succeeds, we then iterate
    # through E_i and E_j simultaneously, checking for elementwise equality.
    #
    # To ascertain union(A_i, i) == union(A_j, j), grossly we perform the same checks:
    # We first check the necessary condition length(union(A_i, i)) == length(union(A_j, j)).
    # If that inexpensive check succeeds, we then check union(A_i, i) and union(A_j, j) for
    # elementwise equality. Here, however, we decompose each of those two checks into
    # simpler separate checks:
    #
    # Noting that A_i and i are disjoint, and similarly A_j and disjoint, the check
    # length(union(A_i, i)) == length(union(A_j, j)) reduces to length(A_i) == length(A_j).
    #
    # The elementwise equality check requires i in A_j, j in A_i, and, having checked
    # length(A_i) == length(A_j) above, that A_i is a subset of union(A_j, j). Note,
    # however, that if our quotient graph is valid, i in A_j implies j in A_i and vice
    # versa. Hence checking the former or latter suffices. Furthermore, note that if
    # E_i and E_j are nonempty and elementwise equal in the quotient graph, i and j
    # must belong to at least one common clique in the elimination graph, in which case
    # we need not check i in A_j or j in A_i, as we already know they must be linked in
    # the elimination graph.
    #
    # We perform the checks in order of increasing expense:
    #
    setpointers = qg.setpointers
    setstorage = qg.setstorage
    headi = setpointers[i]
    headj = setpointers[j]
    # Check length(A_i) == length(A_j)
    lengthAi = setstorage[headi]
    lengthAj = setstorage[headj]
    lengthAi == lengthAj || return false
    DEBUG && println("      Passed A_i and A_j length check")
    # Check length(E_i) == length(E_j)
    lengthEi = setstorage[headi+1]
    lengthEj = setstorage[headj+1]
    lengthEi == lengthEj || return false
    DEBUG && println("      Passed E_i and E_j length check")
    # Check elementwise equality of A_i and A_j, parts fused:
    iinAj = false
    jinAi = false
    firstkAi = headi + 2
    firstkAj = headj + 2
    lastkAi = firstkAi + lengthAi - 1
    lastkAj = firstkAj + lengthAj - 1
    kAi = firstkAi
    kAj = firstkAj
    while true
        kAi > lastkAi && break
        kAj > lastkAj && break
        valAi = setstorage[kAi]
        valAi == j && (jinAi = true; kAi += 1; continue)
        valAj = setstorage[kAj]
        valAj == i && (iinAj = true; kAj += 1; continue)
        valAi == valAj || return false
        kAi += 1; kAj += 1
    end
    DEBUG && println("      Passed A_i and A_j equality check")
    # Arriving at this point ascertains setdiff(A_i, j) = setdiff(A_j, i). Last part:
    (iinAj && jinAi) || (lengthEi != 0) || return false
    # Check elementwise equality of E_i and E_j
    firstkEi = headi + 2 + lengthAi
    firstkEj = headj + 2 + lengthAj
    lastkEi = firstkEi + lengthEi - 1
    lastkEj = firstkEj + lengthEj - 1
    rangeEi = firstkEi:lastkEi
    rangeEj = firstkEj:lastkEj
    for (kEi, kEj) in zip(rangeEi, rangeEj)
        setstorage[kEi] == setstorage[kEj] || return false
    end
    DEBUG && println("      Passed E_i and E_j equality check")
    return true
end

"""
    merge!{T}(qg::QuotientGraph{T}, i::T, j::T)

Transform `QuotientGraph` `qg` such that supervariable `j` merges into supervariable `i`.
"""
function merge!{T<:Integer}(qg::QuotientGraph{T}, i::T, j::T)
    setpointers = qg.setpointers
    setstorage = qg.setstorage
    stackpointers = qg.stackpointers
    headj = setpointers[j]
    setpointers[j] = -i # mark j absorbed into i
    lengthAj = setstorage[headj]
    lengthEj = setstorage[headj+1]    
    # Then we remove all (direct) links from other supervariables to supervariable j
    firstkAj = headj + 2
    lastkAj = firstkAj + lengthAj - 1
    for kAj in firstkAj:lastkAj
        a = setstorage[kAj]
        heada = setpointers[a]
        lengthAa = setstorage[heada] - 1 # correct preemptively
        setstorage[heada] = lengthAa # write correction preemptively
        lengthEa = setstorage[heada+1]
        firstkAa = heada + 2
        firstkEa = heada + 2 + lengthAa # preemptively corrected
        lengthAa = strip!(setstorage, firstkAa, lengthAa+1, j)
        copy!(setstorage, firstkEa, setstorage, firstkEa+1, lengthEa)
        # FUTURE Weaken this copy!
    end
    # Then we remove all links from elements to supervariable j
    firstkEj = headj + 2 + lengthAj
    lastkEj = headj + 2 + lengthAj + lengthEj - 1
    for kEj in firstkEj:lastkEj
        e = setstorage[kEj]
        headLe = setpointers[e]
        lengthLe = setstorage[headLe] - 1 # correct preemptively
        setstorage[headLe] = lengthLe # write correction preemptively
        firstkLe = headLe < stackpointers[1] ? headLe+2 : headLe+1
        strip!(setstorage, firstkLe, lengthLe+1, j)
    end
    # Finally we increase supervariable i's simple-variable count by supervariable j's
    # and mark supervariable j as amalgamated prior to elimination.
    qg.nodestatus[i] = qg.nodestatus[i] + qg.nodestatus[j]
    qg.nodestatus[j] = 0
end

"""
    garbagecollect!{T<:Integer}(qg::QuotientGraph{T})

Garbage collect `QuotientGraph` `qg`'s set storage, freeing up space for new element (L) sets.
"""
function garbagecollect!{T<:Integer}(qg::QuotientGraph{T})
    qgmatorder = qg.matorder
    setstorage = qg.setstorage
    setpointers = qg.setpointers
    stackpointers = qg.stackpointers
    # We garbage collect qg.setstorage in two stages.
    #
    # First we sweep through setpointers. For each live set `i` encountered in setpointers,
    # we retrieve the first entry (either length(A_i) or length(L_i)) from the corresponding
    # section of setstorage, store that entry in the corresponding position in setpointers,
    # and that entry in setstorage with `-i`, indicating the beginning of a live set
    # region for the second stage.
    #
    # Second we sweep through setstorage. For each live set `i` encountered as marked in the
    # preceding stage, we retrieve its first entry from the corresponding location in
    # setpointers, update that setpointer such that it points where we will rewrite the
    # set in setstorage, and then rewrite the first and remaining entries of that set to
    # their new (left-shifted) position in setstorage.
    #
    # Finally we update the stackpointers.
    #
    # First stage, marking active set regions.
    for i in 1:length(setpointers)
        pointeri = setpointers[i]
        if pointeri > 0
            setpointers[i] = setstorage[pointeri] # store first entry in setpointers
            setstorage[pointeri] = -i # mark over first entry
        end
    end
    # Second stage, rewriting, part 1: A/E and inline-L sets
    kread = 1
    kwrite = 1
    while kread < stackpointers[1]
        if setstorage[kread] < 0
            # Beginning of active set region
            seti = -setstorage[kread]; kread += 1
            lengthEi = setstorage[kread]; kread += 1
            lengthAi = setpointers[seti]
            setpointers[seti] = kwrite
            setstorage[kwrite] = lengthAi; kwrite += 1
            setstorage[kwrite] = lengthEi; kwrite += 1
            copy!(setstorage, kwrite, setstorage, kread, lengthAi+lengthEi)
            kwrite += lengthAi + lengthEi
            kread += lengthAi + lengthEi
        else
            kread += 1
        end
    end
    stackpointers[1] = kwrite
    # Second stage, rewriting, part 2: L sets
    while kread < stackpointers[2]
        if setstorage[kread] < 0
            seti = -setstorage[kread]; kread += 1
            lengthLi = setpointers[seti]
            setpointers[seti] = kwrite
            setstorage[kwrite] = lengthLi; kwrite += 1
            copy!(setstorage, kwrite, setstorage, kread, lengthLi)
            kwrite += lengthLi
            kread += lengthLi
        else
            kread += 1
        end
    end
    stackpointers[2] = kwrite
end

"""
    union!(deststorage, setAstorage, setAfirstk, setAlength, setBstorage, setBfirstk, setBlength)

Constructs the union of the sorted vector sets `A` (in `setAstorage`, beginning at
`setAfirstk` and of length `setAlength`) and `B` (analogously defined) and stores that
union in `deststorage`, beginning from the first position in `deststorage`. Returns the
length of the union.
"""
function union!(deststorage,
        setAstorage, setAfirstk, setAlength,
        setBstorage, setBfirstk, setBlength)
    # If both source sets are empty, return an empty set
    setAlength == 0 && setBlength == 0 && return 0
    # If source set A is empty, but source set B is nonempty, write B through
    if setAlength == 0 && setBlength > 0
        copy!(deststorage, 1, setBstorage, setBfirstk, setBlength)
        return setBlength
    end
    # If source set A is nonempty, but source set B is empty, write A through
    if setAlength > 0 && setBlength == 0
        copy!(deststorage, 1, setAstorage, setAfirstk, setAlength)
        return setAlength
    end
    # FUTURE Compactify / clean up control flow above
    # To reach this point, source sets A and B must be nonempty. Merge them properly.
    kA = setAfirstk; endA = setAfirstk + setAlength - 1
    kB = setBfirstk; endB = setBfirstk + setBlength - 1
    kW = 1
    while true
        if setAstorage[kA] < setBstorage[kB]
            deststorage[kW] = setAstorage[kA]
            kA += 1; kW += 1
            if kA > endA # We've exhausted A, but not B
                tailBlength = endB - kB + 1
                copy!(deststorage, kW, setBstorage, kB, tailBlength)
                return kW - 1 + tailBlength
            end
        elseif setAstorage[kA] > setBstorage[kB]
            deststorage[kW] = setBstorage[kB]
            kB += 1; kW += 1
            if kB > endB # We've exhausted B, but not A
                tailAlength = endA - kA + 1
                copy!(deststorage, kW, setAstorage, kA, tailAlength)
                return kW - 1 + tailAlength
            end
        else # setAstorage[kA] == setBstorage[kB]
            deststorage[kW] = setAstorage[kA]
            kA += 1; kB +=1; kW += 1
            if kA <= endA && kB <= endB # We've exhausted neither A nor B
                continue
            end
            if kA > endA && kB <= endB # We've exhausted A, but not B
                tailBlength = endB - kB + 1
                copy!(deststorage, kW, setBstorage, kB, tailBlength)
                return kW - 1 + tailBlength
            end
            if kA <= endA && kB > endB # We've exhausted B, but not A
                tailAlength = endA - kA + 1
                copy!(deststorage, kW, setAstorage, kA, tailAlength)
                return kW - 1 + tailAlength
            end
            if kA > endA && kB > endB # We've exhausted both A and B
                return kW - 1
            end
        end
    end
    # FUTURE Compactify / clean up control flow above
    return # kill type instability warning
end

# FUTURE Rename extendinsert!
"""
    rightshiftinsert!(setstorage, setfirstk, setlength, insertel)

Inserts entry `insertel` into the sorted-vector set beginning at `setfirstk` in `setstorage`
of length `setlength`. Specifically, if the set is empty, this method places the entry in
`setstorage` at `setfirstk`. If the set is not empty, this method shifts all entries in
the set greater than `insertel` one position right, and places `insertel` in the gap
created. This method performs no checking of any form, and assumes `insertel` does not
exist in the original set.
"""
function rightshiftinsert!(setstorage, setfirstk, setlength, insertel)
    setlastk = setfirstk + setlength - 1
    kwrite = setlastk + 1
    kread = setlastk
    
    while true
        if kread < setfirstk # We've exhausted the set
            setstorage[kwrite] = insertel
            break
        end
        if insertel > setstorage[kread]
            setstorage[kwrite] = insertel
            break
        end
        if insertel < setstorage[kread]
            setstorage[kwrite] = setstorage[kread]
            kwrite -= 1
            kread -= 1
            continue
        # elseif insertel == setstorage[kread]
        #     error("equality obtained in rightshiftinsert! insertel $insertel")
        end
    end
        
    return setlength + 1
end

"""
    shiftinsert!(setstorage, firstkS, lengthS, destk, insertval)

Writes set `union(S, {x})` at location `destk` in `setstorage`. Assumes set `S` does not
already contain `x`. Also assumes location `destk` falls to the left of `firstkS`, i.e.
`destk < firstkS`.
"""
function shiftinsert!{T<:Signed}(setstorage::Vector{T}, firstkS::T, lengthS::T, destk::T, insertval::T)
    lastkS = firstkS + lengthS - 1
    kread = firstkS
    kwrite = destk
    # Scan through S till we either exhaust S or encounter an entry greater than insertval,
    # rewriting S beginning at destk as we go.
    while kread <= lastkS && setstorage[kread] < insertval
        setstorage[kwrite] = setstorage[kread]
        kread += 1
        kwrite += 1
    end
    # Either we exhausted S prior to encountering an entry greater than insertval, either
    # because S was empty or no entry greater than p exists in S, or we encountered an
    # entry in S greater than insertval prior to exhausting S. In either case, our next
    # move is to write inserval.
    setstorage[kwrite] = insertval
    kwrite += 1
    # If we did not exhaust S, we need rewrite its remainder. In anticipation
    # of adding summing the entries of S, we ignore special cases that would allow us to
    # avoid this rewrite.
    while kread <= lastkS
        setstorage[kwrite] = setstorage[kread]
        kread += 1
        kwrite += 1
    end
    return lengthS+1
    # FUTURE return value known a priori, shouldn't really compute the increment and return.
    # Instead call sites should be intelligent.
end

function strip!(setstorage, setAfirstk, setAlength, setBfirstk, setBlength, x)
    linearsetdiff!(setstorage, setAfirstk, setAlength, setstorage, setBfirstk, setBlength, x)
end
"""
    linearsetdiff!(setAstorage, setAfirstk, setAlength, setBstorage, setBfirstk, setBlength, x)

Assumes `x` is not in set `B`.
FUTURE Document.
"""
function linearsetdiff!(setAstorage, setAfirstk, setAlength, setBstorage, setBfirstk, setBlength, x)
    # If set A is empty, leave it be
    setAlength == 0 && return setAlength
    # If set B is empty, weaken this to a single-entry setdiff! for x
    setBlength == 0 && return strip!(setAstorage, setAfirstk, setAlength, x)
    #
    # FUTURE Consider grossly simplifying this at expense of some additional ops/branches.
    #
    # This operation consists of two coarse stages:
    #
    # The first coarse stage (1) determines whether there are any entries from union(B, {x})
    # in A, if so drops the least such entry, and, if additional such entries may exist,
    # directs control flow such that the second coarse stage looks for them in the correct
    # parts of union(B, {x}).
    #
    # The second coarse stage (2) writes setdiff(A, union(B, {x})) over A, beginning from
    # the position where the first stage found the first drop-entry.
    #
    # These two coarse stages each consist of three fine stages.
    #
    # The three fine stages of coarse stage one are (1-1), (1-2), and (1-3). Each such fine
    # stage (1-?) either returns, passes control to the next such fine stage (1-?), or
    # skips the other such fine stages (1-?) and passes control to a fine stage in the
    # second coarse stage (2-?). Hence (1-3) only executes after (1-2), and (1-2) only after
    # (1-1).
    #
    # (1-1) Scans for existence of [entries from B less than x] in A. If no such entries
    #       exist in A and no other entries from union(B, {x}) may exist in A, (1-1)
    #       returns. If no such entries exist in A but other entries from union(B, {x}) may
    #       may exist in A, (1-1) passes control to (1-2). If such an entry exists and
    #       additional such entries may exist in A, (1-1) drops the entry and passes
    #       control to (2-1). If such an entry exists, no additional such entries may
    #       exist in A, but other entries from union(B, {x}) may exist in A, (1-1)
    #       drops the entry and passes control to (2-2). If such an entry exists and no
    #       other entries from union(B, {x}) may exist in A, (1-1) returns.
    # (1-2) Scans for existence of x in A. If x does not exist in A and no other entries
    #       from union(B, {x}) may exist in A, (1-2) returns. If x does not exist in A but
    #       [entries from B greater than x] may exist in A, (1-2) passes control
    #       to (1-3). If x exists in A but no other entries from union(B, {x}) may exist in
    #       A, (1-2) drops x and returns. If x exists in A and [entries from B greater than
    #       x] may exist in A, (1-2) drops x and passes control to (2-3).
    # (1-3) Scans for existence of [entries from B greater than x] in A. If no such entries
    #       exist in A, (1-3) returns. If such an entry exists and no additional such
    #       entries may exist in A, (1-3) drops the entry and returns. If such an entry
    #       exists and additional such entries may exist in A, (1-3) drops the entry and
    #       passes control to (2-3).
    #
    # The three fine stages of coarse stage two are (2-1), (2-2), and (2-3). Each such fine
    # stage (2-?) either returns or passes control to the next such fine stage (2-?). The
    # action of fine stages (2-?) correspond roughly correspond to those in (1-?).
    #
    # (2-1) Completes setdiff(A, [entries from B less than x]). If no other entries from
    #       union(B, {x}) may exist in A, (2-1) then returns. Otherwise (2-1) passes control
    #       to (2-2).
    # (2-2) Completes setdiff(A, x). If no [entries from B greater than x] may exist in A,
    #       (2-2) then returns. Otherwise (2-2) passes control to (2-3).
    # (2-3) Completes setdiff(A, [entries from B greater than x]) and returns.
    #
    # Stage 1.
    jumpstage = 0
    # jumpstage functions as a pseudo-goto, identifying the logic that we need jump to next.
    # FUTURE this might be one of those rare cases where real gotos make sense.
    # jumpstage == 0 :: no entry in union(B, {x}) found in A thus far, hence we need proceed
    #                   to the next scanning logic stage (1-?)
    # jumpstage == 1 :: found and dropped [an entry in B less than x] in A, and without
    #                   exhausting B, hence we need proceed to rewriting logic stage (2-1),
    #                   i.e. where additional entries in B less than x may exist in A
    # jumpstage == 2 :: found and dropped [an entry in B less than x] in A, but having
    #                   exhausted B, hence we need proceed to rewriting logic stage (2-2),
    #                   i.e. where no additional entries in B less than x may exist in A
    # jumpstage == 3 :: found and dropped x in A, and without exhausting B, hence we need
    #                   proceed to rewriting logic stage (3-1), i.e. where additional
    #                   entries in B greater than x may exist in A
    kwA = setAfirstk
    krA = setAfirstk
    krB = setBfirstk
    endA = setAfirstk + setAlength - 1
    endB = setBfirstk + setBlength - 1
    setAentry = setAstorage[krA]
    setBentry = setBstorage[krB]
    #
    # Stage (1-1)
    setBentry < x && while true
        if setAentry < setBentry
            krA += 1
            krA > endA && return setAlength # done
            setAentry = setAstorage[krA]
            continue # -> (1-1)
        elseif setAentry > setBentry
            krB += 1
            krB > endB && break # -> (1-2)
            setBentry = setBstorage[krB]
            setBentry >= x && break # -> (1-2)
            continue # -> (1-1)
        else # setAentry == setBentry
            kwA = krA
            krA += 1
            krA > endA && return setAlength - 1 # done
            setAentry = setAstorage[krA]
            krB += 1
            if krB > endB
                jumpstage = 2 # -> (2-2)
            else
                setBentry = setBstorage[krB]
                jumpstage = setBentry < x ? 1 : 2 # -> (2-1), or (2-2)
            end
            break
        end
        # FUTURE Could accelerate the first two cases with inner whiles?
    end
    # Stage (1-2)
    jumpstage == 0 && while true
        if setAentry < x
            krA += 1
            krA > endA && return setAlength # done
            setAentry = setAstorage[krA]
            continue # -> (1-2)
        elseif setAentry > x
            krB > endB && return setAlength # done
            break # -> (1-3)
        else # setAentry == x
            kwA = krA
            krA += 1
            krA > endA && return setAlength - 1 # done
            setAentry = setAstorage[krA]
            if krB > endB
                tailAlength = endA - kwA
                # FUTURE Weaken this copy!
                copy!(setAstorage, kwA, setAstorage, krA, tailAlength)
                return setAlength - 1 # done
            end
            jumpstage = 3 # -> (2-3)
            break
        end
        # FUTURE Could accelerate the first two cases with inner whiles?
    end
    # Stage (1-3)
    jumpstage == 0 && while true
        if setAentry < setBentry
            krA += 1
            krA > endA && return setAlength # done
            setAentry = setAstorage[krA]
            continue # -> (1-3)
        elseif setAentry > setBentry
            krB += 1
            krB > endB && return setAlength # done
            setBentry = setBstorage[krB]
            continue # -> (1-3)
        else # setAentry == setBentry
            kwA = krA
            krA += 1
            krA > endA && return setAlength - 1 # done
            setAentry = setAstorage[krA]
            krB += 1
            if krB > endB
                tailAlength = endA - kwA
                # FUTURE Weaken this copy!
                copy!(setAstorage, kwA, setAstorage, krA, tailAlength)
                return setAlength - 1 # done
            end
            setBentry = setBstorage[krB]
            jumpstage = 3 # -> (2-3)
            break
        end
    end
    # Stage (2-1)
    jumpstage == 1 && while true
        if setAentry < setBentry
            setAstorage[kwA] = setAentry
            kwA += 1
            krA += 1
            krA > endA && return kwA - setAfirstk # done
            setAentry = setAstorage[krA]
            continue # -> (2-1)
        elseif setAentry > setBentry
            krB += 1
            krB > endB && break # -> (2-2)
            setBentry = setBstorage[krB]
            setBentry >= x && break # -> (2-2)
            continue # -> (2-1)
        else # setAentry == setBentry
            krA += 1
            krA > endA && return kwA - setAfirstk # done
            setAentry = setAstorage[krA]
            krB += 1
            krB > endB && break # -> (2-2)
            setBentry = setBstorage[krB]
            continue # -> (2-1)
        end
        # FUTURE Could accelerate the first two cases with inner whiles ?
    end
    # Stage (2-2)
    jumpstage <= 2 && while true
        if setAentry < x
            setAstorage[kwA] = setAentry
            kwA += 1
            krA += 1
            krA > endA && return kwA - setAfirstk # done
            setAentry = setAstorage[krA]
            continue # -> (2-2)
        elseif setAentry > x
            krB > endB && return kwA - setAfirstk # done
            break # -> (2-3)
        else # setAentry == x
            krA += 1
            krA > endA && return kwA - setAfirstk # done
            setAentry = setAstorage[krA]
            if krB > endB
                tailAlength = endA - krA + 1
                # FUTURE Weaken this copy!
                copy!(setAstorage, kwA, setAstorage, krA, tailAlength)
                return kwA - setAfirstk + tailAlength # done
            end
            break # -> (2-3)
        end
        # FUTURE Could accelerate the first two cases with inner whiles ?
    end
    # Stage (2-3)
    while true # jumpstage should be <= 3 here, so no need to check
        if setAentry < setBentry
            setAstorage[kwA] = setAentry
            kwA += 1
            krA += 1
            krA > endA && return kwA - setAfirstk # done
            setAentry = setAstorage[krA]
            continue # -> (2-3)
        elseif setAentry > setBentry
            krB += 1
            krB > endB && break # -> finish
            setBentry = setBstorage[krB]
            continue # -> (2-3)
        else # setAentry == setBentry
            krA += 1
            krA > endA && return kwA - setAfirstk # done
            setAentry = setAstorage[krA]
            krB += 1
            krB > endB && break # -> finish
            setBentry = setBstorage[krB]
            continue # -> (2-3)
        end
        # FUTURE Could accelerate the first two cases with inner whiles ?
    end
    # krB > endB, finish
    tailAlength = endA - krA + 1
    # FUTURE Weaken this copy!
    copy!(setAstorage, kwA, setAstorage, krA, tailAlength)
    return kwA - setAfirstk + tailAlength
end
"""
    linearsetdiff!(setAstorage, setAfirstk, setAlength, setBstorage, setBfirstk, setBlength)

Removes entries in sorted-vector set `B` (stored in `setBstorage` beginning at location
`setBfirstk` and of length `setBlength`) from set `A` (analogously defined), and returns
the length of the modified set `A`.
"""
function linearsetdiff!(setAstorage, setAfirstk, setAlength, setBstorage, setBfirstk, setBlength)
    # If set A is empty, leave it be
    setAlength == 0 && return setAlength
    # If set B is empty, leave A as is
    setBlength == 0 && return setAlength
    #
    # FUTURE Consider grossly simplifying this at expense of some additional ops/branches.
    #
    # This operation consists of two stages. The first stage determines whether there are
    # any entries from B in A, if so drops the least such entry, and, if additional entries
    # from B may exist in A, passes control to the second stage. The second stage writes
    # setdiff(A, B) over A, beginning from the position where the first stage found the
    # first drop-entry.
    kwA = setAfirstk
    krA = setAfirstk
    krB = setBfirstk
    endA = setAfirstk + setAlength - 1
    endB = setBfirstk + setBlength - 1
    setAentry = setAstorage[krA]
    setBentry = setBstorage[krB]
    #
    # Stage 1
    while true
        if setAentry < setBentry
            krA += 1
            krA > endA && return setAlength # done
            setAentry = setAstorage[krA]
            continue # -> S1
        elseif setAentry > setBentry
            krB += 1
            krB > endB && return setAlength # done
            setBentry = setBstorage[krB]
            continue # -> S1
        else # setAentry == setBentry
            kwA = krA
            krA += 1
            krA > endA && return setAlength - 1 # done
            setAentry = setAstorage[krA]
            krB += 1
            if krB > endB
                tailAlength = endA - kwA
                # FUTURE Weaken this copy!
                copy!(setAstorage, kwA, setAstorage, krA, tailAlength)
                return setAlength - 1 # done
            end
            setBentry = setBstorage[krB]
            break # -> S2
        end
        # FUTURE Could accelerate the first two cases with inner whiles ?
    end
    # Stage 2
    while true
        if setAentry < setBentry
            setAstorage[kwA] = setAentry
            kwA += 1
            krA += 1
            krA > endA && return kwA - setAfirstk # done
            setAentry = setAstorage[krA]
            continue # -> S2
        elseif setAentry > setBentry
            krB += 1
            krB > endB && break # -> finish
            setBentry = setBstorage[krB]
            continue # -> S2
        else # setAentry == setBentry
            krA += 1
            krA > endA && return kwA - setAfirstk # done
            setAentry = setAstorage[krA]
            krB += 1
            krB > endB && break # -> finish
            setBentry = setBstorage[krB]
            continue # -> S2
        end
        # FUTURE Could accelerate the first two cases with inner whiles ?
    end
    # krB > endB, finish
    tailAlength = endA - krA + 1
    # FUTURE Weaken this copy!
    copy!(setAstorage, kwA, setAstorage, krA, tailAlength)
    return kwA - setAfirstk + tailAlength
end
"""
    strip!(setstorage, firstkS, lengthS, x)
    
Strips `x` from set `S`. Assumes `S` originally contains `x`.
"""
function strip!(setstorage, firstkS, lengthS, x)
    lastkS = firstkS + lengthS - 1
    kread = firstkS
    # Scan through S till we encounter x
    while setstorage[kread] != x
        kread += 1
    end
    # Skip x, rewrite remainder of S left-shifted
    kwrite = kread
    kread += 1
    while kread <= lastkS
        setstorage[kwrite] = setstorage[kread]
        kwrite += 1
        kread += 1
    end
    return lengthS-1
    # FUTURE Nix the return, unnecessary
end
"""
    stripinsert!(setstorage, firstkS, lengthS, x, y)
    
Strips `x` from set `S` and inserts `y`. Assumes `S` originally contains `x` and does not
contain `y`. Also assumes `x != y`.
"""
function stripinsert!(setstorage, firstkS, lengthS, x, y)
    lastkS = firstkS + lengthS - 1
    kread = firstkS
    if x < y
        # If x is less than y, scan S till we hit x, drop x, proceed rewriting the following
        # entries in S left-shifted till we encounter an entry greater than y or exhaust S,
        # and then write y.
        while setstorage[kread] != x
            kread += 1
        end
        kwrite = kread
        kread += 1
        while kread <= lastkS && setstorage[kread] < y
            setstorage[kwrite] = setstorage[kread]
            kwrite += 1
            kread += 1
        end
        setstorage[kwrite] = y
    else # y < x
        # If y is less than x, scan S till we encounter an entry in S greater than y, cache
        # that entry, write y, and then proceed rewriting the following entries in S
        # (by way of the cache-entry) left-shifted till we encounter x, and drop x.
        while setstorage[kread] < y
            kread += 1
        end
        cachedent = setstorage[kread]
        setstorage[kread] = y
        kread += 1
        while setstorage[kread] != x
            temp = setstorage[kread]
            setstorage[kread] = cachedent
            cachedent = temp
            kread += 1
        end
        setstorage[kread] = cachedent
        kread += 1        
    end
    # Given we've dropped one entry and inserted another, we need not rewrite the tail.
    # Eventually need sum the tail though.
    return
end
"""
    stripinsert!(setstorage, firstkA, lengthA, firstkB, lengthB, x)
    
The logic in this case is much more complex than that where we shift `A` by at least one
location to the left. Hence the separation between this case and `shiftstripinsert!`
Strips set `B` from set `A` and injects `x`. Assumes that `x` was not already in set `A`.
"""
function stripinsert!(setstorage, firstkA, lengthA, firstkB, lengthB, x)
    # If set A is empty, inject x and return
    lengthA == 0 && (setstorage[destk] = x; return one(eltype(setstorage)))
    # If set B is empty, weaken this to a rightshiftinsert! for x
    lengthB == 0 && return rightshiftinsert!(setstorage, firstkA, lengthA, x)
    #
    # This operation is rather complex, so we break it into manageable pieces as follows:
    # To begin, we scan through A till we either (1) exhaust A, (2) exhuast B, (3) find an
    # entry from B in A, or (4) find an entry in A greater than x. If we (1) exhaust A, we
    # write x at the end of A and finish. If we (2) exhaust B, completing the operation
    # reduces to calling rightshiftinsert!(setstorage, firstkremA, lengthremA, x). If we
    # (3) find an entry from B in A, completing the operation reduces to dropping that
    # entry from A and calling shiftstripinsert!(setstorage, firstkremA, lengthremA,
    # kwA, firstkremB, lengthremB, x) where kwA points to the dropped entry's position in
    # A. If we (4) find an entry in A greater than x, see below.    
    #
    # Initial scan
    krA = firstkA
    krB = firstkB
    endA = firstkA + lengthA - 1
    endB = firstkB + lengthB - 1
    entryA = setstorage[krA]
    entryB = setstorage[krB]
    # FUTURE Could break this scan into two parts, avoiding some branches by comparing
    # x with entryB on entryB advancement rather than entryA with both x and entryB.
    # Wait on this, see whether it matters before introducing additional complexity.
    # We do this sort of splitting in the complex linearsetdiff! routine for stripping
    # both a set and an auxiliary entry, see that code for an outline if need be.
    while true
        if entryA > x
            # Jump to code below handling this case
            break
        elseif entryA < entryB # && entryA < x
            krA += 1
            if krA > endA
                # We exhausted A prior to both dropping any entries from A and injecting x.
                # Hence we write x at the end of A and return immediately.
                setstorage[krA] = x
                return lengthA+1
            end
            entryA = setstorage[krA]
        elseif entryA > entryB # && entryA < x
            krB += 1
            if krB > endB
                # We exhausted B prior to both dropping any entries from A and injecting x.
                # Hence completing this operation reduces to a strict insertion of x
                # into the remainder of A.
                firstkremA = krA
                lengthremA = endA - krA + 1
                rightshiftinsert!(setstorage, firstkremA, lengthremA, x)
                return lengthA+1
            end
            entryB = setstorage[krB]
        else # entryA == entryB && entryA < x
            # We are dropping an entry from B in A prior to injecting x. Hence we can drop
            # the entry from A, then complete the operation by calling shiftstripinsert!
            # on the remainder of A, remainder of B, and x.
            lengthpreA = krA - firstkA
            firstkremA = krA + 1
            firstkremB = krB + 1
            lengthremA = endA - krA
            lengthremB = endB - krB
            lengthremA = shiftstripinsert!(setstorage, firstkremA, lengthremA, krA, firstkremB, lengthremB, x)
            return lengthpreA + lengthremA
        end
        # FUTURE Could accelerate the second and third cases with inner whiles?
    end
    # To reach this point, we must have found an entry in A greater than x prior to dropping
    # any entries from A. We proceed as follows: Write x to krA, the location of entryA > x.
    # Advance krA; krA becomes the next-read-and-write position rather than the
    # presently-read position. Rewrite A via the entryA cache till we either (1) exhaust A,
    # (2) exhuast B, or (3) find an entry from B in A. If we (1) exhaust A, write the
    # cached entry and return. If we (2) exhaust B, finish rewriting A by way of the cache
    # entry and then return. If we (3) find an entry from B in A (specifically, we find an
    # entry from B matching the cached entryA which we are presently comparing), then drop
    # the cached entry and call strip! on the remainders of A and B to finish the op.
    setstorage[krA] = x
    krA += 1
    exhaustedA = false
    if krA > endA # Case (1), exhausted A
        exhaustedA = true # handled below
    else
        while true
            if entryA < entryB
                temp = setstorage[krA]
                setstorage[krA] = entryA
                entryA = temp
                krA += 1
                if krA > endA # Case (1), exhausted A
                    exhaustedA = true
                    break
                end
            elseif entryA > entryB
                krB += 1
                if krB > endB # Case (2), exhausted B
                    break # to finish rewriting A via the cache entry below and return
                end
                entryB = setstorage[krB]
            else # entryA == entryB, case (3)
                lengthpreA = krA - firstkA
                firstkremA = krA
                lengthremA = endA - krA + 1
                firstkremB = krB + 1
                lengthremB = endB - krB
                lengthremA = linearsetdiff!(setstorage, firstkremA, lengthremA, setstorage, firstkremB, lengthremB)
                return lengthpreA + lengthremA
            end
        end
    end
    # Code to handle exhaustion of A
    # We need write entryA, but only if entryA does not exist in B.
    exhaustedA && while true
        if entryA < entryB # entryA does not exist in B
            setstorage[krA] = entryA
            return lengthA + 1
        elseif entryA > entryB
            krB += 1
            if krB > endB # exhausted B, entryA does not exist in B
                setstorage[krA] = entryA
                return lengthA + 1
            end
            entryB = setstorage[krB]
        else # entryA == entryB
            return lengthA # entryA exists in B, so we drop this last entry
        end
    end
    # To reach this point, we found an entry in A greater than x prior to dropping any
    # entries from A, injected x, and then exhausted B without dropping any entries from
    # A, i.e. case(2) just above. Here we finish rewriting A vai the cache entry.
    while krA <= endA
        temp = setstorage[krA]
        setstorage[krA] = entryA
        entryA = temp
        krA += 1
    end
    setstorage[krA] = entryA
    return lengthA + 1
end
"""
    shiftstripinsert!(setstorage, firstkA, lengthA, destk, firstkB, lengthB, x)

Strip set `B` from set `A` while inserting `x`, writing the result at `destk` in 
`setstorage`. Assumes `x` is not in set `A`. Assumes `destk` falls to the left of firstkA,
that is `destk < firstkA`.
"""
function shiftstripinsert!(setstorage, firstkA, lengthA, destk, firstkB, lengthB, x)
    # If set A is emtpy, x and return
    lengthA == 0 && (setstorage[destk] = x; return one(eltype(setstorage)))
    # If set B is empty, weaken this to a shiftinsert! for x 
    lengthB == 0 && return shiftinsert!(setstorage, firstkA, lengthA, destk, x)
    # FUTURE shiftinsert! should get rid of return value, this call site should compute
    #
    # FUTURE Consider grossly simplifying this at expense of some additional ops/branches.
    #
    # This operation consists of three stages. In stage (1), we have not exhausted set B
    # nor have we injected entry x. After stage (1), we proceed to either stage (2) if
    # we've injected x but note exhausted set B, or stage (3) if we've exhausted set B
    # but not injected x. In either case we subsequently proceed to stage (3), where we've
    # both exhausted set B and injected x, and we merely need shift the remainder of A.
    kwA = destk
    krA = firstkA
    krB = firstkB
    endA = firstkA + lengthA - 1
    endB = firstkB + lengthB - 1
    entryA = setstorage[krA]
    entryB = setstorage[krB]
    jumpstage = 0
    #
    # Stage 1. Neither B exhausted nor x injected.
    # FUTURE Could break this stage into two parts, avoiding some branches by comparing
    # x with entryB on entryB advancement rather than entryA with both x and entryB.
    # Wait on this, see whether it matters before introducing additional complexity.
    # We do this sort of splitting in the complex linearsetdiff! routine for stripping
    # both a set and an auxiliary entry, see that code for an outline if need be.
    while true
        if entryA > x
            setstorage[kwA] = x
            kwA += 1
            jumpstage = 2
            break # -> S2
        elseif entryA < entryB # && entryA < x
            setstorage[kwA] = entryA
            kwA += 1
            krA += 1
            if krA > endA
                setstorage[kwA] = x
                return kwA - destk + 1 # done
            end
            entryA = setstorage[krA]
            continue # -> S1
        elseif entryA > entryB # && entryA < x
            krB += 1
            if krB > endB
                jumpstage = 3
                break # -> S3
            end
            entryB = setstorage[krB]
            continue # -> S1
        else # entryA == entryB && entryA < x
            krA += 1
            if krA > endA
                setstorage[kwA] = x
                return kwA - destk + 1 # done
            end
            entryA = setstorage[krA]
            krB += 1
            if krB > endB
                jumpstage = 3
                break # -> S3
            end
            entryB = setstorage[krB]
            continue # -> S1
        end
        # FUTURE Could accelerate the second and third cases with inner whiles
    end
    #
    # Stage 2. B not exhausted, x injected.
    jumpstage == 2 && while true
        if entryA < entryB
            setstorage[kwA] = entryA
            kwA += 1
            krA += 1
            krA > endA && return kwA - destk # done
            entryA = setstorage[krA]
            continue # -> S2
        elseif entryA > entryB
            krB += 1
            if krB > endB
                jumpstage = 4
                break # -> S4
            end
            entryB = setstorage[krB]
            continue # -> S2
        else # entryA == entryB
            krA += 1
            krA > endA && return kwA - destk # done
            entryA = setstorage[krA]
            krB += 1
            if krB > endB
                jumpstage = 4
                break # -> S4
            end
            entryB = setstorage[krB]
            continue # -> S2
        end
        # FUTURE Could accelerate the first two cases with inner whiles
    end
    #
    # Stage 3. B exhausted, x not injected.
    if jumpstage == 3
        while entryA < x
            setstorage[kwA] = entryA
            kwA += 1
            krA += 1
            if krA > endA
                setstorage[kwA] = x
                return kwA - destk + 1 # done
            end
            entryA = setstorage[krA]
        end
        # entryA > x
        setstorage[kwA] = x
        kwA += 1
        # -> S4
    end
    #
    # Stage 4. B exhausted and x injected.
    while true
        setstorage[kwA] = entryA
        kwA += 1
        krA += 1
        if krA > endA
            return kwA - destk # done
        end
        entryA = setstorage[krA]
    end
    # FUTURE Simplify this while
end

"""
    copyexcept!(dest, deststart, src, srcstart, copynum, exceptel)

Copy `copynum` elements from `src` beginning with `srcstart` to `dest` beginning at
`deststart`. Do not copy the first instance of element `exceptel` encountered in `src`.
"""
function copyexcept!(dest, destoffset, src, srcoffset, copynum, exceptel)
    kR = srcoffset; endR = srcoffset + copynum - 1
    kW = destoffset
    while kR <= endR
        if src[kR] == exceptel
            kR += 1
            break
        else
            dest[kW] = src[kR]
            kR += 1; kW += 1
        end
    end
    # while kR <= endR
    #     dest[kW] = src[kR]
    #     kR += 1; kW += 1
    # end
    # Instead:
    if kR <= endR
        copy!(dest, kW, src, kR, endR - kR + 1)
    end
    return # kill type instability warning
end

end