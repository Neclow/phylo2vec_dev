include("base.jl")

"""
    get_ancestry(v)
    
Get the "ancestry" of v (see "Returns" paragraph)

v[i] = which BRANCH we do the pairing from

The initial situation looks like this:
                    R
                    |
                    | --> branch 2
                  // \\
    branch 0 <-- //   \\  --> branch 1
                0     1

For v[1], we have 3 possible branches too choose from.
v[1] = 0 or 1 indicates that we branch out from branch 0 or 1, respectively.
The new branch yields leaf 2 (like in ordered trees)

v[1] = 2 is somewhat similar: we create a new branch from R that yields leaf 2
"""
function _vector_to_ancestry(v)
    pairs = [(0, 1)]

    for i in 2:length(v)
        next_leaf = i

        if v[i] <= i - 1
            insert!(pairs, 1, (v[i], next_leaf))
        else
            insert!(pairs, v[i] - length(pairs) + 1, (pairs[v[i]-length(pairs)][1], next_leaf))
        end
    end

    ancestry = zeros(Int, length(pairs), 3)

    parents = Dict{Int,Int}()
    next_parent = length(v) + 1

    for (i, pair) in enumerate(pairs)
        child1, child2 = pair

        parent_child1 = get(parents, child1, child1)
        parent_child2 = get(parents, child2, child2)

        ancestry[i, :] = [parent_child1, parent_child2, next_parent]

        parents[child1] = next_parent
        parents[child2] = next_parent

        next_parent += 1
    end

    ancestry
end


"Build newick tree from `ancestry` with recursion."
function _build_newick(ancestry)
    nmax = size(ancestry, 1)

    function build_newick_inner(p::Int)::String
        c1, c2 = @view ancestry[p-nmax, 1:2]

        left = c1 > nmax ? build_newick_inner(c1) : string(c1)
        right = c2 > nmax ? build_newick_inner(c2) : string(c2)

        return string("(", left, ",", right, ")", p)
    end

    build_newick_inner(ancestry[end, end]) * ";"
end

"Convert vector `v` to Newick representation."
function to_newick(v)
    ancestry = _vector_to_ancestry(v)
    _build_newick(ancestry)
end







"Build newick tree from `ancestry` without recursion."
function _build_newick_old(ancestry)
    c1, c2, p = ancestry[end, :]

    newick = "($c1,$c2)$p;"

    node_idxs = Dict{Int,Int}()
    node_idxs[c1] = 1
    node_idxs[c2] = 2 + digitlen(c1)

    queue = []

    nmax = size(ancestry)[1]

    c1 > nmax && push!(queue, c1)
    c2 > nmax && push!(queue, c2)

    for _ in 2:nmax
        next_parent = pop!(queue)

        c1, c2, p = ancestry[next_parent-nmax, :]

        sub_newick = "($c1,$c2)$p"

        newick = newick[1:node_idxs[p]] * sub_newick * newick[node_idxs[p]+digitlen(p)+1:end]

        node_idxs[c1] = node_idxs[p] + 1
        node_idxs[c2] = node_idxs[c1] + 1 + digitlen(c1)

        c1 > nmax && push!(queue, c1)
        c2 > nmax && push!(queue, c2)
    end

    newick
end

function to_newick_old(v)
    ancestry = _vector_to_ancestry(v)
    _build_newick_old(ancestry)
end