include("base.jl")

function _newick_to_ancestry(newick)
    # pre-allocated ancestry matrix based on number of nodes
    # note: matrix is transposed since Julia operates faster on columns than rows
    ncols = count(==(')'), newick)
    ancestry = fill("", 3, ncols)

    # pot. further speed-ups (to be tested):
    # - replace split with find index of comma, then take substring
    # - replace findfirst, findprev, findnext with loop 

    for col in eachcol(ancestry)
        # get first closing parenthesis...
        close_idx = findfirst(')', newick)
        # ...then go backwards to find matching opening parenthesis
        # removing 4 digits to account for closing brace + 2 numbers (min 1 digit) + comma
        open_idx = findprev('(', newick, close_idx - 4)

        # get c1 and c2 by finding location of comma
        comma_idx = findnext(',', newick, open_idx + 1)
        c1 = @view newick[open_idx+1:comma_idx-1]
        c2 = @view newick[comma_idx+1:close_idx-1]

        # get parent as digits after closing parenthesis but before non-digit character
        p = @view newick[close_idx+1:findnext(!isdigit, newick, close_idx + 1)-1]

        # assign to matrix
        col[1] = c1
        col[2] = c2
        col[3] = p

        # remove processed substring from big string
        @views newick = newick[1:open_idx-1] * newick[close_idx+1:end]
    end

    # transpose matrix back and convert strings to int
    return parse.(Int, ancestry)'
end

function _order_cherries(ancestry)
    # pre-allocate output
    ancestry_sorted = similar(ancestry)
    smallest_child = Dict{Int,Int}()

    # iterate on ancestry from highest parent node to lowest
    for (i, j) in enumerate(sortperm(ancestry[:, end]))
        c1, c2, p = ancestry[j, :]

        # replace value of c1, c2 by value of smallest child connected to node
        c1_smallest_child = get(smallest_child, c1, c1)
        c2_smallest_child = get(smallest_child, c2, c2)

        smallest_child[p] = min(c1_smallest_child, c2_smallest_child)

        # assign to matrix
        ancestry_sorted[i, 1] = c1_smallest_child
        ancestry_sorted[i, 2] = c2_smallest_child
        ancestry_sorted[i, 3] = max(c1_smallest_child, c2_smallest_child)
    end

    return ancestry_sorted
end

function _build_vector(cherries)
    # pre-allocate output
    n = size(cherries, 1)
    v = zeros(Int, n)

    # iterate backwards (not actually needed)
    for i = n:-1:1
        # get cherry
        c1, c2, cmax = cherries[i, :]
        idx = 0

        # find position of cherry where parent is <= current max node 
        # and where current max node is a child
        # big performance boost over using Julia's built-in findfirst
        for j in 1:n
            if cherries[j, end] <= cmax
                idx += 1
                if cherries[j, 1] == cmax || cherries[j, 2] == cmax
                    break
                end
            end
        end

        # assign correct Phylo2Vec number
        if idx == 1
            v[cmax] = min(c1, c2)
        else
            v[cmax] = cmax + idx - 2
        end
    end

    return v
end

function to_vector(newick)
    ancestry = _newick_to_ancestry(newick)
    cherries = _order_cherries(ancestry)
    vector = _build_vector(cherries)

    return vector
end









function _reduce_old(newick)
    newickstr = newick
    ancestry = zeros(Int, 0, 3)

    function do_reduce(newickstr)
        open_idx = 1

        for (i, char) in enumerate(newickstr)
            if char == '('
                open_idx = i + 1
            elseif char == ')'
                # todo max 2 splits (3 elems)
                child1, child2 = split(newickstr[open_idx:i-1], ",")[1:2]
                parent = split(split(newickstr[i+1:end], ",")[1], ")")[1]

                ancestry = vcat(
                    ancestry,
                    [parse(Int, child1) parse(Int, child2) parse(Int, parent)],
                )

                return do_reduce(newickstr[1:open_idx-2] * newickstr[i+1:end])
            end
        end
    end

    do_reduce(newick[1:end-1])

    return ancestry
end

function _order_cherries_old(ancestry)
    ancestry_sorted = ancestry[sortperm(ancestry[:, end]), :]

    small_children = Dict{Int,Int}()

    for i in 1:size(ancestry_sorted, 1)
        c1, c2, p = ancestry_sorted[i, :]

        parent_c1 = get(small_children, c1, c1)
        parent_c2 = get(small_children, c2, c2)

        small_children[p] = min(parent_c1, parent_c2)

        ancestry_sorted[i, :] = [parent_c1, parent_c2, max(parent_c1, parent_c2)]
    end

    return ancestry_sorted
end

function _build_vector_old(cherries)
    v = zeros(Int, size(cherries, 1))

    for i = size(cherries, 1):-1:1
        c1, c2 = @view cherries[i, 1:2]
        c_max = max(c1, c2)

        subset = cherries[cherries[:, end].<=c_max, 1:2]

        idx = findall(subset .== c_max)[1][1]

        if idx == 1
            v[c_max] = min(c1, c2)
        else
            v[c_max] = c_max - 1 + idx - 1
        end
    end

    return v
end