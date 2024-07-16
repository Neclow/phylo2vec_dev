"Get length of digit `d` if represented as a string."
digitlen(d) = ceil(Int, log10(d + 1)) + (d == 0)

"""
    sample_v(nleaves; ordered=false)

Sample a random tree via Phylo2Vec
"""
function sample_v(nleaves; ordered=false)
    if ordered
        [rand(0:i) for i in 0:nleaves-2]
    else
        [rand(0:2i) for i in 0:nleaves-2]
    end
end