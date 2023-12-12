using BlockArrays, LinearAlgebra, Statistics

function barbell(n)
    W1 = Float64.(ones(n,n) - I)
    W = BlockArray{Float64}(undef_blocks, [n,n], [n,n])
    setblock!(W,W1,1,1)
    setblock!(W,W1,2,2)
    setblock!(W,zeros(n,n),1,2)
    setblock!(W,zeros(n,n),2,1)
    W = Array(W)
    W[1,2*n], W[2*n,1] = 1.0, 1.0;

    return  Graph(W)
end

function french_de_groot_lazy(P, x0, T; ϵ=1e-3, target=[])
    N = length(x0)

    x = x0
    t = 0
    P_evolution = (P+I)./2
    targetnorm = norm(target)

    while true
        x = P_evolution*x
        t += 1
        (norm(x-target)/targetnorm ≤ ϵ || t≥T) && break
    end
    return t
end