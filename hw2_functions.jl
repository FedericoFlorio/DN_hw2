using Revise, LinearAlgebra, BlockArrays, Statistics
includet("hw1_functions.jl")

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

exponential(r) = -log(rand())/r     # sampler for exponential distribution of rate r (inverse cumulative method)
∂(A,i) = A.rowval[nzrange(A,i)]     # returns the out neighbours of node i

function french_de_groot(P, x0, T; stats = (t,x)->println("$t $x"), ϵ=1e-6)
    stop(x,X,eps) = (norm(x-X)/norm(X) ≤ eps)
    N = length(x0)

    function simulation(P, x0, T, stats, ϵ)
        x = x0
        t = 0

        while true
            x0 = x
            x = P*x0
            t += 1
            stats(t,x)
            (stop(x,x0,ϵ) || t≥T) && break
        end
        return t
    end
    return simulation(P, x0, T, stats, ϵ)
end

################################################################################################

nstart = 6
nstop = 10
niter = 100
T = 1e3

nn = [2^i for i in nstart:nstop]
tt = Vector{Int}(undef,length(nn))

function do_nothing(t,x)
    return nothing
end

for n in nn
    G = barbell(n)
    times_abs = Vector{Int}(undef,niter)
    x_in = [(mod(i,2)==0 ? 0 : 1) for i in 1:2*n]

    for iter in 1:niter
        times_abs[iter] = french_de_groot(G.P, x_in, T, stats=do_nothing)
    end
    tt[searchsortedfirst(nn,n)] = mean(times_abs)
end

println(nn)
println(tt)