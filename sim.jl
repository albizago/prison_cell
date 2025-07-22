#!/opt/homebrew/bin/julia

# DEPENDENCIES

import Plots
const plt = Plots
import Plots: @gif
import Plots: @animate
plt.gr()
plt.default(fontfamily="Computer Modern",
    titlefont=plt.font(50, "Computer Modern"),
    guidefont=plt.font(44, "Computer Modern"),
    tickfont=plt.font(42, "Computer Modern"),
    legendfont=plt.font(42, "Computer Modern"),
)

import Random as rnd
import Base: ==, hash
using Profile, Dates, DataStructures, LaTeXStrings, Measures

rnd.seed!()

# PARAMETERS

const dim::Int32 = 250 # grid size
const tsteps::Int32 = 100 # duration
const bcs::Bool = false # true = periodic, false = finite
# neighbourhood = true # true = Moore, false = Von Neumann
const simm::Bool = true # true = central squares, false = random condition
const prob::Float64 = 0.9 # portion of C in initial condition. Used if simm=false

b::Float64 = 1.6 # b parameter

const S::Float64, P::Float64, R::Float64 = 0., 0., 1.
global T::Float64 = b

# USER-DEFINED CLASSES
mutable struct Agent 
    # true = C, false = D
    kind::Bool
    prev_kind::Bool
    payoff::Float64

    # Base.show(io::IO, a::Agent) = print(io, a.payoff)
end

function ==(a::Agent, b::Agent)
    return a.kind == b.kind
end

function hash(a::Agent, h::UInt)
    h = hash(a.kind, h)
    h = hash(a.prev_kind, h)
    h = hash(a.payoff, h)
    return h
end

function value(a::Agent)
    # In order: D->D, C->D, D->C, C->C
    return a.prev_kind + 2*a.kind
end 

function get_pay(a::Agent)
    return a.payoff
end 

function rnd_ag()
    rnd_bool = rand() < prob
    Agent(rnd_bool, rnd_bool, 0.)
end

# GRID SETUP

global prisoners = Array{Agent, 2}(undef, dim, dim)

if simm
    for ia in eachindex(prisoners)
        prisoners[ia] = Agent(true, true, 0.)
    end
    half = div(dim, 2)
    for pp in prisoners[half:half+1, half:half+1]
        pp.kind = false
    end
else
    for ia in eachindex(prisoners)
        prisoners[ia] = rnd_ag()
    end
end

global coop = Array{Float64,1}()
# global cooperators = Array{Agent, 2}(undef)

global aff = Array{Float64,1}()

global ent = Matrix{Float64}(undef, 3, 0) # Block entropies

# FUNCTIONS

function neighbours(x,y)
    # central cell is excluded
    if bcs
        nbrs = [prisoners[mod1(i, dim), mod1(j,dim)] for i in y-1:y+1,
                                j in x-1:x+1
                                if (i != y || j != x)]
    else
        nbrs = [prisoners[i, j] for i in max(1, y-1):min(dim, y+1),
                                j in max(1, x-1):min(dim, x+1)
                                if (i != y || j != x)]
    end
    return nbrs
end

function payoff!(ag::Agent, nbrs)
    ag.payoff = 0.
    for nbr in nbrs
        if ag.kind
            ag.payoff += (nbr.kind ? R : S)
        else
            ag.payoff += (nbr.kind ? T : P)
        end
    end
    ag.prev_kind = ag.kind
end

function update!(ag::Agent, nbrs)
    aff_n::Float64 = 0
    change::Bool = ag.kind
    max::Float64 = -Inf
    for nbr in rnd.shuffle(nbrs)
        if nbr.payoff > max
            max = nbr.payoff
            change = nbr.prev_kind
        end
        if nbr.prev_kind == ag.kind
            aff_n += 1 / (dim*dim)
        end
    end
    ag.kind = change
    return aff_n
end

function global_update()
    push!(aff,0)
    push!(coop,0)

    # Calculate payoff
    for ix in eachindex(prisoners)
        if prisoners[ix].kind
            coop[end] += 1 / (dim*dim)
        end
        y, x = Tuple(CartesianIndices(prisoners)[ix])
        payoff!(prisoners[ix], neighbours(x,y))
    end

    # Update strategy
    for id in eachindex(prisoners)
        y, x = Tuple(CartesianIndices(prisoners)[id])
        aff[end] += update!(prisoners[id], neighbours(x,y))
    end
end

function block_entropy()
    global ent
    scales = [3,5,10]
    entropies = [0.,0.,0.]
    for index in 1:3
        scale = scales[index]
        blocks = Vector{Vector{Agent}}()
        for i in 1:(dim-scale+1)
            for j in 1:(dim-scale+1)
                block = prisoners[i:i+scale-1, j:j+scale-1]
                push!(blocks, vec(block))
            end
        end
        blk = counter(blocks)
        tot = length(blocks)
        for p in values(blk)
            entropies[index] += - (1 / scale * scale) * (p / tot) * log2(p / tot)
        end
    end
    ent = hcat(ent, entropies)
end

# DISPLAY

"""
# STATIC
plot = plt.heatmap(value.(prisoners), 
            legend=nothing,
            xticks=nothing, yticks=nothing, 
            aspect_ratio=1, size=(500,500),
            c = :thermal)
plt.display(plot)
"""

"""
# TEST 
for i in 1:5
println("\n")
println(get_pay.(neighbours(11,11)))
println(get_pay.(neighbours(9,9)))
println("\n")
plot = plt.heatmap(value.(prisoners), 
            legend=nothing,
            xticks=nothing, yticks=nothing, 
            aspect_ratio=1, size=(500,500),
            c = :thermal)
plt.display(plot)
global_update()
end
"""

anim = @animate for i ∈ 1:tsteps
    # @time global_update() # monitor speed
    global_update()
    plt.heatmap(value.(prisoners), 
            legend=nothing,
            xticks=nothing, yticks=nothing, 
            aspect_ratio=1, size=(900,900),
            c = :roma)
    if i % 100 == 0
        now_str = Dates.format(now(), "yyyy-mm-dd_HHMMSS")
        plt.savefig("sim/frame_" * now_str * ".png")
    end
    @time block_entropy()
end

# println(coop)
# println(vec(ent))

# OUTPUT

now_str = Dates.format(now(), "yyyy-mm-dd_HHMMSS")

gifname = "sim/anim_" * now_str * ".gif"

plt.gif(anim, gifname, fps = 12)

x_coord = 1:tsteps
plt.plot(x_coord, aff, size=(3000,2100), title="Vicini simili (media)", xlabel="t",margin=40mm,w=8, label=nothing, palette=:berlin)
plt.savefig("sim/aff_" * now_str * ".png")

plt.plot(x_coord, coop, size=(3000,2100), title="Frazione di C", xlabel="t",margin=40mm, w=8, label=nothing, palette=:reds)
plt.savefig("sim/coop_" * now_str * ".png")

plt.plot(x_coord, transpose(ent), size=(3000,2100), title="Block entropies",margin=40mm, w=8, xlabel="t",label=[L"$S_3$" L"$S_5$" L"$S_{10}$"], palette=:Dark2_5)
plt.savefig("sim/ent_" * now_str * ".png")

# println(get_pay.(prisoners))

@profile global_update()

# plt.display(plot)