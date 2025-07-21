#!/opt/homebrew/bin/julia

# DEPENDENCIES

import Plots
const plt = Plots
import Plots: @gif
import Plots: @animate
plt.gr()
import Random as rnd
using Profile

rnd.seed!()

# PARAMETERS

const dim::Int32 = 100 # grid size
const bcs::Bool = true # true = periodic, false = finite
# neighbourhood = true # true = Moore, false = Von Neumann
const simm::Bool = false # true = central squares, false = random condition
const prob::Float64 = 0.9 # portion of C in initial condition. Used if simm=false

const b::Float64 = 1.5 # b parameter

const S::Float64, P::Float64, R::Float64, T::Float64 = 0., 0., 1., b

# USER-DEFINED CLASSES

mutable struct Agent 
    # true = C, false = D
    kind::Bool
    prev_kind::Bool
    payoff::Float64

    # Base.show(io::IO, a::Agent) = print(io, a.payoff)
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
    change::Bool = ag.kind
    max::Float64 = -Inf
    for nbr in rnd.shuffle(nbrs)
        if nbr.payoff > max
            max = nbr.payoff
            change = nbr.prev_kind
        end
    end
    ag.kind = change
end

function global_update()

    push!(coop, 0)

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
        update!(prisoners[id], neighbours(x,y))
    end

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

anim = @animate for i ∈ 1:500
    @time global_update() # monitor speed
    plt.heatmap(value.(prisoners), 
            legend=nothing,
            xticks=nothing, yticks=nothing, 
            aspect_ratio=1, size=(500,500),
            c = :roma)
end

println(coop)

plt.gif(anim, "anim_fps15.gif", fps = 10)

# println(get_pay.(prisoners))

@profile global_update()

# plt.display(plot)