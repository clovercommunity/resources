## Preamble

# only need to do this once
using Pkg
Pkg.add("Plots")
Pkg.add("DifferentialEquations")

# use the "DifferentialEquations" package
using DifferentialEquations
# use the "Plots" package
using Plots

## Core functions

# function which defines the ordinary differential equations
function seir_ode(dY, Y, p, t)
    # unpack the parameters into their own variables
    β, σ, γ, μ = p[1], p[2], p[3], p[4]
    # unpack the state variables into their own variables
    S, E, I, R = Y[1], Y[2], Y[3], Y[4]
    # define the equations-of-state
    dY[1] = μ*(1-S)-β*S*I
    dY[2] = β*S*I-(σ+μ)*E
    dY[3] = σ*E - (γ+μ)*I
    dY[4] = γ*I - μ*R
    return dY
end

# function which integrates the differential equations
function seir(β, σ, γ, μ, tspan=(0.0,100.0); S₀=0.99, E₀=0, I₀=0.01, R₀=0)
    # create a vector containing the parameters
    par         = [β, σ, γ, μ]
    # create a vector of initial conditions
    init        = [S₀, E₀, I₀, R₀]
    # create an "ODEProblem" object using the dynamics, initial conditions,
    # simulation start & end time, and the parameters
    seir_prob   = ODEProblem(seir_ode,init,tspan,par)
    # numerically integrate the solution to the differential equations
    sol         = solve(seir_prob);
    return sol
end

r0(β, σ, γ, μ) = (σ/(σ+μ))*(β/(γ+μ))

## Perform simulations

# define parameters
μ = 0.000028
σ = 0.25
γ = 0.15
β = 100

sol = seir(β, σ, γ, μ)

plot(sol,
    xlabel  = "time (days)",
    ylabel  = "normalized population",
    label   = ["susceptible" "exposed" "infected" "recovered/deceased"],
    title   = "flatten the curve!")
