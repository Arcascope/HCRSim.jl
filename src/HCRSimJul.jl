module HCRSimJul


using Distributions, OrdinaryDiffEq
using .Threads, Plots
using LinearAlgebra
using DataFrames


export CircadianModel, SinglePopModel, @dudt


abstract type CircadianModel end
abstract type LightFunc <: Function end


struct SinglePopModel <: CircadianModel
	θ::AbstractVector
	F!::Function
	γ::Function #Gives the phase from the state variables
	Amp::Function #Gives the amplitude of the state variables
	DLMOObs::Function
	CBTObs::Function
	τDD::Function
	AmpDD::Function
end



macro dudt(Light, rhs, name)

  fn = Symbol(name*"!")

  quote
  function $(esc(fn))(du,u,params,t)
    L(t)=$(esc(Light))(t)
    $(esc(rhs))(du,u,params,t,L)
  end
 end
end




#Shared Utils
include("./Utils/utils.jl")
include("./Utils/CircularStatistics.jl")

#Light Schedules
include("./Light/LightSchedules.jl")



#Models
include("./Models/SinglePopModel.jl")
include("./Models/Forger99.jl")
#include("./Models/SinglePopModelNN.jl")
include("./Models/ModelSolver.jl")



#Plots
include("./Graphics/plotRecipes.jl")















end # module
