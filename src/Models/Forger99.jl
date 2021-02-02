

#=
Implements the Forger 99 Model 

The mininum of x is the CBTmin 
The DLMO time is a back from this.....



=#



export Forger99, forger99rhs!


struct Forger99 <: CircadianModel
	θ::AbstractVector
	F!::Function
	γ::Function #Gives the phase from the state variables
	Amp::Function #Gives the amplitude of the state variables
	DLMOObs::Function
	CBTObs::Function
	τDD::Function
	AmpDD::Function
end


function Forger99(θ::AbstractVector)
	γ(state)=atan(-1*state[2],state[1])
	Amp(state)=sqrt(state[1]^2+state[2]^2)
	CBTObs(u,t,integrator)=u[2] #mininum of x var approx by xc=0.0
	DLMOObs(u,t,integrator)=γ(u)-5π/12.0 
	freerunningτ(x)=24.2
	freerunningAmp(x)=1.0
	Forger99(θ, forger99rhs!, γ, Amp, DLMOObs, CBTObs, freerunningτ, freerunningAmp)
end


#Default constructor for SinglePopModel
function Forger99()
	θ=[24.2,0.23,33.75,0.16,0.0075,0.60,9500.0, 0.55]
	Forger99(θ)
end

""" 
The rhs for the forger 99 model. 

	forger99rhs!(du,u,θ,t,L)

"""
function forger99rhs!(du,u,θ,t,L)

	τx, μ, G, α0, δ, p, I0, k =θ
	x=u[1]
	xc=u[2]
	n=u[3]

	αfunc(t)=α0*(L(t)/I0)^p
	Bhat=G*(1.0-n)*αfunc(t)*(1-0.4*x)*(1-0.4*xc)
	
	du[1]=π/12.0*(xc+Bhat)
    du[2]=π/12.0*(μ*(xc-4.0/3.0*(xc^3))-x*((24.0/(0.99669*τx))^2.0+k*Bhat))
    du[3]=60.0*(αfunc(t)*(1.0-n)- δ*n)
end


"""
Integrate the transients for the forger 99 model returns the ending state

integrate_transients(myLight, model::Forger99; initial=[0.5,0.5,0.0], tstart=0.0, tend=1000.0)

returns: Final State
"""
function integrate_transients(myLight, model::Forger99; initial=[0.5,0.5,0.0], tstart=0.0, tend=1000.0)

	sol=integrate_model(myLight,model,initial, tstart,tend)
	last=copy(sol.u[end])
    return last
end




