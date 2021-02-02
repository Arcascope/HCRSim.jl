

export integrate_model, integrate_transients, integrate_observer



"""
Integrate a circadian model...

integrate_model(myLight::Function, model::CircadianModel, initial::AbstractVector, tstart::AbstractFloat, tend::AbstractFloat)

will return a ode integration object.

"""
function integrate_model(myLight::Function, model::CircadianModel, initial::AbstractVector, tstart::AbstractFloat, tend::AbstractFloat; dt::AbstractFloat=0.1)

	@dudt(myLight, model.F!, "circmodel")
    prob=ODEProblem(circmodel!,initial,(tstart,tend), model.θ)
    sol=solve(prob,RK4(),saveat=dt)
    return sol

end


"""
This is a general integrate transients method for an abstract circadian model 
	integrate_transients(myLight,model::CircadianModel, initial; tstart=0.0, tend=1000.0)
returns the final state
"""
function integrate_transients(myLight::Function, model::CircadianModel, initial::AbstractVector; tstart::AbstractFloat=0.0, tend::AbstractFloat=100*24.0)::AbstractVector

	sol=integrate_model(myLight,model,initial, tstart,tend)
	last=copy(sol.u[end])
	last[2]=last[2] % (2π)
    return last
end





function sensitivity_integration(myLight::Function,model::CircadianModel, initial::AbstractVector, tstart::AbstractFloat, tend::AbstractFloat)
	#=
	    Currently broken because the light function is included in the parameter list,
		need to be able to calculate the gradient of the state variables against the static
		parameters.

		Perhaps this is fixed using the macro method on 6/11/2019
	=#

	@dudt(myLight, model.F!, "lcmodelSens")
    prob=ODEForwardSensitivityProblem(lcmodelSens!,initial,(tstart,tend), model.θ)
    sol_sens=solve(prob, RK4())
	return sol_sens

end

"""
function integrate_observer 

	integrate_observer(myLight, model::CircadianModel, initial::AbstractVector, tstart::AbstractFloat, tend::AbstractFloat, Obs::Function)

The observer function should hit zero at the times of interest, will return a list of times where the observer function 
is zero. The return value is Float32[]


"""
function integrate_observer(myLight, model::CircadianModel, initial::AbstractVector, 
	tstart::AbstractFloat, tend::AbstractFloat, Obs::Function)::AbstractVector
	#=
	Find the CBT crossings. CBT min occurs at ψ=π in the model
	=#

	ObsTimes=[]

	function affect!(integrator)
		push!(ObsTimes, integrator.t)
	end
	

	@dudt(myLight, model.F!, "lcmodel")
	cb=ContinuousCallback(Obs, affect!, nothing, rootfind=true, interp_points=500)
	prob=ODEProblem(lcmodel!,initial,(tstart,tend), model.θ)
    sol=solve(prob, RK4(), callback=cb, dtmax=1.0)
    return ObsTimes

end



function integrate_observer(wdata::DataFrame, model::CircadianModel, Obs::Function; λ::Real=1.0, 
	use_light::Bool=true, waketime::Real=8.0, kwargs...)::AbstractVector
	#=
	Find the CBT crossings. CBT min occurs at ψ=π in the model
	=#

	if (:Lux in names(wdata) && use_light==true)
        Light=interp_flux(wdata.TimeTotal, wdata.Lux)
    else 
        Light=interp_flux(wdata.TimeTotal, λ .*	wdata.Steps)
	end 
	
	initial=ic_from_clocktime(model, wdata.TimeTotal[1]; waketime=waketime)


	return integrate_observer(Light, model, initial, wdata.TimeTotal[1], wdata.TimeTotal[end], Obs; kwargs...)

end




function integrate_strobo(myLight, model::CircadianModel, initial::AbstractVector, tstart::AbstractFloat, tend::AbstractFloat)

	R=[]
	ψ=[]

	function stroboObs(u,t,integrator)
		d=(t % 24.0)-12.0 #crosses zero once a day
	end

	function affect!(integrator)
		push!(R, model.Amp(integrator.u))
		push!(ψ, model.γ(integrator.u[2]))
	end

	@dudt(myLight, model.F!, "lcmodel")
	cb=ContinuousCallback(stroboObs,affect!, nothing, rootfind=true, interp_points=500)
	prob=ODEProblem(lcmodel!,initial,(tstart,tend), model.θ)
    sol=solve(prob,RK4(),callback=cb, dtmax=2.0)
    return [R, ψ]

end



"""
	integrateAmplitude(myLight, initial, tend; pvalues=model_param, tstart=0.0)

Record the amplitude of the model at CBTmin.

Needs to be redone for the actual model
"""
function integrate_amplitude(myLight::Function, model::CircadianModel, initial::AbstractVector, tstart::AbstractFloat, tend::AbstractFloat, Obs::Function)

	R=[]

	function affect!(integrator)
		push!(R, model.Amp(integrator.u[1]))
	end

	@dudt(myLight, model.F!, "lcmodelAmp")
	cb=ContinuousCallback(Obs,affect!, nothing, rootfind=true, interp_points=500)
	prob=ODEProblem(lcmodelAmp!,initial,(tstart,tend), model.θ)
    sol=solve(prob,RK4(), callback=cb, dtmax=2.0)
    return R

end


""" 
This function can be used to get a quick approximation for the circadian phase 
at a certain clock phase (and the initial conditions) to give a model. 

ic_from_clocktime(model::CircadianModel, clock_time::Real; waketime::Real)

returns a set of initial conditions for that clocktime. Waketime is the expected time
they first get light in the morning. Assumes they have 8 hours of dark and 16 hours 
of 150 lux light. 
"""
function ic_from_clocktime(model::CircadianModel, clock_time::Real; waketime::Real=8.0)


	Light(t)=RegularLightSimple(t; Intensity=150.0, wakeUp=waketime, workday=16.0)
	sol=integrate_model(Light, model, [0.5,0.5,0.0], 0.0, 100*24.0)
	sol_map=integrate_model(Light, model, sol.u[end], 0.0, 24.0)
	return sol_map(clock_time)

end 





