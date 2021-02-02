#=
Single Population Model Julia Code for Simulation and Analysis

Some Useful translations for the single population model to biological timepoints
CBT=DLMO+7hrs
CBT=DLMO_mid+2hrs
CBT=circadian phase pi in the model
DLMO=circadian phase 5pi/12=1.309 in the model


=#




export spmodel!


#----------------Constructors-------------------------------------


function SinglePopModel(θ::AbstractVector)
	γ(state)=state[2]
	Amp(state)=state[1]
	CBTObs(u,t,integrator)=angle(exp(im*u[2]+im*π))
	DLMOObs(u,t,integrator)=angle(exp(im*u[2]))-5.0/12.0*π
	freerunningτ(x)=2π/(x[1]+x[2]*sin(x[4])*(1-x[3]/(x[2]*cos(x[4]))))
	freerunningAmp(x)=(1-2*x[3]/(x[2]*cos(x[4])))^(0.25)
	SinglePopModel(θ, spmodel!, γ, Amp, DLMOObs, CBTObs, freerunningτ, freerunningAmp)
end


#Default constructor for SinglePopModel
function SinglePopModel()
	#0.26399938265460443
	θ=[23.8, 0.0635842, 0.024,-0.0931817, 0.3855, 0.195123,-0.0026,-0.957756,0.0400692,0.05,1.5,9325.0,0.0075, 33.75]
	#θnew=[23.808664915808794, 0.06335583729236834, 0.026530273202541896, -0.09895651890451312, 0.39026016291832744, 0.20684411280731035, -0.005195405934242825, -0.9412929484477952, 0.017166548345910106, 0.04977857816684212, 1.49854218815377, 9325.000915697265, 0.00687102399745335, 33.75]
	SinglePopModel(θ)
end

function SinglePopModel(γ::AbstractFloat)
	ω0=rand(Normal(23.8, γ))
	θ=[ω0, 0.0635842, 0.024,-0.0931817, 0.3855, 0.195123,-0.0026,-0.957756,0.0400692,0.05,1.5,9325.0,0.0075, 33.75]
	SinglePopModel(θ)
end

function SinglePopModel(params::Dict)
	SinglePopModel(build_params_sp(params))
end



function spmodel!(du,u,params,t,L)

	τ,K,γ, β1, A1, A2, βL1, βL2, σ, αconst, p, I0, δ, G=params

	R=u[1]
	Ψ=u[2]
	n=u[3]

	α0(t)=αconst*abs(L(t))^p/(abs(L(t))^p+I0)
	B(t)=33.75*(1.0-n)*α0(t)

	du[1]=-1.0*γ*R+K*cos(β1)/2.0*R*(1.0-R^4)+B(t)*(A1*0.5*(1.0-R^4)*cos(Ψ+βL1)+A2*0.5*R*(1.0-R^8)*cos(2.0*Ψ+βL2))
	du[2]=2π/τ+K/2.0*sin(β1)*(1+R^4)+B(t)*(σ-A1*0.5*(R^3+1.0/R)*sin(Ψ+βL1)-A2*0.5*(1.0+R^8)*sin(2.0*Ψ+βL2))
	du[3]=60.0*(α0(t)*(1.0-n)-δ*n)

end


"""
integrate_transients(myLight,model::SinglePopModel; initial=[1.0,0.0,0.0], tstart=0.0, tend=24.0*100)

Integrates the model for a long time to get rid of any transients
"""
function integrate_transients(myLight,model::SinglePopModel; initial=[1.0,0.0,0.0], tstart=0.0, tend=24.0*100)

	sol=integrate_model(myLight,model,initial, tstart,tend)
	last=copy(sol.u[end])
	last[2]=last[2] % (2π)
    return last
end




"""
This utility function returns a parameter vector for the single pop model given
a dictionary of parameter values. For any values not specified in the dict
will be filled in from the the optimal parameter set as defined in the default constructor.

function build_params_sp(paramsDict::Dict)::AbstractFloat

"""
function build_params_sp(paramsDict::Dict)::AbstractVector

    m=SinglePopModel()
    p1=deepcopy(m.θ)

    pnames=["τ","K","γ", "β1", "A1", "A2", "βL1", "βL2", "σ", "αconst", "p", "I0", "δ", "G"]
    ParamIdxDict=Dict(zip(pnames, 1:length(pnames)))

    for p in keys(paramsDict)
        p1[ParamIdxDict[p]]=paramsDict[p]
    end

    return(p1)

end
