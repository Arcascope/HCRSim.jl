

#=
Implements the Hilaire VDP Model from 2007



=#



struct Hilaire2007Model <: CircadianModel
	θ::AbstractVector
	F!::Function
	Amp::Function
	γ::Function #Gives the phase from the state variables
	DLMOObs::Function
	CBTObs::Function
	τDD::Function
	AmpDD::Function
end

function Hilaire2007Model(θ::AbstractVector)
	γ(state)=state[2]
	CBTObs(u,t,integrator)=angle(exp(im*u[2]+im*π))
	DLMOObs(u,t,integrator)=angle(exp(im*u[2]))-5.0/12.0*π
	freerunningτ(x)=24.2
	freerunningAmp(x)=1.0
	Hilaire2007Model(θ, spmodelNN!, γ, DLMOObs, CBTObs, freerunningτ, freerunningAmp)
end


#Default constructor for SinglePopModel
function Hilaire2007Model()
	#0.26399938265460443
	θ=[23.8, 0.0635842, 0.024,-0.0931817, 0.3855, 0.195123,-0.0026,-0.957756,0.0400692,0.05,1.5,9325.0,0.0075, 33.75]
	θnn=zeros(200)
	Hilaire2007Model(vcat(θ, θnn))
end

function Hilaire2007Model(γ::AbstractFloat)
	τ=rand(Normal(23.8, γ))
	θ=[τ, 0.0635842, 0.024,-0.0931817, 0.3855, 0.195123,-0.0026,-0.957756,0.0400692,0.05,1.5,9325.0,0.0075, 33.75]
	θnn=zeros(200)
	Hilaire2007Model(vcat(θ, θnn))
end





function hilaire2007!(du,u,θ,t,L,SW)

	τc,G,k,μ,β,q,ρ,I0,p,a0=θ

	x=u[1]
	xc=u[2]
	n=u[3]

	αNonPhotic(I)=a0*(I/I0)^p*(I/(I+100)

    C = t % 24;
    phi_xcx = -2.98;
    phi_ref = 0.97;
    CBTmin = phi_xcx + phi_ref;
    CBTmin = CBTmin*24/(2*π);
    psi_cx = C - CBTmin;
    psi_cx = psi_cx % 24;

    Bh = G * (1 - n) * αNonPhotic(I)
    B = Bh * (1 - .4 * x) * (1 - .4 * xc);

    #Subtract from 1 to make the sign work
    #From St. Hilaire (2007): sigma equals either 1 (for sleep/rest) or 0 (for wake/activity),
    sigma = 1 - SW(t);
    if (sigma < 1/2)
        sigma = 0;
    else
        sigma = 1;
    end

    Nsh = ρ*(1/3 - sigma);

    if (psi_cx > 16.5 && psi_cx < 21)
        Nsh = rho*(1/3);
    else
        Ns = Nsh*(1 - tanh(10*x));
    end

    du[1] = π / 12.0 * (xc + μ*((1/3)*x + (4/3)*x^3 - 256/105*x^7.0) + B + Ns);
    du[2] = π / 12.0 * (q*B*xc - x*((24/(0.99729*τc))^2 + k*B));
    du[3] = 60.0 * (αNonPhottic(abs(L(t))) * (1.0 - n) -β * n);

end



function integrateTransients(myLight,model::Hilaire2007Model; initial=[1.0,0.0,0.0], tstart=0.0, tend=1000.0)

	sol=integrateModel(myLight,model,initial, tstart,tend)
	last=copy(sol.u[end])
	last[2]=last[2] % (2π)
    return last
end
