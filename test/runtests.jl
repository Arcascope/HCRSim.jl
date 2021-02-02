
using HCRSimJul
using Test

isaNumVec(v)=all(isa.(v, Number))

LightFuncs=[HCRSimJul.RegularLight, HCRSimJul.SlamShift, HCRSimJul.ShiftWorkLight, HCRSimJul.SocialJetLag]
models=[SinglePopModel(), Forger99()]



@testset "Integration_Transients" begin

	for f in LightFuncs
		for m in models
			@test isaNumVec(integrate_transients(f, m)) #test is all returned values are numbers
		end
    end
end; 



@testset "Integration Models" begin 

	for f in LightFuncs
		for m in models
			@test isaNumVec(integrate_model(f,m, [0.70,1.0,0.0], 0.0, 500.0)[end]) #test is all returned values are numbers
		end
	end
	
end; 



@testset "Integration Observer DLMO" begin 

	for f in LightFuncs
		for m in models
			@test isaNumVec(integrate_observer(f,m, [0.70, 0.0,0.0], 0.0, 500.0, m.DLMOObs)) #test is all returned values are numbers
		end
	end
end; 



@testset "Integration Observer CBT" begin 

	for f in LightFuncs
		for m in models
			@test isaNumVec(integrate_observer(f,m, [0.70, 0.0,0.0], 0.0, 500.0, m.CBTObs)) #test is all returned values are numbers
		end
    end
end; 


@testset "Parameter Updates SinglePopModel" begin 

	p1=SinglePopModel(Dict("τ"=> 25.0, "γ"=> 10.0, "σ"=> -5.0))
	@test p1.θ[1] ≈	 25.0 
	@test p1.θ[3] ≈ 10.0 
	@test p1.θ[9] ≈ -5.0

end;




@testset "Entrainment Test" begin 

	for m in models
		dlmo=integrate_observer(HCRSimJul.RegularLight, m, [0.70, 1.0,0.0], 0.0, 1000.0, m.DLMOObs) #test is all returned values are numbers
		@test isapprox(diff(dlmo)[end], 24.0; atol=1e-2)
	end
end;


@testset "Model Comparison Tests CBT/DLMO" begin 

	cbt1=integrate_observer(HCRSimJul.RegularLight, models[1], [0.70, 1.0,0.0], 0.0, 1000.0, models[1].CBTObs)
	cbt2=integrate_observer(HCRSimJul.RegularLight, models[2], [0.70, 1.0,0.0], 0.0, 1000.0, models[2].CBTObs)

	@test isapprox((cbt1[end] - cbt2[end])%24.0, 0.0; atol=1.0)

	dlmo1=integrate_observer(HCRSimJul.RegularLight, models[1], [0.70, 1.0,0.0], 0.0, 1000.0, models[1].DLMOObs)
	dlmo2=integrate_observer(HCRSimJul.RegularLight, models[2], [0.70, 1.0,0.0], 0.0, 1000.0, models[2].DLMOObs)

	@test isapprox((dlmo1[end] - dlmo2[end])%24.0, 0.0; atol=1.0)

end 


@testset "ic finder" begin 

	for m in models 
		@test isaNumVec(HCRSimJul.ic_from_clocktime(m, 8.0))
	end

	
end 




