

export plot_actogram, plot_actogram!, plot_phase_marker!, dynamics_plot

#=
Plots.pyplot()

default(titlefont = (8, :black),
legendfontsize = 6,
guidefont = (6, :black),
tickfont = (6, :black))
=#

function rectangle(w, x, y)
    Shape(x .+ [0, w, w, 0], y .+ [0, 0, 1, 1])
end


function getRect(timeon, timeoff, num_days)
    bottom_x = mod(timeon, 24.0)
    bottom_y = floor(timeon / 24.0) #changed this from num_days- to just numdays
    r1 = rectangle(timeoff - timeon, bottom_x, bottom_y)
    r2 = rectangle(timeoff - timeon, bottom_x + 24.0, bottom_y)
    return (r1, r2)
end 





""" 
This function will create an actogram plot from the given "light" function 
    plot_actogram!(pl, LightFun, tstart, tend; threshold = 10.0) 
The threshold sets the difference between the light (yellow) and dark (black) 

"""
function plot_actogram!(pl::Plots.Plot, LightFun::Function, tstart::AbstractFloat, tend::AbstractFloat; threshold::AbstractFloat = 10.0, kwargs...)

    num_days = ceil((tend - tstart) / 24.0) - 1
    ynum_days=floor(tend/24.0)+1

    #Add the light schedule to the plot
    p = plot!(pl, 0:48, 0:48, label = "", yflip = true,
        xtickfontsize = 6, ytickfontsize = 6,
        xguidefontsize = 6, yguidefontsize = 6,
        titlefontsize = 8, alpha=0.0,
        kwargs...
    ) #set up a plot with the right axes
    yflip!(p, true)
    ylims!(p, 0, ynum_days)
    xlims!(p, 0, 48)
    yticks!(p, 1:ynum_days-1)
    xlabsdp = string.(vcat(collect(0:3:24), collect(3:3:24)))
    xticks!(0:3:48, xlabsdp)
    xlabel!("CT")
    ylabel!("Day")
    #title!("Actogram")

    timedata=tstart:0.25:tend
    lightdata=zeros(length(timedata))

    tsample=tstart:0.05:tend
    lightsample = LightFun.(tsample)

    for (k,tval) in enumerate(timedata)
        #take a sample
        lvals=LightFun.(tval-0.125:0.01:tval+0.125)
        lightdata[k]=maximum(lvals)
    end


    lightsOn = false

    if (lightdata[1] > threshold)
        lightsOn = true
        lightStart = timedata[1]
    else
        darkOn = true
        darkStart = timedata[1]
    end

    dayCounter = floor(timedata[1] / 24.0)

    for i = 1:length(lightdata)
        currentDay = floor(timedata[i] / 24.0)

        if (currentDay ≠ dayCounter)
            dayCounter = currentDay
            if (lightsOn)
                r1, r2 = getRect(lightStart, timedata[i], num_days)
                plot!(p, r1, color = :yellow, label = "")
                plot!(p, r2, color = :yellow, label = "")
                if (i + 1 < length(timedata))
                    lightStart = timedata[i+1]
                end
            else
                r1, r2 = getRect(darkStart, timedata[i], num_days)
                plot!(p, r1, color = :black, label = "")
                plot!(p, r2, color = :black, label = "")
                if (i + 1 < length(timedata))
                    darkStart = timedata[i+1]
                end
            end
        end

        if (lightdata[i] < threshold && lightsOn)
            r1, r2 = getRect(lightStart, timedata[i-1], num_days)
            plot!(p, r1, color = :yellow, label = "")
            plot!(p, r2, color = :yellow, label = "")
            lightsOn = false
            darkOn = true
            darkStart = timedata[i]
        end

        if (!(lightsOn) && lightdata[i] >= threshold)
            lightsOn = true
            lightStart = timedata[i]
            darkOn = false
            r1, r2 = getRect(darkStart, timedata[i-1], num_days)
            plot!(p, r1, color = :black, label = "")
            plot!(p, r2, color = :black, label = "")
        end
    end

    return p
end


function plot_actogram(LightFun::Function, tstart::AbstractFloat, tend::AbstractFloat; threshold::AbstractFloat = 10.0, kwargs...)

    pl=plot() 
    return plot_actogram!(pl, LightFun, tstart, tend; threshold = threshold, kwargs... )
end 





function plot_actogram!(pl::Plots.Plot, wdata::DataFrame; use_light::Bool=true, threshold::AbstractFloat=10.0, kwargs...)

    if (:Lux in names(wdata) && use_light==true)
        Light=interp_flux(wdata.TimeTotal, wdata.Lux)
    else 
        Light=interp_flux(wdata.TimeTotal, wdata.Steps)
    end 

    return plot_actogram!(pl, Light, wdata.TimeTotal[1], wdata.TimeTotal[end]; threshold=threshold, kwargs...)

end 


function plot_actogram(wdata::DataFrame; use_light::Bool=true, threshold::AbstractFloat=10.0, kwargs...)

    pl=plot()
    return plot_actogram!(pl, wdata; use_light=use_light, threshold=threshold, kwargs...)

end 


"""

plot_phase_marker!(pl,dlmo_times; lineplot=false, kwargs...)

Add a phase marker plot (scatter) to the given axis. This typically would 
be something like DLMO, CBTmin etc. 

Some common keyword arguments would be markersize (default=2) and lineplot::Bool whether to 
add a line to the plot
"""
function plot_phase_marker!(pl::Plots.Plot,dlmo_times::AbstractVector; markersize=2, lineplot::Bool=false, color=:blue, kwargs...)

    xvals = deepcopy(dlmo_times)
    yvals = deepcopy(dlmo_times)

    xvals = xvals .% (24.0)
    yvals = floor.(yvals / 24.0) .+ 0.5

    scatter!(
        pl,
        xvals,
        yvals,
        label = "",
        color = color,
        markersize = markersize;
        kwargs...
    )

    scatter!(
        pl,
        xvals .+ 24.0,
        yvals,
        label = "",
        color=color,
        markersize = markersize;
        kwargs...
    )
    if (lineplot)
        plot!(pl, xvals, yvals, label = "", color= color)
        plot!(pl, xvals .+ 24.0, yvals, label = "", color=color)
    end
    return pl
end



"""
Make a plot of the dynamics in phase space of the singlepop
also has side plots of the time dynamics of the amplitude and
phase.

"""
function dynamics_plot(myLight::Function, model::CircadianModel, initial::AbstractVector, tstart::AbstractFloat, tend::AbstractFloat;savefig=false)

    lay = @layout [ a{0.6w} [b; c]]

	sol=integrate_model(myLight, model, initial, tstart, tend)

	lightColor=[]
	ls=[]
	for tt in sol.t
		if myLight(tt)>10.0
			push!(lightColor, :gold)
			push!(ls, :solid)
		else
			push!(lightColor, :black)
			push!(ls, :dash)
		end
	end


    pl1=plot(sol[2,:], sol[1,:],lw=1,legend=false, linestyle=:solid, color=lightColor, proj = :polar, showaxis=false, grid=:none)
	scatter!(pl1, [0.0, 0.0], color=:darkgreen, proj=:polar, marker=:cross)

	plot!(pl1, 5π/12*ones(11), 0.0:0.1:1.0, proj=:polar, color=:blue, lw=2.0)
	annotate!(pl1, 5π/12, 1.02, text("DLMO",8,:center, :blue))

	plot!(pl1, π*ones(11), 0.0:0.1:1.0, proj=:polar, color=:red, lw=2.0)
	annotate!(pl1, 98π/100, 0.95, text("CBTmin",8,:center, :red))

	ylims!(pl1, (0,1.1))

    #radius plot
    pl2=plot(sol.t,sol[1,:], lw=1, linestyle=:solid, legend=false, color=lightColor, grid=:none)
	title!("Circadian Amplitude")
	ylabel!("Amplitude")
	xlabel!("Hrs")
	ylims!(pl2, (0,1))
	xlims!(pl2, (-2, tend))
	xticks!(pl2, 0:24:tend)



    #angle plot
	ψ=sin.(sol[2,:])
    pl3=plot(sol.t,ψ, lw=1, linestyle=:solid, legend=false, color=lightColor, grid=:none)
	title!("Circadian Phase")
	ylabel!("Sine Phase")
	xlabel!("Hrs")
	xticks!(pl3, 0:24:tend)
	ylims!(pl3, (-1, 1))
	#yticks!(pl3, 0.0:π/2:2π, ["0", "π/2", "π", "3π/2", "2π"])


    pl=plot(pl1, pl2, pl3, layout=lay)
    if (savefig)
        savefig("./Images/spdynamicsPlot.eps")
    end
    display(pl)
    return pl

end





function plot_circle!(pl::Plots.Plot, ϕ::AbstractVector)

    Z=sum(exp.(im .* p2))/length(p2)
	R=abs(Z)
	ψ=angle(Z)

    Rrange=range(0,stop=R,length=10)

    #Add unit circle to the axis
    θ=0.0:0.01:2π
	plot!(pl, θ, ones(length(θ)), proj=:polar, lims=(0,1.05), label="", color=:black, showaxis=false)
	scatter!(pl, ϕ, ones(length(p2)), proj=:polar, label="", color=:red, markersize=5)

    #Add DLMO marker
	plot!(pl, 5π/12*ones(11), 0.0:0.1:1.0, proj=:polar, color=:darkgreen, lw=2.0, label="")
	annotate!(pl, 5π/12, 1.02, text("DLMO",8,:center, :darkgreen))

    #Order Parameter
	plot!(pl, ψ*ones(10), Rrange, proj=:polar, color=:red, lw=2.0, label="R=$(round(R, digits=3))")

    return pl

end



"""
function mae_plot(dlmo_actual::AbstractVector, dlmo_lco::AbstractVector; kwargs...)

returns a plot
"""
function mae_plot(dlmo_actual::AbstractVector, dlmo_lco::AbstractVector; kwargs...)

    dlmo_actual .%= 24.0
	dlmo_lco .%= 24.0

	mae1=mean(abs_hour_diff.(dlmo_lco, dlmo_actual))
    #Calculate the percentage under 1 hr 
    errors_under_1hr=sum(abs_hour_diff.(dlmo_lco, dlmo_actual) .<= 1.0)

    println("The MAE (hrs) in DLMO prediction is: $(mae1) for model mae1, with $(errors_under_1hr) under 1 hour error out of $(length(dlmo_actual))")

	#Make the plot range from from -12 to 12
    dlmo_lco=cut_phases_12.(dlmo_lco)
    dlmo_actual=cut_phases_12.(dlmo_actual)

    pl=scatter(dlmo_actual, dlmo_lco, legend=:topleft; kwargs...)
    plot!(pl, -12:1:12, -12:1:12, ls=:dash, color=:gray, lw=2.0, label="")
	ylabel!("Model Prediction (hrs)")
	xlabel!("Experimental DLMO (hrs)")
	return(pl)

end 


"""
function mae_plot!(pl, dlmo_actual::AbstractVector, dlmo_lco::AbstractVector; kwargs...)

returns a plot
"""
function mae_plot!(pl, dlmo_actual::AbstractVector, dlmo_lco::AbstractVector; kwargs...)

    dlmo_actual .%= 24.0
	dlmo_lco .%= 24.0

	mae1=mean(abs_hour_diff.(dlmo_lco, dlmo_actual))
    #Calculate the percentage under 1 hr 
    errors_under_1hr=sum(abs_hour_diff.(dlmo_lco, dlmo_actual) .<= 1.0)

    println("The MAE (hrs) in DLMO prediction is: $(mae1) for model mae1, with $(errors_under_1hr) under 1 hour error out of $(length(dlmo_actual))")

	#Make the plot range from from -12 to 12
    dlmo_lco=cut_phases_12.(dlmo_lco)
    dlmo_actual=cut_phases_12.(dlmo_actual)

    scatter!(pl, dlmo_actual, dlmo_lco, legend=:topleft; kwargs...)
    plot!(pl, -12:1:12, -12:1:12, ls=:dash, color=:gray, lw=2.0, label="")
	ylabel!("Model Prediction (hrs)")
	xlabel!("Experimental DLMO (hrs)")
	return(pl)

end 


