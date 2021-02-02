
#=
This file defines some circular stats functions which are useful when analyzing circadian data.




=#


"""
Find the mean angle for a collection of angles (in radians)
circular_mean(a::AbstractVector)
"""
function circular_mean(a::AbstractVector)

        return angle(sum(1.0/length(a)*exp.(im .*a)))

end

"""
Find the circular correlation coeff between two sets of angles in radians

circular_correlation(a1::AbstractVector, a2::AbstractVector)

https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Circular_Data_Correlation.pdf
"""
function circular_correlation(a1::AbstractVector, a2::AbstractVector)

    @assert all(a1 .< 2π)
    @assert all(a2 .< 2π)

    ψ1=circular_mean(a1)
    ψ2=circular_mean(a2)

    t1=sum(sin.(a1 .- ψ1).^2)
    t2=sum(sin.(a2 .- ψ2).^2)

    rc=sum(sin.(a1 .- ψ1) .* sin.(a2 .- ψ2))/sqrt(t1*t2)

    return rc

end


function meanClockTime(data::AbstractVector)

	#convert to radians
	θ=data .* (1.0/24) .* (2π) 

	ψ=circular_mean(θ)

	if (ψ<0.0)
		ψ+=2π
	end

	return ψ*12.0/π

end


function clocktime_to_radian(d::AbstractVector)

    d_converted=d .* π/12.0

    return d_converted

end


function radian_to_clocktime(d::AbstractVector)

    clocktimes=(d .* 12.0/π) .% 24.0

    for c in eachindex(clocktimes)
        if clocktimes[c]<0.0
            clocktimes[c]+=24.0
        end
    end

    return clocktimes
end




"""
function abs_hour_diff(x,y)

Find the difference in hours between
two clock times (wrapped)
"""
function abs_hour_diff(x,y)
    a1=min(x,y)
    a2=max(x,y)
    s1=a2-a1
    s2=24.0+a1-a2
    return(min(s1,s2))
end



"""
Function to make the branch cut for the DLMO times be at 12 instead of 24.
This is better because lots of DLMOs are near midnight, but many fewer are near
noon.

cut_phases_12(timept)
Can be vectorized using a dot
"""
function cut_phases_12(p)

    while (p<0.0)
        p+=24.0
    end

    p=p % 24.0

    if p>12
        return p-24.0
    else
        return p
    end

end




rad2hr(x)=x*12.0/π
hr2rad(x)=x*π/12.0


"""
    cdist(a, b, hours::Bool=false) -> angle
Return the angular distance from `a` to `b` (b - a) in the forward
direction; hence if `b` is 'behind' `a`, `distance` will be negative.
Angles are confined to be the smaller possible angle, so are in the range
[-π:π], or [-180°:180°].
Angles are in radians, unless `hours` == true.
"""
cdist(a, b) = mod(b - a + pi, 2π) - pi






"""
    cmean(a::Array, hours::Bool=false) -> mean
Return the mean angle from the set of angles in `a`.
Angles are in radians, unless `hours` == true.
"""
cmean(a) = atan(sum(sin.(a)), sum(cos.(a)))
cmean(a, hours::Bool) = hours ? rad2deg(cmean(deg2rad.(a))) : cmean(a)



"""
    cmedian(a::AbstractVector; radians::Bool=false) -> median
Return the median angle from the set of angles `a`.
Angles are in radians, unless `hours` == true. Then we assume
the inputs are 24.0 periodic
"""
function cmedian(a::AbstractVector; hours::Bool=false)

    medians = Vector{Float64}()
    A = Float64.(sort(a))
    n = length(A)
    radians && (A .= clocktime_to_radian.(A))

    # Compute trial bisectors, which are either the points themselves or adjacent means
    p = Array{Float64}(undef, n)
    if iseven(n)
        for i = 1:n
            j = (i + 1 <= n) ? i + 1 : 1
            p[i] = cmean([A[i], A[j]])
        end
    else
        p[:] .= A[:]
    end
    # Try all possible diameters
    for i = 1:n
        # Count points on either side of diameter
        n_plus = sum(cdist.(p[i], A) .> 0)
        n_minus = sum(cdist.(p[i], A) .< 0)
        if n_plus == n_minus
            # Determine which side of circle is correct direction by counting the number
            # within pi/2 of each of the two opposing possible medians
            if sum(abs.(cdist.(p[i], A)) .<= pi/2) > sum(abs.(cdist.(p[i], A)) .> pi/2)
                push!(medians, A[i])
            else
                push!(medians, (A[i] + 2pi)%2pi - pi)
            end
        end
    end
    # If there is more than one median, take the mean thereof (Otieno & Anderson-Cook
    # (2003), Journal of Modern Applied Statistical Methods, 2(1), 168-176)
    median = if length(medians) > 1
        cmean(medians)
    elseif length(medians) == 1
        medians[1]
    else
        error("Zero medians found.  Are data axial but axial!=true?")
    end
    hours && (median = radian_to_clocktime(median))

    return median
end



"""
    cresultant(a, hours=false) -> R
    cresultant(a, w, hours=false) -> Rc
Return the resultant vector length, `R`, from a set of angles, `a`.
If data are binned by binwidth `w`, an unbiased estimate of `R`, `Rc` is returned
when `w` is supplied.
Angles are in radians, unless `hours` == true.
"""
cresultant(a, hours::Bool) = hours ? sqrt(sum(sin.(hr2rad.(a)))^2 + sum(cos.(hr2rad.(a)))^2)/length(a) : sqrt(sum(sin.(a))^2 + sum(cos.(a))^2)/length(a)

#cresultant(a, w::Real; hours::Bool=false) = hours ? cresultant(a, true)*hr2rad(w)/(2sin(hr2rad(w)/2)) : cresultant(a, false)*w/(2*sin(w/2))



"""
    cstd(a, hours=false) -> σ
Return the standard deviation, `σ`, for a set of angles `a`.
Angles are in radians, unless `hours == true`.
"""
cstd(a; hours::Bool=false) = sqrt(-2*log(cresultant(a, hours)))





"""
    cvariance(a, hours=false) -> σ²
Return the circular variance, `σ²`, of a set of angles `a`.
Angles are in radians, unless `hours == true`.
"""
function cvariance(a; hours::Bool=false)
    a_mean = cmean(a, hours)
    if hours
        1 - sum(cos.(hr2rad.(a .- a_mean)))/length(a)
    else
        1 - sum(cos.(a .- a_mean))/length(a)
    end
end
