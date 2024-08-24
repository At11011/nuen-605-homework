# Author: Nathaniel Thomas
# Date: August 19th, 2024
# Class: NUEN 604 - 600c

using Plots
using StatsBase
using PrettyTables
using Measurements
using Printf
using Unitful
using Symbolics
using LsqFit
using DataFrames
using GLM

gr()

function problem1()
    data = [133, 132, 131, 134, 135, 137, 131, 137, 134, 136, 138, 129, 133, 133, 134, 134, 137, 132, 134, 131]

    exp_sum = sum(data)
    exp_mean = exp_sum/length(data)

    header = (["Count", "Occurances", "Probability (F(x))"])
    occurance_map = countmap(data)
    count = collect(keys(occurance_map))
    frequency = collect(values(occurance_map))
    norm_prob = frequency ./ sum(frequency)
    table_data = sortslices(hcat(count, frequency, norm_prob), dims=1, by=row -> row[1])

    println("Problem 1:")
    println("(i) Sum:\t$exp_sum counts")
    println("(ii) Mean:\t$exp_mean counts")
    println("(iii) The Gaussian model is appropriate in this case since the mean is more than 25.")

    pretty_table(
    table_data,
    formatters = ft_printf(["%d", "%d", "%1.2f"], [1, 2, 3]),
    header = header,
    header_crayon = crayon"yellow bold",
    tf = tf_unicode_rounded
    )

    b_range = range(minimum(data), maximum(data) + 1, maximum(data) - minimum(data) + 1)
    bar(count, norm_prob, xtick = minimum(data):maximum(data), title="(iv) Histogram of data", xaxis = "conting measurement", yaxis = "probability", label = "Data")
    vline!([exp_mean], label="Experimental Mean ($exp_mean)", lw=2)
    savefig("./output/homework1/problem_1_plot.svg")

    calc_exp_mean = sum(norm_prob .* count)
    println("(v) Calculated experimental mean (x̄ = ∑x⋅F(x)):\t$calc_exp_mean")
    println("The calculated mean ", calc_exp_mean == exp_mean ? "does" : "does not"," match the experimental mean.")
    println("(vi) The data seems to roughly match a Gaussian distribution, but it is a very rough fit. ",
    "There is a peak around the experimental mean and relatively small tails. ",
    "To make the distribution smoother, more measurements are needed.")
end


function problem2()
    data = [12548, 12725, 12684, 12987, 12784]
    # histogram(data)
    # gui()
    println("\nProblem 2:")
    println("(i) In this case, since the mean is greater than 25, the Gaussian model is best.")
    println("(ii) When considering each measurement individually, the standard deviation can be calculated by taking"
                * " the square root of the mean, in this case, the value of each measurement.")
    std_dev = sqrt.(data)
    totals = data .± std_dev
    println("(iv) Values ± Standard deviations:")
    [println("\t$(@sprintf("%d", Measurements.value(val))) ± $(@sprintf("%1.2f", Measurements.uncertainty(val))) counts") for val in totals]
    
    duration = 120u"s"
    count_rate = totals ./ duration

    println("(iv) Count rates:")
    println(count_rate)
    [println("\t$(@sprintf("%d", Measurements.value(ustrip(val)))) ± $(@sprintf("%1.2f", Measurements.uncertainty(ustrip(val)))) counts/s") for val in count_rate]

    return nothing
end

function problem3()
    @variables x, y, z, σx, σy, σz
    
    function error_prop(u, subx, suby, subz)
        du²dx² = (Differential(x)(u))^2 
        du²dy² = (Differential(y)(u))^2
        du²dz² = (Differential(z)(u))^2
        println("\tdu²dx²:\t\t$(expand_derivatives(substitute(du²dx², Dict([x => subx, y => suby, z => subz]))))")
        println("\tdu²dy²:\t\t$(expand_derivatives(substitute(du²dy², Dict([x => subx, y => suby, z => subz]))))")
        println("\tdu²dz²:\t\t$(expand_derivatives(substitute(du²dz², Dict([x => subx, y => suby, z => subz]))))")
        return du²dx² * σx^2 +
               du²dy² * σy^2 +
               du²dz² * σz^2
    end

    println("\nProblem 3:")
    println("(i):")
    u = x - y
    println("\tExpression:\tu = $u")
    σᵤ = error_prop(u,x, y, z)
    println("\tError Equation:\tσᵤ = √($(expand_derivatives(σᵤ)))")

    println("(ii):")
    u = x / y
    println("\tExpression:\tu = $u")
    σᵤ = error_prop(u, x, y, z)
    println("\tError Equation:\tσᵤ = √($(expand_derivatives(σᵤ)))")

    @variables k
    println("(iii):")
    u = k*x 
    println("\tExpression:\tu = $u")
    σᵤ = error_prop(u, x, y, z)
    println("\tError Equation:\tσᵤ = √($(expand_derivatives(σᵤ)))")
    
    @variables A₀, λ, t
    println("(iv):")
    u = x*exp(-y*z) 
    println("\tExpression:\tA = $(substitute(u, Dict([x => A₀, y => λ, z => t])))")
    σᵤ = error_prop(u, A₀, λ, t)
    println("\tError Equation:\tσA = √($(expand_derivatives(σᵤ)))")
    
    @variables M₁, M₂, M₁₂
    println("(iv):")
    u = (x + y - z) /(2*x*y) 
    println("\tEquation:\tτ = $(substitute(u, Dict([x => M₁, y => M₂, z => M₁₂])))")
    σᵤ = error_prop(u, M₁, M₂, M₁₂)
    println("\tError Equation:\tστ = √($(expand_derivatives(σᵤ)))")
end 

function problem4()
    distance = (10:10:100)u"mm"
    counts = [2071, 1611, 1241, 1040, 856, 652, 599, 492, 443, 386]
    error = sqrt.(counts)
    totals = counts .± error
    time = 30u"s"
    println("\nProblem 4")
    println("(i) The inverse square law states that the intensity of a source is inversely proportional to the square of the distance from the source.")
    count_rate = totals ./ time
    scatter(
        distance,
        count_rate,
        xaxis = "Distance",
        yaxis = "Count rate (count/s)",
        title = "(ii) Count rate vs distance",
        label = "Data",
        markershape = :circle
    )
    # Develop a fit
    x = 1 ./ (ustrip(distance).^ 2)
    y = ustrip(Measurements.value.(count_rate))
    model = lm(@formula(y ~ x), DataFrame(x = x, y = y))
    offset = coef(model)[1]u"s^-1"
    coefficient = coef(model)[2]u"mm^2*s^-1"
    
    r_squared = r2(model)
    fit_data = offset .+ coefficient ./ distance.^2
    plot!(
        distance,
        fit_data,
        label = "Fit\ny = $(@sprintf "%0.3f + %0.3f / x²" ustrip(offset) ustrip(coefficient))\n(R² = $(@sprintf "%0.5f" r_squared))",
    )

    savefig("./output/homework1/problem_4_plot.svg")
    println("(iii) Based on this fit, the data does not follow the inverse square law.")
end 

function problem5()
    prob = 0.7
    total = 34
    partial = 28
    println("\nProblem 5:")
    println("(i) In this case, the experimental mean is unknown, and the success/failure is determined by a successful or unsuccessful measurement. This type of binary criteria lends itself to the binomial distribution.")
    part_prob = factorial(big(total), big(total - partial)) / factorial(big(partial)) * prob ^ partial * (1 - prob) ^ (total - partial)
    println("(ii) The probability is $(@sprintf "%0.5f" part_prob)")
    avg_success = prob * total
    std_dev = sqrt(avg_success * (1 - prob))
    println("(iii) The predicted average number of successes is $(@sprintf "%0.1f. The standard deviation is %0.5f." avg_success std_dev)")
end

function problem6()
    time = 60u"s"
    count_back = 81
    count_total = 567
    count_source = count_total - count_back
    println("\nProblem 6:") 
    println("(i) Net count:\t$count_source")
    println("(ii) Count rates:\n\tBackground:\t$(ustrip(count_back / time)) counts/s\n\tWith source:\t$(ustrip(count_total/time)) counst/s")
    count_back_unc = sqrt(count_back)
    count_total_unc = sqrt(count_total)
    count_source_unc = sqrt(count_back_unc^2 + count_total_unc^2)
    println("(iii) Net Count Rate:\t$((count_source ± count_source_unc) / ustrip(time)) counts/s")
end

function problem7()
    count = (1/(0.05/100))^2
    frac_std = sqrt(count)/count * 100
    println("\nProblem 7:")
    println("(i) $count counts are needed for $(@sprintf "%0.3f" frac_std)% fractional standard deviation.")
end

function problem8()
    count1 = 897
    count2 = 920
    count1_err = sqrt(count1)
    count2_err = sqrt(count2)
    average = (count1 + count2)/2

    expr = (x + y) / 2
    println("\nProblem 8:")
    println("(i) Average and Error:\t$(average ± sqrt(count1_err^2 + count2_err^2)/2)")
end

function problem9()
    yield = 85u"percent"
    time = 60u"s"
    total_rate = 360.06u"s^-1"
    source_rate = 279.85u"s^-1"
    background_rate = (total_rate * time - source_rate * time) / time
    background_rate = background_rate ± √background_rate * 1u"s^-(1/2)"
    total_rate = total_rate ± √total_rate* 1u"s^-(1/2)"
    source_rate = total_rate - background_rate

    A₀ = (2.74e5 ± 2.74e4)u"Bq"
    λ = 0.74/24u"hr"
    t = 6*24u"hr" 
    A = A₀*exp(-λ*t)
    efficiency = source_rate / (yield * A)    

    println("Problem 9:")
    println("(i):\n\tActivity:\t$A")
    println("\tEfficiency:\t$(uconvert(u"percent", efficiency))")

    println("(ii)")

    efficiency = 25u"percent"
    Lc = 2.33 * Measurements.uncertainty(background_rate * time)
    min_count_rate = (total_rate * time - background_rate * time) /time
    min_activity = min_count_rate / (yield * efficiency)
    println("\tMinimum Detectable Activity: $(uconvert(u"Bq", min_activity))")
end

problem1()
problem2()
problem3()
problem4()
problem5()
problem6()
problem7()
problem8()
problem9()
