using GLM, DataFrames, Plots
pgfplotsx()

function linear_regression(x,y)
    df = DataFrame(X = x, Y = y)
    model = lm(@formula(Y ~ 0 + X), df)
    x_range = range(minimum(x), maximum(x), length=100)
    predictions = predict(model, DataFrame(X = x_range))
    return x_range, predictions, coef(model)[1]
end

x = [200,   100,    20,     4   ]
y = [1100,  520,    110,    20  ]
x_range, predictions, slope = linear_regression(x,y)

plot(
    title = "Plot of Pre-Amp Height vs Pulser Pulse Height" *
            "\nfor Attenuation of 1, 2, 10, and 50",
    xaxis = "Pulser Height (mv)",
    yaxis = "Pre-Amp Height (mv)"
)
scatter!(x, y, label = "Data" )
plot!(
    x_range,
    predictions,
    color=:red,
    linestyle=:dash,
    label="Fit (slope = $(round(slope, sigdigits = 2)))", linewidth=1
)
savefig("~/Code/LaTeX/Fall 2024/NUEN 605/NUEN 605 Worksheet 1/plot_3.2.3.tikz")

x = 100:100:1000
y = [0.244, 0.468, 0.704, 0.96, 1.2, 1.42, 1.64, 1.94, 2.14, 2.36]
x_range, predictions, slope = linear_regression(x,y)

plot(
    title="Plot of Peak-to-Peak vs Pulse-Height Setting",
    xaxis = "Pulse-Height Settings (#/1000)",
    yaxis = "Peak-to-Peak (V)"
)
scatter!(x, y, label = "Data")
plot!(
    x_range,
    predictions,
    color=:red,
    linestyle=:dash,
    label="Fit (slope = $(round(slope, sigdigits = 2)))", linewidth=1
)
savefig("~/Code/LaTeX/Fall 2024/NUEN 605/NUEN 605 Worksheet 1/plot_3.2.4b.tikz")
