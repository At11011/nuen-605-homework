# Author: Nathaniel Thomas
# Date: August 19th, 2024
# Class: NUEN 604 - 600c

using Plots
using StatsBase
using PrettyTables

unicodeplots()

# Problem 1
data = [133, 132, 131, 134, 135, 137, 131, 137, 134, 136, 138, 129, 133, 133, 134, 134, 137, 132, 134, 131]
occurance_map = countmap(data) # iii
exp_sum = sum(data) # i
println("Sum:\t$exp_sum")
exp_mean = exp_sum/length(data) # ii 
println("Mean:\t$exp_mean")
println("The Gaussian model is appropriate in this case since the mean is more than 25.") # iii
b_range = range(minimum(data), maximum(data) + 1, maximum(data) - minimum(data) + 1)

frequency = collect(values(occurance_map))
table_data = sort(hcat(collect(keys(occurance_map)), frequency), dims = 1)
header = (["Count", "Occurances"])
pretty_table(
    table_data,
    formatters = ft_printf("%5.2f", 2:4),
    header = header,
    header_crayon = crayon"yellow bold",
    tf = tf_unicode_rounded
) # iv

histogram(data, bins = b_range) # iv
gui()
