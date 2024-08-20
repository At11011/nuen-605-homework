# Author: Nathaniel Thomas
# Date: August 19th, 2024
# Class: NUEN 604 - 600c

using Plots

# Problem 1
data = [133, 132, 131, 134, 135, 137, 131, 137, 134, 136, 138, 129, 133, 133, 134, 134, 137, 132, 134, 131]
println("Unique data: ", unique(data))
exp_sum = sum(data) # i
println("Sum:\t$exp_sum")
exp_mean = exp_sum/length(data) # ii 
println("Mean:\t$exp_mean")
println("The Gaussian model is appropriate in this case since the mean is more than 25.")
b_range = range(minimum(data), maximum(data) + 1, maximum(data) - minimum(data) + 1)
histogram(data, bins = b_range) # iii
gui()
