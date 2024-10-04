include("analysis-functions.jl")
using CSV
using Tables
using DataFrames
using LsqFit
using Dates
using OnlineLogBinning
using PyPlot

##### CONSTANTS #####

run = 1
update_type = "Wolff" # or Glauber or Swendsen-Wang
start = "Sequential" # or Random (new hot start at every T)
L = 50
temp_high = 2.5
temp_low = 0.1
step = -0.1
EQ_NO = 1000 * L
MEAS_NO = 10000
MEAS_GAP = 1
J0 = 1
alpha = 0.5
lambda = L

##### SET-UP #####

if !isdir("./Analysis/run$run")
    dirpath = mkdir("./Analysis/run$run")
    filename = "analysis-features.txt"
    fpath = joinpath(dirpath, filename)
    open(fpath, "w") do file
        write(file, "Measurement gap: $MEAS_GAP")
    end
end

##### MAIN CODE #####

custom = true
if custom
    test_temp = 4.0
    max_binsize = 50

    energy, magnetization_sq, stiffness_sum_sq = EM_importer("./Data/run$run/meas-$test_temp.h5")

    energy = Vector([energy[i] for i in range(1, length(energy)) if i>EQ_NO]) # could also use condition && i%50==0
    magnetization_sq = Vector([magnetization_sq[i] for i in range(1, length(magnetization_sq)) if i>EQ_NO])
    stiffness_sum_sq = Vector([stiffness_sum_sq[i] for i in range(1, length(stiffness_sum_sq)) if i>EQ_NO])
    
    energy_errors = binning(energy, max_binsize)
    magnetization_errors = binning(sqrt.(magnetization_sq), max_binsize)
    stiffness_sum_sq = binning(sqrt.(stiffness_sum_sq), max_binsize)
    # energy_errors = binning_logarithmic(energy, 0, 3)
    # magnetization_errors = binning_logarithmic(sqrt.(magnetization_sq), 0, 3)
    # stiffness_sum_sq = binning_logarithmic(sqrt.(stiffness_sum_sq), 0, 3)

    bin_levels = Vector(range(1, max_binsize))

    # d3 = DataFrame(b = bin_levels, E_err = energy_errors, M_err = magnetization_errors, Y_err = stiffness_sum_sq)
    # CSV.write("Analysis/run$run/binning_errors_$test_temp.csv", d3)
end

bacc = BinningAccumulator()
println(push!(bacc, energy))
# println(bacc[level = 0])
# println(bacc[level = 1])
# println(bacc[level = 2])
# println(bacc[level = 3])
# println(bacc[level = 4])
# println(bacc[level = 5])

result = fit_RxValues(bacc)
println(result)
println(autocorrelation_time(result))

errors = [std_error(bacc; level = i) for i = 0:12]
levels = 0:12

x = plot(levels, errors, linestyle="", marker="x")
savefig("./Plots/run17/binning-errors/packagebinning_$test_temp.png")


