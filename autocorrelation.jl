include("analysis-functions.jl")
using CSV
using Tables
using DataFrames
using LsqFit
using Dates
using BinningAnalysis

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

temps = Vector(range(start=temp_high, stop=temp_low; step=step))
# temps = Vector([0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5])

autocorrelation_calculation = false
tau_cutoff = 100

if autocorrelation_calculation
    println(Dates.Time(Dates.now()))

    tau_integrated_EvT = []
    tau_integrated_MvT = []

    for T in temps
        energy, magnetization_sq, stiffness_sum_sq = EM_importer("./Data/run$run/meas-$T.h5")

        energy = energy[(EQ_NO+1):end]./(L*L)
        magnetization = (sqrt.(magnetization_sq[(EQ_NO+1):end]))./(L*L)

        println("Computing autocorrelations for $T K .....")
        autocorrelation_E_0 = autocorrelation_t(energy, 0)
        autocorrelation_M_0 = autocorrelation_t(magnetization, 0)
        autocorrelation_E = [autocorrelation_t(energy, t) for t = 1:(length(energy)-1)]
        autocorrelation_M = [autocorrelation_t(magnetization, t) for t = 1:(length(magnetization)-1)]
        println("Computation finished")

        tau_integrated_E = integrared_autocorrelation_time(abs.(autocorrelation_E), autocorrelation_E_0, tau_cutoff)
        tau_integrated_M = integrared_autocorrelation_time(abs.(autocorrelation_M), autocorrelation_M_0, tau_cutoff)

        push!(tau_integrated_EvT, tau_integrated_E)
        push!(tau_integrated_MvT, tau_integrated_M)
    end

    d2 = DataFrame(T = temps, tau_int_E = tau_integrated_EvT, tau_int_M = tau_integrated_MvT)
    CSV.write("Analysis/run$run/autocorrelation_time.csv", d2)

    println(Dates.Time(Dates.now()))
end