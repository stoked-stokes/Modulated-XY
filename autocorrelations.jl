include("analysis-functions.jl")
using BinningAnalysis
using BinningAnalysis: correlation, unbinned_tau
using GLMakie
using CSV
using DataFrames
using LsqFit
using Dates

##### CONSTANTS #####

run = 28
update_type = "Wolff" # or Glauber or Swendsen-Wang
start = "Sequential" # or Random (new hot start at every T)
L = 50
temp_high = 1.8
temp_low = 0.05
step = -0.05
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

if ARGS[1] == "direct"

    if !isdir("./Analysis/run$run/direct-binning")
        dirpath = mkdir("./Analysis/run$run/direct-binning")
        dirpath = mkdir("./Analysis/run$run/direct-binning/correlations")
    end

    M = 1000
    tau_EvT = []
    tau_MvT = []

    for T in temps
        energy, magnetization_sq, stiffness_sum_sq = EM_importer("./Data/run$run/meas-$T.h5")
        energy = energy[(EQ_NO+1):end]
        magnetization = (sqrt.(magnetization_sq[(EQ_NO+1):end]))

        println("Computing autocorrelations for $T K .....")
        correlations_E = [correlation(energy, k) for k in 0:M]
        correlations_M = [correlation(magnetization, k) for k in 0:M]
        println("Computation finished")

        dT = DataFrame(index = 1:length(correlations_E), correlationsE = correlations_E, correlationsM = correlations_M)
        CSV.write("./Analysis/run$run/direct-binning/correlations/$T-chi.csv", dT)

        push!(tau_EvT, unbinned_tau(correlations_E))
        push!(tau_MvT, unbinned_tau(correlations_M))
    end

    d = DataFrame(T = temps, taus_E = tau_EvT, taus_M = tau_MvT)
    CSV.write("./Analysis/run$run/direct-binning/tau-T.csv", d)

elseif ARGS[1] == "full"

    if !isdir("./Analysis/run$run/full-binning")
        dirpath = mkdir("./Analysis/run$run/full-binning")
        dirpath = mkdir("./Analysis/run$run/full-binning/binsize")
    end

    tau_EvT = []
    tau_MvT = []

    for T in temps
        energy, magnetization_sq, stiffness_sum_sq = EM_importer("./Data/run$run/meas-$T.h5")
        energy = energy[(EQ_NO+1):end]
        magnetization = (sqrt.(magnetization_sq[(EQ_NO+1):end]))

        println("Binning data for $T K .....")
        FBE = FullBinner(energy)
        taus_E = all_taus(FBE)
        FBM = FullBinner(magnetization)
        taus_M = all_taus(FBE)
        println("Binning finished")

        dT = DataFrame(binsize = 1:length(taus_E), taus_E = taus_E, taus_M = taus_M)
        CSV.write("./Analysis/run$run/full-binning/binsize/$T-tau-binsize.csv", dT)

        push!(tau_EvT, tau(FBE))
        push!(tau_MvT, tau(FBM))
    end

    d = DataFrame(T = temps, taus_E = tau_EvT, taus_M = tau_MvT)
    CSV.write("./Analysis/run$run/full-binning/tau-T.csv", d)

elseif ARGS[1] == "log"

    if !isdir("./Analysis/run$run/log-binning")
        dirpath = mkdir("./Analysis/run$run/log-binning")
        dirpath = mkdir("./Analysis/run$run/log-binning/binsize")
    end

    tau_EvT = []
    tau_MvT = []

    for T in temps
        energy, magnetization_sq, stiffness_sum_sq = EM_importer("./Data/run$run/meas-$T.h5")
        energy = energy[(EQ_NO+1):end]
        magnetization = (sqrt.(magnetization_sq[(EQ_NO+1):end]))

        println("Binning data for $T K .....")
        LBE = LogBinner(energy)
        taus_E = all_taus(LBE)
        LBM = LogBinner(magnetization)
        taus_M = all_taus(LBM)
        println("Binning finished")

        dT = DataFrame(binsize = 1:length(taus_E), taus_E = taus_E, taus_M = taus_M)
        CSV.write("./Analysis/run$run/log-binning/binsize/$T-tau-binsize.csv", dT)

        push!(tau_EvT, tau(LBE))
        push!(tau_MvT, tau(LBM))
    end

    d = DataFrame(T = temps, taus_E = tau_EvT, taus_M = tau_MvT)
    CSV.write("./Analysis/run$run/log-binning/tau-T.csv", d)

elseif ARGS[1] == "jackknife"

    if !isdir("./Analysis/run$run/jackknife")
        dirpath = mkdir("./Analysis/run$run/jackknife")
    end

end