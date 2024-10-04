include("analysis-functions.jl")
using CSV
using Tables
using DataFrames
using LsqFit
using Dates

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

temps = Vector(range(start=temp_high, stop=temp_low; step=step))

observable_calculation = true
if observable_calculation
    EvT = []
    MvT = []
    CvvT = []
    XvT = []
    YvT = []

    for T in temps
        println("Computing average energy at temperature $T .....")

        energy, magnetization_sq, stiffness_sum_sq = EM_importer("./Data/run$run/meas-$T.h5")

        energy = Vector([energy[i] for i in range(1, length(energy)) if i>EQ_NO && i%MEAS_GAP==0])
        magnetization_sq = Vector([magnetization_sq[i] for i in range(1, length(magnetization_sq)) if i>EQ_NO && i%MEAS_GAP==0])
        stiffness_sum_sq = Vector([stiffness_sum_sq[i] for i in range(1, length(stiffness_sum_sq)) if i>EQ_NO && i%MEAS_GAP==0])

        Eps = avg_energy_ps(energy, L)
        Mps = avg_mag_ps(magnetization_sq, L)
        Cv = specific_heat(energy, T, L)
        Chi = mag_susceptibility(magnetization_sq, T, L)
        Yu = helicity_modulus(stiffness_sum_sq, T, L, mean(energy))

        push!(EvT, Eps)
        push!(MvT, Mps)
        push!(CvvT, Cv)
        push!(XvT, Chi)
        push!(YvT, Yu)

        println("Computation finished.")
    end

    T_diffs = [temps[i+1] - temps[i] for i=1:(length(temps)-1)]
    E_diffs = [EvT[i+1] - EvT[i] for i=1:(length(EvT)-1)]
    derivative = [E_diffs[i] / T_diffs[i] for i=1:(length(EvT) - 1)]
    push!(derivative, 0.0)

    d1 = DataFrame(T = temps, Eps = EvT, Mps = MvT, Cv = CvvT, dEdT = derivative, Chi = XvT, Yu = YvT)
    CSV.write("Analysis/run$run/measurements.csv", d1)
end

