include("simulation-functions.jl")
using Random
using Dates

##### CONSTANTS #####

run = 27
update_type = "Wolff" # or Glauber or Swendsen-Wang
start = "Random" # or Random (new hot start at every T)
L = 10
temp_high = 2.5
temp_low = 0.05
step = -0.05
EQ_NO = 100 * L
MEAS_NO = 1000

##### SET-UP #####

dirpath = mkdir("./Data/run$run")
filename = "simulation-features.txt"
fpath = joinpath(dirpath, filename)
open(fpath, "w") do file
    write(file, "System size: $L\n")
    write(file, "Update type: $update_type\n")
    write(file, "Start: $start\n")
    write(file, "Temperature, high: $temp_high\n")
    write(file, "Temperature, low: $temp_low\n")
    write(file, "Temperature, stepsize: $step\n")
    write(file, "Equilibriation steps: $EQ_NO\n")
    write(file, "Measurement steps: $MEAS_NO")
end

##### SIMULATION #####

println(Dates.Time(Dates.now()))

temps = Vector(range(start=temp_high, stop=temp_low; step=step))

for T in temps
    println("Temperature = ", T)
    spins = random_config(L)

    clus_fl = DumpFile("./Data/run$run/clus-$T.h5")
    meas_fl = DumpFile("./Data/run$run/meas-$T.h5")
    config_fl = DumpFile("./Data/run$run/config-$T.h5")

    println("Equilibriation stage started .....")

    for i = 1:EQ_NO
        cluster_size = wolff_update!(T, spins)

        E = measure_energy(spins)
        M = measure_mag_sq(spins)
        Y = measure_stiffness_sum_sq(spins)
        vec = [E, M, Y]

        dump!(clus_fl, "cluster", cluster_size)
        dump!(meas_fl, "vector", vec)
        dump!(config_fl, "matrix", spins)

        println("Equilibrium update $i complete .....")
    end

    println("Equilibriation stage ended ($EQ_NO steps)")

    println("Measurement stage started .....")

    for i = 1:MEAS_NO
        cluster_size = wolff_update!(T, spins)

        E = measure_energy(spins)
        M = measure_mag_sq(spins)
        Y = measure_stiffness_sum_sq(spins)
        vec = [E, M, Y]

        dump!(clus_fl, "cluster", cluster_size)
        dump!(meas_fl, "vector", vec)
        dump!(config_fl, "matrix", spins)

        println("Measurement update $i complete .....")
    end

    println("Measurement stage ended ($MEAS_NO steps)")
end

println(Dates.Time(Dates.now()))