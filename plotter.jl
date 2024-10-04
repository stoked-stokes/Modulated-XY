include("analysis-functions.jl")
using GLMakie
using CSV
using DataFrames

run = 20
L = 16
EQ_NO = 500 * L

##### PLOTTING #####

if ARGS[1] == "direct-binning"

    d = CSV.read("./Analysis/run$run/direct-binning/tau-T.csv", DataFrame)

    f = Figure(size = (800, 500))

    ax1 = Axis(f[1, 1], xlabel = "Temperature", ylabel = "Autocorrelation time")
    scatter!(ax1, d.T, d.tau_E, color = :red, label = "τᴱ")
    axislegend()

    ax2 = Axis(f[1, 2], xlabel = "Temperature", ylabel = "Autocorrelation time")
    scatter!(ax2, d.T, d.tau_M, color = :blue, label = "τᴹ")
    axislegend()

    save("./Plots/run$run/tau-integrated/tau-T-direct.png", f)

elseif ARGS[1] == "full-binning"

    d = CSV.read("./Analysis/run$run/full-binning/tau-T.csv", DataFrame)

    f = Figure(size = (800, 500))

    ax1 = Axis(f[1, 1], xlabel = "Temperature", ylabel = "Autocorrelation time")
    scatter!(ax1, d.T, d.taus_E, color = :red, label = "τᴱ")
    axislegend()

    ax2 = Axis(f[1, 2], xlabel = "Temperature", ylabel = "Autocorrelation time")
    scatter!(ax2, d.T, d.taus_M, color = :blue, label = "τᴹ")
    axislegend()

    save("./Plots/run$run/tau-integrated/tau-T-full.png", f)

    for T = 1.5:1.5
        d2 = CSV.read("./Analysis/run$run/full-binning/binsize/$T-tau-binsize.csv", DataFrame)

        g = Figure(size = (800, 500))

        ax3 = Axis(g[1, 1], xlabel = "Binsize", ylabel = "Autocorrelation time")
        scatter!(ax3, d2.binsize, d2.taus_E, color = :red, label = "τᴱ")
        axislegend()

        ax4 = Axis(g[1, 2], xlabel = "Binsize", ylabel = "Autocorrelation time")
        scatter!(ax4, d2.binsize, d2.taus_M, color = :blue, label = "τᴹ")
        axislegend()

        save("./Plots/run$run/tau-integrated/full-tau-binsize/$T-tau-binsize.png", g)
    end

elseif ARGS[1] == "log-binning"

    d = CSV.read("./Analysis/run$run/log-binning/tau-T.csv", DataFrame)

    f = Figure(size = (800, 500))

    ax1 = Axis(f[1, 1], xlabel = "Temperature", ylabel = "Autocorrelation time")
    scatter!(ax1, d.T, d.taus_E, color = :red, label = "τᴱ")
    axislegend()

    ax2 = Axis(f[1, 2], xlabel = "Temperature", ylabel = "Autocorrelation time")
    scatter!(ax2, d.T, d.taus_M, color = :blue, label = "τᴹ")
    axislegend()

    save("./Plots/run$run/tau-integrated/tau-T-log.png", f)

    for T = 1.5:1.5
        d2 = CSV.read("./Analysis/run$run/log-binning/binsize/$T-tau-binsize.csv", DataFrame)

        g = Figure(size = (800, 500))

        ax3 = Axis(g[1, 1], xlabel = "Binsize", ylabel = "Autocorrelation time")
        scatter!(ax3, d2.binsize, d2.taus_E, color = :red, label = "τᴱ")
        axislegend()

        ax4 = Axis(g[1, 2], xlabel = "Binsize", ylabel = "Autocorrelation time")
        scatter!(ax4, d2.binsize, d2.taus_M, color = :blue, label = "τᴹ")
        axislegend()

        save("./Plots/run$run/tau-integrated/log-tau-binsize/$T-tau-binsize.png", g)
    end

elseif ARGS[1] == "jackknife"

elseif ARGS[1] == "boltzmann"

    # for T = 0.1:0.1
    #     d3 = CSV.read("./Analysis/run$run/measurement/distribution/$T-boltzmann-distribution.csv", DataFrame)

    #     h = Figure(size = (800, 500))

    #     ax5 = Axis(h[1, 1], xlabel = "Energy", ylabel = "Number of entries")
    #     hist!(ax5, log.(d3.weight).*(-T), bins = 15, color = :red, label = "Values")
    #     axislegend()

    #     save("./Plots/run$run/boltzmann-distribution/$T-number-energy.png", h)
    # end

    for T = 4.0:4.0
        energy, magnetization_sq, stiffness_sum_sq = EM_importer("./Data/run$run/meas-$T.h5")
        energy = energy[(EQ_NO+1):end]

        h = Figure(size = (800, 500))

        ax5 = Axis(h[1, 1], xlabel = "Energy", ylabel = "Number of entries")
        hist!(ax5, energy/L^2, color = :red, label = "Values")
        axislegend()

        save("./Plots/run$run/boltzmann-distribution/$T-number-energy.png", h)
    end

end