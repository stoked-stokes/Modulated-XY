using Random
using Statistics
using Distributions
using Dumper
using GLMakie
using DataStructures

##### ANALYSIS FUNCTIONS #####

function EM_importer(filename)

    fl = DumpFile(filename)

    vecs = fl["vector"]
    energy = vecs[:,1,:]
    magnetization_sq = vecs[:,2,:]
    stiffness_sum_sq = vecs[:,3,:]

    return energy, magnetization_sq, stiffness_sum_sq
end

function avg_energy_ps(e_series, size)
   
    return mean(e_series) / (size * size)
end

function avg_mag_ps(m_sq_series, size)
   
    return mean(sqrt.(m_sq_series)) / (size * size)
end

function specific_heat(e_series, T, size)

    # println(e_series.^2)
    # println(mean(e_series)^2)
    # println(mean(e_series.^2))

    return (mean(e_series.^2) - mean(e_series)^2) / (T^2*size^4)
end

function mag_susceptibility(m_sq_series, T, size)

    return (mean(m_sq_series) - mean(sqrt.(m_sq_series))^2) / (T*size^4)
end

function mag_accumulant(m_sq_series)

    return mean(m_sq_series.^2) / mean(m_sq_series)
end

function helicity_modulus(stiff_sum_sq_series, T, size, avg_energyps)
    ### measure the helicity modulus (spin stiffness) of the given configuration
    ### config must be a 2D array

    avg_stiff_sum_sq = mean(stiff_sum_sq_series)
    
    return (-0.5) * avg_energyps - (1 / T) * (avg_stiff_sum_sq / size^2)
end

function variance(series)

    r = length(series)
    avg_series = sum(series) / r

    s = 0

    for i = 1:r
        s += (series[i] - avg_series)^2
    end

    return (1 / (r - 1)) * s
end

function autocorrelation_t(series, t)

    points = length(series) - t
    avg_series = sum(series) / length(series)
    sums = 0

    for i = 1:points
        sums += (series[i] - avg_series) * (series[i + t] - avg_series)
    end

    return (sums / points)
end

function integrared_autocorrelation_time(auto_series, C_0, cutoff)

    sum = 0

    for i = 1:cutoff
        sum += auto_series[i]
    end

    return 1 + (2 / C_0) * sum
end

function binning_logarithmic(series, level, cutoff)
    ### To get the original series at the first level, start from level = 0
    ### Works only for length(series) = 2^m, for integer m

    if length(series) < 2^cutoff
        return println("Series is too small for this level of binning")
    end

    if level == cutoff
        return []
    end

    error_series = []

    std_err = std(series) / sqrt(length(series))
    append!(error_series, std_err)

    if length(series) % 2 != 0
        series = series[(begin+1):end]
    end

    series = (series[1:2:end] + series[2:2:end]) / 2
    sub_error_series = binning_logarithmic(series, level+1, cutoff)
    error_series = vcat(error_series, sub_error_series)

    return error_series
end

function binning(series, max_binsize)

    M = length(series)
    errors = []

    for k = 1:max_binsize
        M_k = floor(Int, M / k)
        remainder_k = M % k
        quotient_k = div(M, k)

        binned_series_k = mean(reshape(series[(remainder_k+1):end], :, quotient_k), dims=1)
        
        std_err_k = std(binned_series_k) / sqrt(M_k)
        append!(errors, std_err_k)
    end

    return errors
end

function boltzmann(e_series, T)
    beta = 1 / T
    weights = exp.(- beta .* e_series)

    return weights ./ maximum(weights)
end

##### TEST CODE #####
