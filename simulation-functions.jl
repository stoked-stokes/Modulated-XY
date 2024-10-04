using Random
using Distributions
using Dumper
using DataStructures
using Statistics

##### CONSTANTS #####
J = 1
nn = [(0, 1), (-1, 0), (0, -1), (1, 0)]

##### SIMULATION FUNCTIONS #####

function random_config(size)
    ### returns a random configuration of spins on a lattice of size L x L

    return rand(Uniform(-pi, pi), (size, size))
end

function index_check(i, j, k)
    value = (1 <= i + k[1] <= Ly) && (1 <= j + k[2] <= Lx)

    return value
end

function periodic_chainer(size, index, k)
    if index[1] + k[1] == 0
        return (size, index[2]+k[2])
    elseif index[1] + k[1] == (size+1)
        return (1, index[2]+k[2])
    elseif index[2] + k[2] == 0
        return (index[1]+k[1], size)
    elseif index[2] + k[2] == (size+1)
        return (index[1]+k[1], 1)
    else
        return (index[1]+k[1], index[2]+k[2])
    end
end

function x_lattice_coupling(curr_site, bond, size, base_coupling, alpha, wavelength, phi)
    wave_vec = 2 * pi / wavelength

    if bond[2] == -1
        curr_site = periodic_chainer(size, curr_site, bond)
    end

    return base_coupling + alpha * sin(wave_vec * (curr_site[2] - 1) + phi)
end

function measure_energy(config, base_coupling, alpha, wavelength, phi=0)
    ### measure the energy of the given configuration with periodic boundary conditions
    ### config must be a 2D array

    size = Int(sqrt(length(config)))

    e = 0

    for i= 1:size, j = 1:size
        r = config[i, j]

        for k in nn
            if k[2] == 0
                s = config[periodic_chainer(size, (i, j), k)[1], periodic_chainer(size, (i, j), k)[2]]
                e += -1 * base_coupling * cos(r - s)
            else
                s = config[periodic_chainer(size, (i, j), k)[1], periodic_chainer(size, (i, j), k)[2]]
                e += -1 * x_lattice_coupling([i, j], k, size, base_coupling, alpha, wavelength, phi) * cos(r - s)
            end
        end
    end

    return e / 2
end

function measure_mag_sq(config)
    ### measure the square of norm of magnetization squared of the given configuration
    ### config must be a 2D array

    magx = sum(cos.(config))
    magy = sum(sin.(config))

    mag_sq = magx^2 + magy^2
    return mag_sq
end

function measure_stiffness_sum_sq(config)
    ### measure the helicity modulus (spin stiffness) of the given configuration
    ### config must be a 2D array

    size = Int(sqrt(length(config)))
    arb_unit_vec = (0, 1)

    sum = 0
    for i= 1:size, j = 1:size
        r = config[i, j]

        for k in nn
            s = config[periodic_chainer(size, (i, j), k)[1], periodic_chainer(size, (i, j), k)[2]]
            sum += sin(r - s) * ((arb_unit_vec[1] * k[1]) + (arb_unit_vec[2] * k[2]))
        end
    end
    
    return sum^2
end

function random_spin_update(size)
    ### returns a randomly chosen lattice site and random angle for update

    y_index = rand(1:size)
    x_index = rand(1:size)
    angle = rand(Uniform(-pi, pi))

    return (y_index, x_index), angle
end

function reflection_across_plane(theta, arg)
    ### Reflects the 2D vector with argument theta across the plane orthogonal to the 2D vector with argument arg
    ### Returns the argument of the reflected vector

    reflected_x_coordinate = cos(theta) - 2 * cos(theta - arg) * cos(arg)
    reflected_y_coordinate = sin(theta) - 2 * cos(theta - arg) * sin(arg)

    reflected_angle = atan(reflected_y_coordinate, reflected_x_coordinate)
    return reflected_angle
end

function wolff_update!(current_config, temp, base_coupling, alpha, wavelength, phi=0)
    ### performs one update of the wolff cluster algorithm for markov chain monte carlo

    size = Int(sqrt(length(current_config)))
    beta = 1 / temp

    arg = rand(Uniform(-pi, pi))
    # println("arg: ", arg)
    seed = vec(rand(1:size, (1, 2)))
    # println("seed: ", seed)
    # println("old seed angle: ", current_config[seed[1], seed[2]])
    current_config[seed[1], seed[2]] = reflection_across_plane(current_config[seed[1], seed[2]], arg)
    # println("new seed angle: ", current_config[seed[1], seed[2]])

    unvisited = Deque{Vector{Int64}}()
    push!(unvisited, seed)
    # println("unvisited length: ", length(unvisited))
    cluster_size = 1
    while !isempty(unvisited)
        site = pop!(unvisited)
        # println("site: ", site)
        for k in nn
            nbr = periodic_chainer(size, site, k)
            # println("neighbour: ", nbr)
            # println("old neighbour angle: ", current_config[nbr[1], nbr[2]])

            if k[1] == 0
                alignment = 2 * beta * x_lattice_coupling(site, k, size, base_coupling, alpha, wavelength, phi) * cos(arg - current_config[seed[1], seed[2]]) * cos(arg - current_config[nbr[1], nbr[2]])
            else
                alignment = 2 * beta * base_coupling * cos(arg - current_config[site[1], site[2]]) * cos(arg - current_config[nbr[1], nbr[2]])
            end
            probability = 1 - exp(min(0, alignment))
            # println("probability: ", probability)

            if rand() <= probability
                current_config[nbr[1], nbr[2]] = reflection_across_plane(current_config[nbr[1], nbr[2]], arg)
                # println("new neighbour angle: ", current_config[nbr[1], nbr[2]])
                pushfirst!(unvisited, vec([nbr[1], nbr[2]]))
                cluster_size += 1
            end
        end
        # println("unvisited length: ", length(unvisited))
    end

    return cluster_size
end

##### TEST CODE #####

testing = false

if testing
    Random.seed!(10)

    len = 20
    spins = random_config(len)
    display(spins)

    x_bonds = Vector([Vector{Float64}(undef, len) for _ = 1:len])
    foo = Vector([Vector{Float64}(undef, len) for _ = 1:len])

    for k = 1:len
        for j = 1:len
            x_bonds[k][j] = x_lattice_coupling([k, j], (0, 1), len, 1, 1, len, 0)
            foo[k][j] = 0.0
        end
    end

    tries = 100
    for i = 1:tries
        old_spins = deepcopy(spins)
        x = wolff_update!(spins, 0.4, 1, 1, len)

        for k = 1:len
            for j = 1:len
                foo[k][j] = foo[k][j] + abs(spins[k, j] - old_spins[k, j])
            end
        end
    end

    rung_avg = mean(foo, dims=1)[1]
    println(rung_avg)

    E = measure_energy(spins, 1, 0.5, len)
    println(E / len^2)
end