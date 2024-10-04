import h5py
import numpy as np
import matplotlib.pyplot as plt
import csv

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
points = int(((temp_high - temp_low)/(-step))+1)
temps = np.linspace(temp_high, temp_low, points)
J0 = 1
alpha = 0.5
lamda = L

##### TIME SERIES #####

energy_var = np.array([])

time_series = False
if time_series:
    for T in temps:
        T = round(T, 2)
        filename = f"./Data/run{run}/meas-{T}.h5"

        print(f"Importing time series data for T = {T} ...")
        with h5py.File(filename, "r") as file:
            vecs = file["vector"]
            energy = np.array(vecs[:,0])
            magnetization_sq = np.array(vecs[:,1])
            stiffness_sum_sq = np.array(vecs[:,2])

            file.close()
        magnetization = np.sqrt(magnetization_sq)
        print("Import finished ...")

        energy_var = np.append(energy_var, np.mean(np.square(energy[(EQ_NO+1):] / L**2)) - np.mean(energy[(EQ_NO+1):] / L**2)**2)
        
        print(f"Plotting time series data for T = {T} ...")
        time_domain = range(1, len(energy)+1)

        fig, ax = plt.subplots()
        ax.plot(time_domain, energy / L**2, marker='x', color='red', label='Energy per spin')
        ax.axvline(x=EQ_NO, color='black', linestyle='--')
        ax.set_title(f'T = {T}')
        ax.set_ylabel('Energy')
        ax.set_xlabel('Monte Carlo time (steps)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./Plots/run{run}/time-series/Evt@{T}.png')
        # plt.show()

        # fig, ax = plt.subplots()
        # ax.plot(time_domain, energy**2 / L**4, marker='x', color='red', label='Energy squared per spin')
        # ax.axvline(x=EQ_NO, color='black', linestyle='--')
        # ax.set_title(f'T = {T}')
        # ax.set_ylabel(r'E^2')
        # ax.set_xlabel('Monte Carlo time (steps)')
        # plt.legend()
        # plt.tight_layout()
        # plt.savefig(f'./Plots/run{run}/energy^2/Evt@{T}.png')
        # plt.show()

        fig, ax = plt.subplots()
        ax.plot(time_domain, magnetization / L**2, color='blue', marker='x', label='Magnetization per spin')
        ax.axvline(x=EQ_NO, color='black', linestyle='--')
        ax.set_title(f'T = {T}')
        ax.set_ylabel('Magnetization')
        ax.set_xlabel('Monte Carlo time (steps)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./Plots/run{run}/time-series/Mvt@{T}.png')
        # plt.show()

        fig, ax = plt.subplots()
        ax.plot(time_domain, stiffness_sum_sq, color='green', marker='x', label='Stiffness sum')
        ax.axvline(x=EQ_NO, color='black', linestyle='--')
        ax.set_title(f'T = {T}')
        ax.set_ylabel('Stiffness sum')
        ax.set_xlabel('Monte Carlo time (steps)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./Plots/run{run}/time-series/Yvt@{T}.png')
        # plt.show()

        print("Plotting finished ...")

##### MEASUREMENTS #####

# print(energy_var)

measurements = True
if measurements:
    T = np.array([])
    Eps = np.array([])
    Mps = np.array([])
    Cv = np.array([])
    dEdT = np.array([])
    Chi = np.array([])
    Yu = np.array([])

    print("Importing measurement data ...")
    with open(f'Analysis/run{run}/measurements.csv', 'r') as csvfile:
        csv_reader = csv.reader(csvfile)
        headers = next(csvfile)

        for row in csv_reader:
            T = np.append(T, float(row[0]))
            Eps = np.append(Eps, float(row[1]))
            Mps = np.append(Mps, float(row[2]))
            Cv = np.append(Cv, float(row[3]))
            dEdT = np.append(dEdT, float(row[4]))
            Chi = np.append(Chi, float(row[5]))
            Yu = np.append(Yu, float(row[6]))
    print("Importing finished ...")

    # fit_E = True
    # if fit_E:
        

    print("Plotting measurement data ...")
    fig, ax = plt.subplots()
    ax.plot(T, Eps, color='red', marker='x', label='Energy per spin')
    ax.set_title(f'{EQ_NO} steps + {MEAS_NO} steps')
    ax.set_ylabel(r'$E$')
    ax.set_xlabel('Temperature')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./Plots/run{run}/against-temp/EvT.png')
    # plt.show()

    fig, ax = plt.subplots()
    ax.plot(T, Mps, color='blue', marker='x', label='Magnetization per spin')
    ax.set_title(f'{EQ_NO} steps + {MEAS_NO} steps')
    ax.set_ylabel(r'$|M|$')
    ax.set_xlabel('Temperature')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./Plots/run{run}/against-temp/MvT.png')
    # plt.show()

    fig, ax = plt.subplots()
    ax.plot(T, Cv, color='green', marker='x', label='Specific heat')
    # ax.set_yscale('log')
    ax.set_title(f'{EQ_NO} steps + {MEAS_NO} steps')
    ax.set_ylabel(r'$C_v$')
    ax.set_xlabel('Temperature')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./Plots/run{run}/against-temp/CvvT.png')
    # plt.show()

    # new_step = np.abs(step) / 2
    # T_new = np.array([(i + new_step) for i in temps[:-1]])

    # fig, ax = plt.subplots()
    # ax.plot(T_new, dEdT[:-1], color='purple', marker='x', label='Derivative of Energy')
    # # ax.plot(T, Cv*L**2*T, color='green', marker='x', label='Specific heat')
    # # ax.set_yscale('log')
    # ax.set_title(f'{EQ_NO} steps + {MEAS_NO} steps')
    # ax.set_ylabel(r'$\frac{dE}{dT}$')
    # ax.set_xlabel('Temperature')
    # plt.legend()
    # plt.tight_layout()
    # plt.savefig(f'./Plots/run{run}/against-temp/dEdT.png')
    # # plt.show()

    fig, ax = plt.subplots()
    ax.plot(T, Chi, color='black', marker='x', label='Magnetic Susceptibility')
    ax.set_title(f'{EQ_NO} steps + {MEAS_NO} steps')
    ax.set_ylabel(r'$\chi_M$')
    ax.set_xlabel('Temperature')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./Plots/run{run}/against-temp/XvT.png')
    # plt.show()

    fig, ax = plt.subplots()
    ax.plot(T, Yu / L**2, color='orange', marker='x', label='Helicity Modulus')
    ax.plot(T, (2/np.pi)*T, color='brown', linestyle='--')
    ax.set_ylim(0.0, 1.0)
    ax.set_title(f'{EQ_NO} steps + {MEAS_NO} steps')
    ax.set_ylabel(r'$\Upsilon$')
    ax.set_xlabel('Temperature')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./Plots/run{run}/against-temp/YvT.png')
    # plt.show()
    print("Plotting finished ...")

##### CLUSTER SIZES #####

clusters = False
if clusters:
    avg_cluster_size = np.array([])
    error_avg_cluster_size = np.array([])

    for T in temps:
        T = round(T, 2)
        filename = f"./Data/run{run}/clus-{T}.h5"

        print(f"Importing cluster size data for T = {T} ...")
        with h5py.File(filename, "r") as file:
            clus = file["cluster"]
            cluster_sizes = np.array(clus[:])

            file.close()
        print("Import finished ...")

        avg_cluster_size = np.append(avg_cluster_size, np.mean(cluster_sizes))
        error_avg_cluster_size = np.append(error_avg_cluster_size, np.std(cluster_sizes))

        print(f"Plotting clsuter size data for T = {T} ...")
        time_domain = range(1, len(cluster_sizes)+1)

        fig, ax = plt.subplots()
        ax.plot(time_domain, cluster_sizes, marker='x', color='red', label='Cluster size')
        ax.axvline(x=EQ_NO, color='black', linestyle='--')
        ax.set_title(f'T = {T}')
        ax.set_ylabel('Cluster size')
        ax.set_xlabel('Monte Carlo time (steps)')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'./Plots/run{run}/cluster-sizes/@{T}.png')
        # plt.show()
    
    fig, ax = plt.subplots()
    ax.plot(temps, avg_cluster_size, marker='x', color='blue', label='Avg. c')
    # ax.errorbar(temps, avg_cluster_size, yerr = error_avg_cluster_size, marker='o', color='blue', capsize=2.0, label='Avg. c')
    # ax.axvline(x=EQ_NO, color='black', linestyle='--')
    # ax.set_title(f'T = {T}')
    ax.set_ylabel('Average cluster size')
    ax.set_xlabel('Temperature')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./Plots/run{run}/cluster-sizes/avg_cluster_size.png')
    # plt.show()

##### AUTOCORRELATION TIME #####

autocorrelation_time = False
if autocorrelation_time:
    T = np.array([])
    tau_int_E = np.array([])
    tau_int_M = np.array([])

    print("Importing measurement data ...")
    with open('Analysis/autocorrelation_time.csv', 'r') as csvfile:
        csv_reader = csv.reader(csvfile)
        headers = next(csvfile)

        for row in csv_reader:
            T = np.append(T, float(row[0]))
            tau_int_E = np.append(tau_int_E, float(row[1]))
            tau_int_M = np.append(tau_int_M, float(row[2]))
    print("Importing finished ...")

    print("Plotting autocorrelation times ...")
    fig, ax = plt.subplots()
    ax.plot(T, tau_int_E, color='blue', marker='x', label=r'$\tau_{int}$ for E')
    ax.plot(T, tau_int_M, color='orange', marker='^', label=r'$\tau_{int}$ for M')
    ax.set_title(f'{EQ_NO} steps + {MEAS_NO} steps')
    ax.set_ylabel('Autocorrelation time')
    ax.set_xlabel('Temperature')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./Plots/run{run}/tau-integrated/tau-int.png')
    # plt.show()
    print("Plotting finished ...")

# l_domain = range(1, (length(energy)-1))

# plot_normalized_autocorrelation_graphs = false
# if plot_normalized_autocorrelation_graphs
#     z = plot(l_domain, autocorrelation_E./autocorrelation_E_0, marker="x", color="red", label="Normalized autocorrelation, E")
#     legend()
#     savefig("./SimulationVisualize/run$run/autocorrelations/ACEvt@$T.png")
#     close("all")

#     z = plot(l_domain, autocorrelation_M./autocorrelation_M_0, marker="x", color="blue", label="Normalized autocorrelation, M")
#     legend()
#     savefig("./SimulationVisualize/run$run/autocorrelations/ACMvt@$T.png")
#     close("all")
# end

##### BINNING #####

binning_plots = False
if binning_plots:
    test_temp = 2.45
    max_binsize = 50

    b = np.array([])
    E_err = np.array([])
    M_err = np.array([])
    Y_err = np.array([])

    print("Importing binning errors data ...")
    with open(f'Analysis/run{run}/binning_errors_{test_temp}.csv', 'r') as csvfile:
        csv_reader = csv.reader(csvfile)
        headers = next(csvfile)

        for row in csv_reader:
            b = np.append(b, float(row[0]))
            E_err = np.append(E_err, float(row[1]))
            M_err = np.append(M_err, float(row[2]))
            Y_err = np.append(Y_err, float(row[3]))
    print("Importing finished ...")

    print("Plotting measurement data ...")
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(b, E_err, color='red', marker='x', label='Energy errors')
    ax.plot(b, M_err, color='blue', marker='x', label='Magnetization errors')
    ax.plot(b, Y_err, color='green', marker='x', label='Spin stiffness errors')
    ax.set_title(f'{EQ_NO} steps + {MEAS_NO} steps, {max_binsize} bins, T = {test_temp}')
    ax.set_ylabel(r'$\Delta_O$')
    ax.set_xlabel('Bin size')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'./Plots/run{run}/binning-errors/EMYvb_{test_temp}.png')
    plt.show()

    print("Plotting finished ...")


