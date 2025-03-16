import time
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool

from ising_simulation import simulate
from syk_simulation import (
    calculate_entropy,
    calculate_otoc,
    calculate_specific_heat,
    run_multiple_syk_simulations,
)

st.set_page_config(
    layout="wide",
    initial_sidebar_state="auto",
    page_title="Quantum Models Simulation",
    page_icon="🧊",
)

st.title("Quantum Models Simulation")
st.sidebar.title("Simulation Options")

# Add select box for choosing between Ising and SYK models
model = st.sidebar.selectbox(
    "Choose the model to simulate:", ("3D Ising Model", "SYK Model")
)

if model == "3D Ising Model":
    N = st.sidebar.number_input("Lattice size (N)", min_value=5, max_value=50, value=20)
    J = st.sidebar.number_input(
        "Interaction strength (J)", min_value=0.1, max_value=10.0, value=1.0
    )
    h = st.sidebar.number_input(
        "External magnetic field (h)", min_value=0.0, max_value=10.0, value=0.1
    )
    T_min = st.sidebar.number_input(
        "Minimum temperature (T_min)", min_value=0.1, max_value=10.0, value=1.0
    )
    T_max = st.sidebar.number_input(
        "Maximum temperature (T_max)", min_value=T_min, max_value=10.0, value=4.0
    )
    T_step = st.sidebar.number_input(
        "Temperature step (T_step)", min_value=0.1, max_value=1.0, value=0.1
    )
    steps = st.sidebar.number_input(
        "Number of Monte Carlo steps",
        min_value=1000,
        max_value=1000000,
        value=100000,
        step=1000,
    )
    equilibration_steps = st.sidebar.number_input(
        "Number of equilibration steps",
        min_value=0,
        max_value=100000,
        value=10000,
        step=1000,
    )

    if st.button("Run simulation"):
        temperatures = np.arange(T_min, T_max + T_step, T_step)
        magnetizations = []
        energies = []
        susceptibilities = []
        heat_capacities = []

        simulation_args = [
            (N, J, h, T, steps, equilibration_steps) for T in temperatures
        ]

        # Use parallelization to speed up simulations
        start = time.time()
        progress_bar = st.progress(0, text="Simulating...")
        results = []
        with Pool() as p:
            for i, result in enumerate(p.imap_unordered(simulate, simulation_args)):
                results.append(result)
                progress_bar.progress((i + 1) / len(simulation_args))
        end = time.time()
        st.write("Simulation time: {:.2f} seconds".format(end - start))

        for T, (magnetization, energy, susceptibility, heat_capacity) in zip(
            temperatures, results
        ):
            magnetizations.append(magnetization)
            energies.append(energy)
            susceptibilities.append(susceptibility)
            heat_capacities.append(heat_capacity)

        [col1, col2] = st.columns(2)
        [col3, col4] = st.columns(2)

        # Magnetization plot
        fig, ax = plt.subplots()
        ax.plot(
            temperatures,
            magnetizations,
            marker="o",
            linestyle="-",
            label="Magnetization",
        )
        ax.set_xlabel("Temperature")
        ax.set_ylabel("Magnetization")
        ax.set_title("Magnetization vs. Temperature")
        ax.legend()
        ax.grid()
        with col1:
            st.pyplot(fig)

        # Energy plot
        fig, ax = plt.subplots()
        ax.plot(temperatures, energies, marker="o", linestyle="-", label="Energy")
        ax.set_xlabel("Temperature")
        ax.set_ylabel("Energy")
        ax.set_title("Energy vs. Temperature")
        ax.legend()
        ax.grid()
        with col2:
            st.pyplot(fig)

        # Susceptibility plot
        fig, ax = plt.subplots()
        ax.plot(
            temperatures,
            susceptibilities,
            marker="o",
            linestyle="-",
            label="Susceptibility",
        )
        ax.set_xlabel("Temperature")
        ax.set_ylabel("Susceptibility")
        ax.set_title("Susceptibility vs. Temperature")
        ax.legend()
        ax.grid()
        with col3:
            st.pyplot(fig)

        # Heat capacity plot
        fig, ax = plt.subplots()
        ax.plot(
            temperatures,
            heat_capacities,
            marker="o",
            linestyle="-",
            label="Heat Capacity",
        )
        ax.set_xlabel("Temperature")
        ax.set_ylabel("Heat Capacity")
        ax.set_title("Heat Capacity vs. Temperature")
        ax.legend()
        ax.grid()
        with col4:
            st.pyplot(fig)


elif model == "SYK Model":
    N = st.sidebar.slider("Number of Majorana fermions (N):", 4, 10, 4, step=2)
    T_values = np.linspace(0.1, 10, 100)
    num_samples = st.sidebar.slider("Number of simulations for averaging:", 1, 100, 10)
    beta = st.sidebar.slider(
        "Inverse temperature (beta)", min_value=0.1, max_value=10.0, value=1.0, step=0.1
    )
    J = st.sidebar.slider(
        "Coupling strength (J)", min_value=0.1, max_value=10.0, value=1.0, step=0.1
    )
    time_steps = st.sidebar.slider(
        "Number of time steps", min_value=1, max_value=1000, value=10, step=1
    )

    if st.button("Plot scrambling effect"):
        times, oto_correlations = calculate_otoc(N, J, num_samples, time_steps)
        plt.figure(figsize=(10, 6))
        plt.plot(times, oto_correlations, marker="o")
        plt.xlabel("Time")
        plt.ylabel("OTOC")
        plt.title("Scrambling effect in SYK model")
        plt.grid()
        st.pyplot(plt)

    if st.button("Run SYK Model Simulation"):
        eigenvalues_list, eigenvectors_list = run_multiple_syk_simulations(
            N, num_samples
        )

        entropies_list = [
            calculate_entropy(eigenvalues, T_values) for eigenvalues in eigenvalues_list
        ]
        specific_heats_list = [
            calculate_specific_heat(eigenvalues, T_values)
            for eigenvalues in eigenvalues_list
        ]

        average_entropy = np.mean(entropies_list, axis=0)
        average_specific_heat = np.mean(specific_heats_list, axis=0)
        average_eigenvalues = np.mean(eigenvalues_list, axis=0)
        average_eigenvectors = np.mean(eigenvectors_list, axis=0)

        fig, ax = plt.subplots()
        ax.plot(average_eigenvalues, "o")
        ax.set(
            xlabel="Eigenvalue index",
            ylabel="Energy",
            title="Energy levels of the SYK model",
        )
        st.pyplot(fig)

        st.write("Entropy vs. Temperature:")
        fig, ax = plt.subplots()
        ax.plot(T_values, average_entropy)
        ax.set(xlabel="Temperature", ylabel="Entropy", title="Entropy vs. Temperature")
        plt.close(fig)
        st.pyplot(fig)

        st.write("Specific Heat vs. Temperature:")
        fig, ax = plt.subplots()
        ax.plot(T_values, average_specific_heat)
        ax.set(
            xlabel="Temperature",
            ylabel="Specific Heat",
            title="Specific Heat vs. Temperature",
        )
        plt.close(fig)
        st.pyplot(fig)
