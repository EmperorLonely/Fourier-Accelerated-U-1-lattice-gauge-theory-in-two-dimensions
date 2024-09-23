# Fourier accelerated U(1) lattice gauge theory in two dimensions. 

## Hybrid Monte Carlo program that computes the expectation value of Wilson loops for discretized U(1) theory in two dimensions, that utilizes Fourier acceleration to speed up convergence.

The goal of this project is to test the potential of a numerical technique called Fourier acceleration at combating critical slowing down, a long-standing challenge in lattice gauge simulations, and determine whether the accelerated algorithm provides any substantial advantage over the standard one.

In lattice simulations, a framework which enables the numerical exploration of quantum field theories, stochastic methods such as Markov Chain Monte Carlo Methods are the most popular algorithms of choice. When the lattice spacing tends to zero, or a critical point is approached, a phenomenon known as critical slowing down emerges where the time it takes for the simulation to converge diverges to infinity. Fourier acceleration aims to mitigate this by speeding up the slow modes and slowing down the fast modes of the plane wave solution to the theory.

In this project we test Fourier acceleration for U(1) lattice gauge theory in two dimensions using the Hybrid Monte Carlo algorithm (AKA Hamiltonian Monte Carlo.)

The bulk of this repository consists of the Fourier accelerated HMC program used in this project. 

# All necessary background theory, and project details are provided in the following report: 
https://drive.google.com/file/d/1TSl3n5FtMvDfum35ly6HxnjRsgCGJFML/view?usp=sharing

The report provides all necessary background information including a description on the HMC algorithm, and Fourier acceleration, as well as pseudocode for the standard HMC algorithm. Below 


