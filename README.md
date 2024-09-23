# **Fourier Accelerated U(1) Lattice Gauge Theory in Two Dimensions**


## **Overview**
This repository contains a **Hybrid Monte Carlo (HMC) program** that computes the expectation value of Wilson loops for discretized U(1) gauge theory in two dimensions. The program utilizes **Fourier acceleration** to speed up convergence and explores its potential to combat **critical slowing down**, a significant issue in lattice gauge theory simulations.

**Objective:**  
Test the effectiveness of Fourier acceleration in U(1) lattice gauge theory simulations to determine whether it provides a substantial advantage over standard algorithms in combating critical slowing down.

### **Key Concepts**
- **Lattice Gauge Theory**: A numerical approach for exploring quantum field theories.
- **Critical Slowing Down**: A phenomenon where the time required for a simulation to converge grows exponentially as lattice spacing approaches zero or as the critical point is approached.
- **Fourier Acceleration**: A technique to accelerate slow modes in the simulation to reduce the impact of critical slowing down.
- **Hybrid Monte Carlo (HMC)**: A Markov Chain Monte Carlo method that uses Hamiltonian dynamics to propose new configurations.

## **Project Features**
- **Hybrid Monte Carlo Algorithm**: Implements the HMC algorithm to simulate U(1) gauge theory.
- **Fourier Acceleration**: Introduces Fourier acceleration to mitigate critical slowing down by modifying the dynamics of the algorithm.
- **Wilson Loops Computation**: The primary observable calculated to measure the effectiveness of the algorithm.

- ## **Background and Theory**

If you're unfamiliar with the concepts of **lattice gauge theory**, **Hybrid Monte Carlo (HMC)**, or **Fourier acceleration**, you can find a detailed explanation in the following project report:

[📄 **Full Report on HMC and Fourier Acceleration**](https://drive.google.com/file/d/1TSl3n5FtMvDfum35ly6HxnjRsgCGJFML/view?usp=sharing)

The report includes:
- **Detailed descriptions** of the HMC algorithm and Fourier acceleration.
- **Pseudocode** for the standard HMC algorithm.
- **Discussion** on the numerical tests performed for U(1) lattice gauge theory.
