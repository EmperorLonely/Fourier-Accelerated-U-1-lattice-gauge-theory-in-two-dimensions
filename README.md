# **Fourier Accelerated U(1) Lattice Gauge Theory in Two Dimensions**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

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

## **Installation and Usage**

### **Prerequisites**
- Python 3.x or C++ (depending on your implementation)
- Required libraries: `NumPy`, `SciPy`, `Matplotlib` (for Python implementation)
  
### **Steps to Install**
1. Clone this repository:
    ```bash
    git clone https://github.com/yourusername/U1-lattice-gauge-theory.git
    ```
2. Navigate into the project directory:
    ```bash
    cd U1-lattice-gauge-theory
    ```
3. Install the required dependencies:
    ```bash
    pip install -r requirements.txt
    ```

### **Running the Program**
Run the main script to compute the Wilson loops for a given lattice size and gauge coupling:
```bash
python main.py --lattice_size 10 --gauge_coupling 1.0
