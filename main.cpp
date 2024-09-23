#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <numeric>
#include <chrono>
#include <fftw3.h>
#include <complex>
#include <fstream>
#include <iomanip>
#include <string>

using Complex = std::complex<double>;

const int L = 200; //Width of the lattice
const int N = L*L; //Total number of lattice sites
double step_size = 0.027; // Leapfrog integration time step size
const int num_steps = 20; // Nr. of leapfrog steps
const int num_samples = 50000; // Total nr. of Monte Carlo trajectories / samples
const int therm = 5000; // Nr. of samples thermalized
double beta = 19; //Inverse temperature, controls interaction strength
int length = 8; //Wilson loop length
int width = 8; //Wilson loop width
double M = 0.9; //Acceleration parameter 

//
std::vector<double> plaq(N, 0);  
std::vector<double> lat(2*N, 0);
std::vector<double> q(2*N, 0);
std::vector<double> p(2*N, 0);
std::vector<double> force(2*N, 0);
std::vector<double> samples;
std::vector<double> dH;
std::vector<double> A(N, 0);


//Used later for generating random values
std::mt19937 gen(12345);
std::uniform_real_distribution<> dis(-M_PI, M_PI);
std::uniform_real_distribution<> metropolis(0.0, 1.0);
std::normal_distribution<> normal_dist(0.0, 1.0);

// Declare FFTW plans globally to reuse them
fftw_plan forward_plan_2d, backward_plan_2d;

// Declare global reusable complex vectors for FFTs
std::vector<Complex> p_0(N, 0), p_1(N, 0), q_0(N, 0), q_1(N, 0);

// Initialize FFTW plans for efficient reuse.
// FFTW plans are computationally expensive to create, so they are initialized once
// and reused throughout the simulation to optimize performance.
void init_fftw_plans() {
    forward_plan_2d = fftw_plan_dft_2d(L, L, reinterpret_cast<fftw_complex*>(p_0.data()), reinterpret_cast<fftw_complex*>(p_0.data()), FFTW_FORWARD, FFTW_ESTIMATE);
    backward_plan_2d = fftw_plan_dft_2d(L, L, reinterpret_cast<fftw_complex*>(p_0.data()), reinterpret_cast<fftw_complex*>(p_0.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
}

void destroy_fftw_plans() {
    fftw_destroy_plan(forward_plan_2d);
    fftw_destroy_plan(backward_plan_2d);
    fftw_cleanup();
}

// Constructs the mass matrix for Fourier acceleration in momentum space.
// Formula accounts for Fourier modes in both x and y directions.
void construct_a(std::vector<double> &A) {
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            // Fourier acceleration factor in x and y directions 
            A[x * L + y] = 1.0 / ((1.0 - M) + M * 4.0 * (sin(M_PI * double(x) / double(L)) * sin(M_PI * double(x) / double(L)) + sin(M_PI * double(y) / double(L)) * sin(M_PI * double(y) / double(L))));
        }
    }
}

// Ensures Hermitian symmetry in the momentum-space lattice.
// This is required to guarantee real values when performing the inverse FFT.
void make_hermitian(std::vector<Complex> &p) {
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            if (x == 0 || x == L / 2) {
                if (y == 0 || y == L / 2) {
                    p[x * L + y] = {p[x * L + y].real(), 0};
                } else {
                    // Apply symmetry conditions across the lattice
                    p[x * L + y] = {p[x * L + y].real(), p[L * x + (L - y)].imag() };
                    p[y * L + x] = {p[(y) * L + x].real() , p[(L - y) * L + x].imag() };
                    p[x * L + (L - y)] = std::conj(p[x * L + y]);
                    p[(L - y) * L + x] = std::conj(p[y * L + x]);
                }
            } else if (y != 0 && y != L / 2) {
                p[x * L + y] = {p[x * L + y].real() , p[(L - x) * L + (L - y)].imag() };
                p[(L - x) * L + (L - y)] = std::conj(p[x * L + y]);
            }
        }
    }
}

//Performs Fast Fourier Transform on vector p
void forward_fft(std::vector<Complex> &p) {
    fftw_execute_dft(forward_plan_2d, reinterpret_cast<fftw_complex*>(p.data()), reinterpret_cast<fftw_complex*>(p.data()));
}

//Performs inverse Fast Fourier Transform on vector p
void backward_fft(std::vector<Complex> &p) {
    fftw_execute_dft(backward_plan_2d, reinterpret_cast<fftw_complex*>(p.data()), reinterpret_cast<fftw_complex*>(p.data()));
    for (auto &i : p) {
        i /= L * L;
    }
}

//Generate initial fictitious momentum distribution, starting in momentum space
void sample_p(std::vector<double> &p, std::vector<double> &A) {
    for (int i = 0; i < N; i++) {
        double variance = sqrt(N / A[i]);
        std::normal_distribution<> normal_dist(0, variance);
        p_0[i] = {normal_dist(gen)/sqrt(2), normal_dist(gen)/sqrt(2)};
        p_1[i] = {normal_dist(gen)/sqrt(2), normal_dist(gen)/sqrt(2)};
    }

    make_hermitian(p_0);
    make_hermitian(p_1);

    //Transform to position space
    backward_fft(p_0);
    backward_fft(p_1);

    for (int i = 0; i < N; i++) {
        p[i] = p_0[i].real();
        p[i + N] = p_1[i].real();
    }
}

//Perform leapfrog update on the gauge field links
void update_q(std::vector<double> &q, std::vector<double> &p, const std::vector<double> &A, double step_size) {
    for (int i = 0; i < N; i++) {
        p_0[i] = p[i];
        p_1[i] = p[i + N];
        q_0[i] = q[i];
        q_1[i] = q[i+N];
    }

    forward_fft(p_0);
    forward_fft(p_1);
    forward_fft(q_0);
    forward_fft(q_1);

    for (int i = 0; i < N; i++) {
        q_0[i] += step_size*A[i]*p_0[i];
        q_1[i] += step_size*A[i]*p_1[i];
    }
    backward_fft(q_0);
    backward_fft(q_1);

    for (int i = 0; i < N; i++) {
        q[i] = q_0[i].real();
        q[i+N] = q_1[i].real();
    }
}

//Compute kinetic energy 
double compute_ke(std::vector<double> &p, std::vector<double> &A) {
    for (int i = 0; i < N; i++) {
        p_0[i] = p[i];
        p_1[i] = p[i+N];
    }

    forward_fft(p_0);
    forward_fft(p_1);

    double sum = 0;
    for (int i = 0; i < N; i++) {
        sum += A[i]*(std::norm(p_0[i]));
        sum += A[i]*(std::norm(p_1[i]));
    }
    return 0.5*sum/(N);
}

//Completely randomize a lattice
void hot_start(std::vector<double> &lat) {
    for (int i = 0; i < 2*N; i++) {
        lat[i] = dis(gen);
    }
}

//Calculate the link angle variable for all lattice sites / plaquettes and store into a vector
void construct_plaquette(std::vector<double> &lat, std::vector<double> &P) {
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            P[x*L + y] = lat[x*L + y] + lat[((x+1)%L)*L + y + N] - lat[x*L + (y+1)%L] - lat[N + x*L + y];
        }
    }
}

//Compute Wilson action of lattice with plaquettes P
double compute_action(const std::vector<double> &P) {
    double S = 0.0;
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            S += -beta * cos(P[x*L + y]);
        }
    }
    return S;
}

//Compute Hamiltonian
double compute_h(const std::vector<double> &plaq, const std::vector<double> &p) {
    double H = compute_action(plaq);
    for (int i = 0; i < 2*N; i++) {
        H += 0.5 * p[i] * p[i];
    }
    return H;
}

//Compute gradient of action
void compute_gradient(std::vector<double> &force, const std::vector<double> &P) {
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            force[x*L + y] = beta * (sin(P[x*L + y]) - sin(P[x*L + (y - 1 + L) % L]));
            force[N + x*L + y] = beta * (-sin(P[x*L + y]) + sin(P[((x - 1 + L) % L)*L + y]));
        }
    }
}

//Compute the expectation value of Wilson loops over the lattice
//Length and width of Wilson loops defined at the top
double wilson_loop(const std::vector<double> &lat) {
    double W = 0.0;
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            double S = 0.0;
            for (int i = 0; i < width; i++) {
                S += lat[((x + i) % L)*L + y];
                S -= lat[((x + i) % L)*L + (y + length) % L];
            }
            for (int i = 0; i < length; i++) {
                S -= lat[N + x*L + (y + i) % L];
                S += lat[N + ((x + width) % L)*L + (y + i) % L];
            }
            W += cos(S);
        }
    }
    return W / (L * L);
}

// Leapfrog integration algorithm for the Hybrid Monte Carlo update.
// The algorithm alternates between updating the momentum (p) and the gauge field (q).
// The force is computed from the gradient of the action, and the momentum is updated in half-steps.
void leap_frog(std::vector<double> &q, std::vector<double> &p, std::vector<double> A) {
    construct_plaquette(q, plaq);
    compute_gradient(force, plaq);

    // Initial half-step momentum update
    for (int i = 0; i < 2*N; i++) {
        p[i] -= 0.5 * step_size * force[i];
    }

    // Perform num_steps leapfrog steps
    for (int step = 0; step < num_steps - 1; step++) {
        update_q(q, p, A, step_size);

        construct_plaquette(q, plaq);
        compute_gradient(force, plaq);

        // Full-step momentum update
        for (int i = 0; i < 2*N; i++) {
            p[i] -= step_size * force[i];
        }
    }

     // Final half-step momentum update
    update_q(q, p, A, step_size);

    construct_plaquette(q, plaq);
    compute_gradient(force, plaq);

    for (int i = 0; i < 2*N; i++) {
        p[i] -= 0.5 * step_size * force[i];
    }
}

//
bool hmc_update(std::vector<double> &lat) {
    q = lat;
    construct_plaquette(q, plaq);
    sample_p(p, A);

    double H_old = compute_action(plaq) + compute_ke(p, A);
    leap_frog(q, p, A);
    construct_plaquette(q, plaq);
    double H_new = compute_action(plaq) + compute_ke(p, A);

    double d_H = H_old - H_new;
    dH.push_back(exp(d_H));

    if (metropolis(gen) < exp(d_H)) {
        lat = q;
        return true;
    } else {
        return false;
    }
}

int main() {
    // Initialize FFTW plans and global vectors for FFTs
    init_fftw_plans();
    
    // Filename to store results, based on lattice size, beta, length, and M
    std::string file_name = std::to_string(L) + "_" + std::to_string(int(beta)) + "_" + std::to_string(length) + "_" + std::to_string(M) + ".txt";

    // Construct the mass matrix for Fourier acceleration
    construct_a(A);

    // Measure the time taken for the simulation
    auto start = std::chrono::high_resolution_clock::now();

    // Run Monte Carlo updates until acceptance rate is within the desired range
    // More sophisticated methods to convergence within accepted acceptance range
    // should be considered.
    while (true) {
        double accept = 0;
        hot_start(lat);// Randomize the lattice
        samples.clear(); // Clear samples instead of redeclaring

        // Loop over the number of Monte Carlo samples
        for (int i = 0; i < num_samples; i++) {
            // Check acceptance rate and adjust step size if needed
            if (i % 1000 == 0) {
                if (accept / i < 0.60 && beta != 0) {
                    accept = (accept*num_samples)/i;
                    break;
                } else if (accept / i > 0.85 && beta != 0) {
                    accept = (accept*num_samples)/i;
                    break;
                }
            }

            // Perform HMC update and collect Wilson loop sample
            if (hmc_update(lat) == true) {
                accept += 1;
            }
            samples.push_back(wilson_loop(lat));
        }

        // Adjust step size based on acceptance rate
        if (accept / num_samples > 0.60 && accept / num_samples < 0.85 || beta == 0) {
            std::cout << "Step_size & Acceptance rate: " <<  step_size << " " << accept / num_samples  << std::endl;
            break;
        } else if (accept / num_samples < 0.60) {
            std::cout << "Too large " << accept / num_samples << std::endl;
            step_size *= 0.98;
        } else {
            std::cout << "Too small " << accept / num_samples << std::endl;
            step_size *= 1.03;
        }
    }

    // Measure elapsed time and print
    auto end = std::chrono::high_resolution_clock::now();  // End timing
    std::chrono::duration<double> elapsed = end - start;  // Calculate elapsed time
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;  // Print elapsed time

    // Output results to a file
    std::ofstream outfile(file_name);
    outfile << std::fixed << std::setprecision(16);

    for (double i : samples) {
        outfile << i << std::endl;
    }

    outfile.close();



    // Clean up FFTW plans before exiting
    destroy_fftw_plans();
    
    return 0;
}
