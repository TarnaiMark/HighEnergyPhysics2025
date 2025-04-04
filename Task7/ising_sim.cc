// 2D Ising model — compute cumulant ratios C4/C2 and X4/X2 - 3*X2
// TASK: Complete the missing parts (marked with TODO) to implement the simulation

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <iomanip>

class IsingModel {
public:
    IsingModel(int L, double T);
    void run(int steps);                            // Run Monte Carlo steps (without measurement)
    void computeRatios(double &C42, double &R42);   // Measure and compute ratios

private:
    int L;                                          // Lattice linear size (L x L)
    double T;                                       // Temperature
    std::vector<std::vector<int>> lattice;          // 2D lattice of spins
    std::mt19937 rng;                               // Random number generator
    std::uniform_real_distribution<double> dist;

    void initialize();                              // Initialize the lattice
    void metropolisStep();                          // Perform one Metropolis step
    int totalMagnetization();                       // Compute total magnetization
};

// Constructor: set lattice size and temperature
IsingModel::IsingModel(int L, double T) : L(L), T(T), dist(0.0, 1.0) {
    rng.seed(std::random_device{}());
    initialize();
}

// Initialize the lattice: all spins up (+1)
void IsingModel::initialize() {
    lattice.resize(L, std::vector<int>(L, 1)); // 
}

// Perform a full Metropolis update sweep (L*L random spin flips)
void IsingModel::metropolisStep() {
    std::uniform_int_distribution<int> siteDist(0, L - 1);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i = 0; i < L * L; ++i) {
        int x = siteDist(rng);
        int y = siteDist(rng);

        // TODO: Compute the energy change ΔE for flipping spin at (x, y)
        // Use periodic boundary conditions: (x+1)%L, (x-1+L)%L, etc.
        // ΔE = 2 * s_xy * sum of neighbors
        
        int s = lattice[x][y];
        int neighborSum = lattice[x][(y + 1) % L] + lattice[x][(y - 1 + L) % L] + lattice[(x - 1 + L) % L][y] + lattice[(x + 1) % L][y];
        double d_E = 2 * s * neighborSum;

        // TODO: Accept the flip with probability min(1, exp(-ΔE / T))
        // If accepted, flip the spin: lattice[x][y] *= -1;
        if (d_E < 0 || dist(rng) < exp(-d_E / T)) {
            lattice[x][y] *= -1; // Flip the spin
        }
    }
}

// Compute total magnetization: sum of all spins
int IsingModel::totalMagnetization() {
    int M = 0;
    
    // TODO: Sum all values in the lattice and return total magnetization
    for(int i = 0; i < L ; ++i){
        for(int j = 0; j < L ; ++j){
            M += lattice[i][j];
        }
    }

    return M;
}

// Measure magnetization statistics and compute cumulant ratios
void IsingModel::computeRatios(double &C42, double &R42) {
    
    int samples = 100000;
    std::vector<double> mvals;
    double M1_sum = 0.0;
    for (int i = 0; i < samples; ++i) {
        metropolisStep();
        int M = totalMagnetization();
        double m = static_cast<double>(M) / (L * L); // rescaling may not be needed for this type of C_{n}
        mvals.push_back(m);
        M1_sum += m;
    }
    double M1 = M1_sum / samples;
    double  M2 = 0.0, M3 = 0.0, M4 = 0.0;
    
    // TODO: Compute raw moments M1–M4
    for (int i = 0; i < samples; ++i) {
        M2 += mvals[i] * mvals[i];
        M3 += mvals[i] * mvals[i] * mvals[i];
        M4 += mvals[i] * mvals[i] * mvals[i] * mvals[i];
    }
    M2 /= samples;
    M3 /= samples;
    M4 /= samples;

    // TODO: Compute cumulants C2 and C4
    
    double C1 = M1;
    double C2 = M2 - M1 * M1;
    double C4 = M4 - 4 * M3 * M1 + 6 * M2 * M1 * M1 - 3 * M1 * M1 * M1 * M1;

    //double R42 = ...


    C42 = C4/C2;
    R42 = C42 - 3 * C2;
}

// Run given number of sweeps (without measurement)
void IsingModel::run(int steps) {
    for (int i = 0; i < steps; ++i)
        metropolisStep();
}

int main() {
    int L = 30;
    std::cout << "# T\tC4/C2\tX4/X2 - 3*X2" << std::endl;

    double Tmin = 0.1;
    double Tmax = 5.0;
    double dTdefault = 0.1;

    // Estimate total number of steps for progress bar
    int totalSteps = 0;
    for (double T = Tmin; T < Tmax; ) {
        double dT = (T > 2 && T < 2.6) ? dTdefault / 2 : dTdefault;
        T += dT;
        ++totalSteps;
    }

    
    int currentStep = 0;
    double dT = dTdefault;

    // Main simulation loop
    for (double T = Tmin; T < Tmax; ) {
        dT = (T > 2 && T < 2.6) ? dTdefault / 2 : dTdefault;

        //Using the written function
        IsingModel model(L, T);
        model.run(100000);

        double C42 = 0, R42 = 0;
        model.computeRatios(C42, R42);

        std::cout << T << "\t" << C42 << "\t" << R42 << std::endl;

        // --- PROGRESS BAR ---
        ++currentStep;
        double percent = 100.0 * currentStep / totalSteps;
        int barWidth = 40;
        std::cerr << "\rProgress: [";
        int pos = static_cast<int>(barWidth * percent / 100.0);
        for (int i = 0; i < barWidth; ++i)
            std::cerr << (i < pos ? "=" : " ");
        std::cerr << "] " << std::fixed << std::setprecision(1) << percent << "%";
        std::cerr.flush();

        T += dT;
    }

    std::cerr << std::endl; // Final newline after progress bar
    return 0;
}