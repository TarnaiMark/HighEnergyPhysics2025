from scipy.optimize import root
import numpy as np
import matplotlib.pyplot as plt

def field_eqs(vars, mu):
    m02 = 2.3925e-4 # GeV^2
    lambda1 = -1.6738
    lambda2 = 23.5078
    c1 = 1.3086 # GeV
    h_N = 5.128e-3 # GeV^3
    h_S = 1.472e-2 # GeV^3
    gF = 4.5708
    M0 = 0.3511 # GeV
    phi_N, phi_S = vars

    # Equation for phi_N
    mN = gF/ 2 * phi_N
    gamma_N = mu/mN

    U_N = m02*phi_N + lambda1*(phi_N**2 + phi_S**2)*phi_N + lambda2/2*phi_N**3 - c1/np.sqrt(2.) * phi_N*phi_S - h_N
    term1_N = 2*gamma_N*(2./3. * gamma_N**2 - 1.) * (gamma_N**2 - 1.)**(1./2.) * (5./2. - gamma_N**2)
    term2_N = 4./3. * gamma_N * (gamma_N**2 - 1.)**(1./2.) - 1
    term3_N = 1./(gamma_N + (gamma_N**2 - 1.)**(-1./2.)) * (gamma_N + (gamma_N**2 - 1.)**(-1./2.) * gamma_N**2)
    log_term_N = -4. * (np.log(mN/M0) + np.log(gamma_N + (gamma_N**2 - 1.)**(1./2.)))
    eq_N = U_N + 3.*gF*mN**3/(8.*np.pi)*(term1_N + term2_N + term3_N + log_term_N)

    # Equation for phi_S
    U_S = m02*phi_S + lambda1*(phi_N**2 + phi_S**2)*phi_S + lambda2/2*phi_S**3 - c1/(2.*np.sqrt(2.)) * phi_N**2 - h_S
    mS = gF/np.sqrt(2.) * phi_S
    gamma_S = mu/mS

    U_S = m02*phi_S + lambda1*(phi_N**2 + phi_S**2)*phi_S + lambda2/2*phi_S**3 - c1/np.sqrt(2.) * phi_N**2 - h_S
    term1_S = 2*gamma_S*(2./3. * gamma_S**2 - 1.) * (gamma_S**2 - 1.)**(1./2.) * (5./2. - gamma_S**2)
    term2_S = 4./3. * gamma_S * (gamma_S**2 - 1.)**(1./2.) - 1
    term3_S = 1./(gamma_S + (gamma_S**2 - 1.)**(-1./2.)) * (gamma_S + (gamma_S**2 - 1.)**(-1./2.) * gamma_S**2)
    log_term_S = -4. * (np.log(mS/M0) + np.log(gamma_S + (gamma_S**2 - 1.)**(1./2.)))
    eq_S = U_S + 3.*gF*mS**3/(np.sqrt(2.)*8.*np.pi)*(term1_S + term2_S + term3_S + log_term_S)

    return [eq_N, eq_S]

# Define values for mu
mu_vals = np.linspace(0,1,50000)

# Solve the system
phi_N_vals,phi_S_vals = [],[]

for mu_val in mu_vals:
    if mu_val == mu_vals[0]:
        initial_guess = [0.1411, 0.1416]

    sol = root(field_eqs, initial_guess, args=(mu_val))

    # Collecting the results
    phi_N_vals.append(sol.x[0])
    phi_S_vals.append(sol.x[1])
    initial_guess = sol.x
    #print(f"mu: {mu_val}, Solution: {sol.x}, Success: {sol.success}, Message: {sol.message}")

# Making plots
plt.figure(figsize=(10, 6))
plt.plot(np.array(mu_vals)*1000, np.array(phi_N_vals)*1000, label="$\\phi_N(\\mu_q)$")
plt.plot(np.array(mu_vals)*1000, np.array(phi_S_vals)*1000, label="$\\phi_S(\\mu_q)$")
plt.title("The $\\mu_f$ dependence of $\\phi_N$ and $\\phi_S$", fontsize=20)
plt.xlabel("$\\mu_W$ [MeV]", fontsize=17)
plt.ylabel("$\\phi_N/S$ [MeV]", fontsize=17)
plt.legend()
plt.savefig("mu_dependence.png")
plt.show()