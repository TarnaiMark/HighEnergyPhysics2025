import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Data loader function (using Pandas library)
def pandasLoadandRaa(file, skiprows=28):
    print("Processing: ", file)
    data_all = pd.read_csv(file,skiprows=skiprows,delim_whitespace=True,header=None,names=['a','id','0','1','px','py','pz','e'])
    data = data_all.dropna()
    N_evt = (len(data_all) - len(data))/2.0
    return data, N_evt

# File names of Au+Au and p+p collisions from different centrality classes
pp_fnames = [
    "/project/femtoscopy/data/raa_hijing/pp_200_200000_ev_0",
    "/project/femtoscopy/data/raa_hijing/pp_200_200000_ev_1",
    "/project/femtoscopy/data/raa_hijing/pp_200_200000_ev_2",
    "/project/femtoscopy/data/raa_hijing/pp_200_200000_ev_3",
    "/project/femtoscopy/data/raa_hijing/pp_200_200000_ev_4",
    "/project/femtoscopy/data/raa_hijing/pp_200_200000_ev_5",
    
]
# "/project/femtoscopy/data/raa_hijing/pp_200_200000_ev_6",
#    "/project/femtoscopy/data/raa_hijing/pp_200_200000_ev_7"

auau_fnames = [
    "/project/femtoscopy/data/raa_hijing/AuAu_200_0_5_cent_1000_ev",
    "/project/femtoscopy/data/raa_hijing/AuAu_200_5_10_cent_1000_ev",
    "/project/femtoscopy/data/raa_hijing/AuAu_200_10_20_cent_1000_ev",
    "/project/femtoscopy/data/raa_hijing/AuAu_200_20_30_cent_1000_ev",
    "/project/femtoscopy/data/raa_hijing/AuAu_200_30_40_cent_1500_ev",
    "/project/femtoscopy/data/raa_hijing/AuAu_200_40_50_cent_1000_ev"
]

n_coll = [1065.4, 845.4, 602.6, 373.8, 219.8, 120.3]
n_part = [351.4, 299.0, 234.6, 166.6, 114.2, 74.4]

# Initialize lists for RAA values
RAA_ph = []
yerr_ph = []

RAA_pipos = []
yerr_pipos = []

RAA_pizero = []
yerr_pizero = []

RAA_piminus = []
yerr_piminus = []


# Looping through all centrality classes
for i in range(len(auau_fnames)):
    # Loading the data
    data_pp, N_evt_pp = pandasLoadandRaa(pp_fnames[i])
    data_auau, N_evt_au = pandasLoadandRaa(auau_fnames[i])

    # Calculate the number of photons (and other particles) from the initial Au+Au and p+p collisions
    ph_numau = len(data_auau[(data_auau['id'] == 22) & (data_auau['0'] == 0)])
    ph_numpp = len(data_pp[(data_pp['id'] == 22) & (data_pp['0'] == 0)])
    
    # Calcuaate the number of the three different pions (does not need to be direct)
    pipos_numau = len(data_auau[(data_auau['id'] == 211)])
    pizero_numau = len(data_auau[(data_auau['id'] == 111)])
    piminus_numau = len(data_auau[(data_auau['id'] == -211)])
    
    pipos_numpp = len(data_pp[(data_pp['id'] == 211)])
    pizero_numpp = len(data_pp[(data_pp['id'] == 111)])
    piminus_numpp = len(data_pp[(data_pp['id'] == -211)])
    
    # Calculate RAA and its error with basic error propagation
    Raa_ph = (ph_numau / N_evt_au) / (ph_numpp / N_evt_pp * n_coll[i])
    RAA_ph.append(Raa_ph)
    yerr_ph.append(Raa_ph * np.sqrt((np.sqrt(ph_numau)/ph_numau)**2 + (np.sqrt(ph_numpp)/ph_numpp)**2))
    
    # Calculate RAA for pions
    Raa_pipos = (pipos_numau / N_evt_au) / (pipos_numpp / N_evt_pp * n_coll[i])
    RAA_pipos.append(Raa_pipos)
    yerr_pipos.append(Raa_pipos * np.sqrt((np.sqrt(pipos_numau)/pipos_numau)**2 + (np.sqrt(pipos_numpp)/pipos_numpp)**2))
    
    Raa_pizero = (pizero_numau / N_evt_au) / (pizero_numpp / N_evt_pp * n_coll[i])
    RAA_pizero.append(Raa_pizero)
    yerr_pizero.append(Raa_pizero * np.sqrt((np.sqrt(pizero_numau)/pizero_numau)**2 + (np.sqrt(pizero_numpp)/pizero_numpp)**2))
    
    Raa_piminus = (piminus_numau / N_evt_au) / (piminus_numpp / N_evt_pp * n_coll[i])
    RAA_piminus.append(Raa_piminus)
    yerr_piminus.append(Raa_piminus * np.sqrt((np.sqrt(piminus_numau)/piminus_numau)**2 + (np.sqrt(piminus_numpp)/piminus_numpp)**2))

# Plotting the results
plt.plot([0,360],[1,1],'k')
plt.errorbar(n_part, RAA_ph, yerr=yerr_ph, capsize=4.5, marker = 'x', label=r'$\gamma$', linestyle='', fillstyle='none')
plt.errorbar(n_part, RAA_pipos, yerr=yerr_pipos, capsize=4.5, marker = 'x', label=r'$\pi^+$', linestyle='', fillstyle='none')
plt.errorbar(n_part, RAA_pizero, yerr=yerr_pizero, capsize=4.5, marker = 'x', label=r'$\pi^0$', linestyle='', fillstyle='none')
plt.errorbar(n_part, RAA_piminus, yerr=yerr_piminus, capsize=4.5, marker = 'x', label=r'$\pi^-$', linestyle='', fillstyle='none')


plt.xlim(0.0,360.0)
plt.grid(True)
plt.xlabel(r'N$_{part}$')
plt.ylabel(r'R$_{AA}$')
plt.legend()
plt.savefig("R_aa.png")

