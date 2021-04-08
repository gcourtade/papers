'''
This script does the following:
1. Reads an enemble of PDB files and experimental PREs
2. Finds the coordinates of the active site histidines
3. Places a copper atom between the His
4. Finds all the backbone amide H atoms
5. Measures the distances between the H amides and the metal
6. Measures the angles between the H amides and the metal
7. Calculate PREs and saves them
8. Plots them together with experimental PREs (H_0p05_PREs.dat)

Part of the script uses code
available from https://github.com/KULL-Centre/DEERpredict

LICENSE: GNU General Public License v3.0

Gaston Courtade - March 2020

'''

import numpy as np
import matplotlib.pyplot as plt
from biopandas.pdb import PandasPdb
import pandas as pd
import math


# user input
# remember to also change the filename variable
# in the run_on_ensemble function
df_ref = pd.read_csv('H_0p05_PREs.dat', sep='\t')
prefix = 'CuBlAA10'
outfile = '{}_calc_PREs_GC.dat'.format(prefix)
other_His = 90 #90 for bla
plot=True
tmix = 10e-7 # mixing time of INEPT transfer (e.g. 10 ms)
R2dia = 12.5 # R2diamagnetic, (e.g 12.5 s-1)
tau_r = 10.2e-9 # protein rotational correlation time (e.g. 10.2 ns)
tau_i = 5e-10 # correlation time of the internal motion (e.g. 500 ps)
B0 = 600 # spectrometer field in MHz (e.g. 600 MHz)
g = 2.13 # electron g-factor (e.g. 2.13)

#read and normalize experimental reference data
ref_res = df_ref['res'].values
ref_err = df_ref['error'].values
ref_data = df_ref['ratio'].values/np.nanmax(df_ref['ratio'].values)

# constant definitions
angstrom_to_meters = 1e-8
mhz_to_rads = 2.0 * math.pi * 1e6
rads_to_mhz = (1/mhz_to_rads)



def run_on_ensemble(nr_of_states, prefix, tmix, R2dia, tau_r, tau_i, B0, g):
    distances_list = []
    H_coordinates = []
    Cu_coordinates = []
    for i in range(nr_of_states):
        filename='5lw4_state{:02d}.pdb_100ps.pdb'.format(i+1)
        ppdb = PandasPdb()
        ppdb.read_pdb(filename)
        #define the coordinates of the metal site
        His1_coord = get_coord(ppdb, 1, 'ND1')
        other_His_coord = get_coord(ppdb, other_His, 'NE2')
        metal_coord = (His1_coord + other_His_coord)/2

        #calcualte distances between the metal site and all amide H
        df = ppdb.df['ATOM']
        distances = ppdb.distance(xyz=metal_coord, records=('ATOM',)).to_frame()
        distances.columns = ['distances']
        df_distances = pd.concat([df,distances], axis=1, join_axes=[df.index])
        amide_H = df_distances.loc[df_distances['atom_name']== 'H']

        residues = np.asarray(amide_H['residue_number'])
        distances = np.asarray(amide_H['distances'])*angstrom_to_meters #in meters
        distances_list.append(distances)
        amide_H_coord = get_coord(ppdb, 'all', 'H')
        amide_H_coord = np.column_stack((amide_H_coord[0],amide_H_coord[1],amide_H_coord[2]))
        H_coordinates.append(amide_H_coord)
        Cu_coordinates.append(metal_coord)

    H_coordinates = np.average(np.array(H_coordinates), axis=0)
    Cu_coordinates = np.average(np.array(Cu_coordinates), axis=0)
    h_probe_vector = H_coordinates - Cu_coordinates
    vect = np.sum(h_probe_vector*angstrom_to_meters)
    cos = vect/np.average(distances_list, axis=0)
    r6_dist = np.average(distances_list, axis=0)**(-6)
    r3_dist = np.average(distances_list, axis=0)**(-3)    
    gamma_2_PREs = gamma_2(r6_dist, r3_dist,cos, tau_r, 10e-7, tau_i, 0.5, B0, g)
    I_ratio = calc_I_ratio(gamma_2_PREs, tmix, R2dia)
    return residues,I_ratio


def get_coord(ppdb, residue_number, atom_name):
    coord = []
    if residue_number=='all':
        coord.append(ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name'] == atom_name]['x_coord'])
        coord.append(ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name'] == atom_name]['y_coord'])
        coord.append(ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name'] == atom_name]['z_coord'])
    else:
        
        coord.append(float(ppdb.df['ATOM'][ppdb.df['ATOM']['residue_number'] == residue_number][ppdb.df['ATOM']['atom_name'] == atom_name]['x_coord']))
        coord.append(float(ppdb.df['ATOM'][ppdb.df['ATOM']['residue_number'] == residue_number][ppdb.df['ATOM']['atom_name'] == atom_name]['y_coord']))
        coord.append(float(ppdb.df['ATOM'][ppdb.df['ATOM']['residue_number'] == residue_number][ppdb.df['ATOM']['atom_name'] == atom_name]['z_coord']))
    coord = np.asarray(coord)
    return coord

def chi2(ref,calc):
    return np.nansum( (ref-calc)**2 / ref )


def J_SBMF(r6_dist, tau_r, tau_s, tau_i, wh, S2):
    #tau_c = 1/((1/tau_r) + (1/tau_s))
    tau_c = tau_r
 #   tau_t = 1/((1/tau_r) + (1/tau_s) + (1/tau_i))
    tau_t = 1/((1/tau_r) + (1/tau_i))
    return r6_dist * ( (S2*tau_c/(1+(wh*tau_c)**2)) + ((1-S2)*tau_t/(1+(wh*tau_c)**2)) )

def gamma_2(r6_dist, r3_dist, cos, tau_r, tau_s, tau_i, s, wh_mhz, g):
    mu0 = 1.257e-6 # permeability of free space (m kg s-2 Å-2)
    muB = 9.274e-24 # Bohr magneton (J T-1)
    lambda_h = 2.675e8 # proton gyromagnetic ratio (rad s-1 T-1)
    wh = mhz_to_rads * wh_mhz
    S2_radial = (r6_dist**-1)*(r3_dist**2)

    S2_angular = (((3 / 2) * np.power(cos, 2)) - 0.5)
    
    S2 = S2_angular * S2_radial

    return (1/15)*(mu0/(4*math.pi))*(lambda_h**2)*(g**2)*(muB**2)*s*(s+1)*(4*J_SBMF(r6_dist, tau_r, tau_s, tau_i, 0, S2) + 3*J_SBMF(r6_dist, tau_r, tau_s, tau_i, wh, S2))

def calc_I_ratio(gamma_2_PREs, mixing_time, R2dia):
    return (R2dia/(R2dia + gamma_2_PREs))*np.exp(-1*gamma_2_PREs*mixing_time)



# calculate PREs and save them
aashift = 31 #first residue in BlLPMO10A is 32
calc_res, calc_data = run_on_ensemble(20, prefix, tmix, R2dia, tau_r, tau_i, B0, g)
with open(outfile, 'w') as f:
    #f.write('res\texp\terror\tcalc\n')
    #file written as: res exp error calc
    for i in range(1,173):
        c = np.nan
        r = np.nan
        e = np.nan
        if i in calc_res:
            idx_calc = list(calc_res).index(i)
            c = calc_data[idx_calc]
            if i in ref_res:
                idx_ref = list(ref_res).index(i)
                r = ref_data[idx_ref]
                e = ref_err[idx_ref]
        f.write(f'{i+aashift}\t{r}\t{e}\t{c}\n')

#calculate reduced chi2
chi2 = chi2(ref_data, calc_data[1:])
print('chi2r', chi2)


#plot results
if plot:
    fsize=16
    data = np.loadtxt(outfile)

    plt.figure(figsize=(7.08,2.75), dpi=300)
    plt.plot(data[:,0], data[:,1], color='k', linewidth=1, label='Exp')
    plt.plot(data[:,0], data[:,3], color='r', linewidth=1, label='Calc')
    plt.fill_between(data[:,0], data[:,1]-data[:,2],data[:,1]+data[:,2], color='grey', alpha=0.6)
    plt.xlabel('Amino acids', fontsize=fsize)
    plt.xlim(1+aashift+1,172+aashift+1)
    plt.ylabel('I$_{Cu}$/I$_{apo}$', fontsize=fsize)
    plt.tight_layout()
    plt.savefig('PRE_results_v2.png')
    plt.show()




