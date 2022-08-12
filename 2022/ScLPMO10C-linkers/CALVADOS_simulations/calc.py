from analyse import *
import MDAnalysis
import time
import os
import glob
import sys
import psutil
import logging
import shutil
import h5py
import pandas as pd
import numpy as np
import mdtraj as md
from mdtraj import element
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
import time

proteins = pd.read_pickle('proteins.pkl')

def calcRs(traj):
    pairs = traj.top.select_pairs('all','all')
    d = md.compute_distances(traj,pairs)
    dmean = d.mean(axis=0)
    ij = np.array(range(1,traj.n_atoms))
    diff = [x[1]-x[0] for x in pairs]
    dij = np.empty(0)
    for i in ij:
        dij = np.append(dij,dmean[diff==i].mean())
    return ij, dij

def calcCs(t):
    bonds = t.xyz[:,1:,:] - t.xyz[:,:-1,:]
    l2 = np.mean(np.linalg.norm(bonds, axis=2))**2 #number of frames, beads, xyz (Cartesian coords)
    dot = np.einsum('ijk,ilk->ilj',bonds,bonds)
    dmean = dot.mean(axis=0)
    N = t.n_atoms
    diff = np.abs(np.arange(N-1)-np.arange(N-1)[:,np.newaxis])
    dij = np.empty(0)
    for j in range(1,N-2):
        dij = np.append(dij,dmean[diff==j].mean())
    return np.arange(t.n_atoms-2), np.insert(dij/l2,0,1)

def calcRg(df,t,fasta):
    masses = df.loc[fasta,'MW'].values
    masses[0] += 2
    masses[-1] += 16
    cm = np.sum(t.xyz*masses[np.newaxis,:,np.newaxis],axis=1)/masses.sum()
    si = np.linalg.norm(t.xyz - cm[:,np.newaxis,:],axis=2)
    rgarray = np.sqrt(np.sum(si**2*masses,axis=1)/masses.sum())
    rg = rgarray.mean()
    rgE = np.std(rgarray)/np.sqrt(rgarray.size)
    return rgarray, rg, rgE

residues = pd.read_csv('residues.csv',index_col='one')
t0 = time.time()

data = {}
f = lambda x,R0,v : R0*np.power(x,v)
for name in proteins.index:
    traj = md.load_dcd("{:s}/{:s}.dcd".format(name,name),"{:s}/{:s}.pdb".format(name,name))
    rgarray, _, _ = calcRg(residues,traj,proteins.loc[name].fasta)
    traj = traj[5000:]

    f = lambda x,k : np.exp(-x/k)
    ij, cij = calcCs(traj)
    popt, pcov = curve_fit(f,ij,cij)
    lp = 0.38*popt[0]
    lp_E = 0.38*pcov[0,0]**.5

    f = lambda x,R0,v : R0*np.power(x,v)
    ij,dij = calcRs(traj)
    popt, pcov = curve_fit(f,ij[ij>10],dij[ij>10],p0=[.4,.5])
    nu = popt[1]
    nu_E = pcov[1,1]**.5

    f = lambda x, v : lp*np.power(x,v)
    popt, pcov = curve_fit(f,ij[ij>10],dij[ij>10],p0=[.5])
    nu_lp = popt[0]
    nu_lp_E = pcov[0]**.5

    _, rg, rg_E = calcRg(residues,traj,proteins.loc[name].fasta)
    data[name] = {'ij':ij, 'dij':dij, 'nu':nu, 'nu_E':nu_E, 'lp':lp, 'lp_E':lp_E,
            'rg':rg, 'rg_E':rg_E, 'rgarray':rgarray, 'nu_lp':nu_lp, 'nu_lp_E':nu_lp_E}

df = pd.DataFrame(data=data).T

df.to_pickle('calc.pkl')
print('Timing {:.3f}'.format(time.time()-t0))
