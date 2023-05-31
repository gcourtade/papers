import pandas as pd
import numpy as np
import mdtraj as md
import itertools

def initProteinsDF():
    df_seq = pd.read_csv('linker_sequences.csv',index_col='seq_id')
    proteins = pd.DataFrame(index=df_seq.index,columns=['eps_factor','temp','pH','ionic','fasta'])
    for seq_id in df_seq.index:
        proteins.loc[seq_id] = dict(eps_factor=0.2,temp=298,pH=6.5,ionic=0.1,fasta=list(df_seq.loc[seq_id][0]))
    return proteins

def genParamsLJ(df,name,prot):
    fasta = prot.fasta.copy()
    r = df.copy()
    types = list(np.unique(fasta))
    sigmamap = pd.DataFrame((r.sigmas.values+r.sigmas.values.reshape(-1,1))/2,
                            index=r.sigmas.index,columns=r.sigmas.index)
    lambdamap = pd.DataFrame((r.lambdas.values+r.lambdas.values.reshape(-1,1))/2,
                             index=r.lambdas.index,columns=r.lambdas.index)
    lj_eps = prot.eps_factor*4.184
    # Generate pairs of amino acid types
    pairs = np.array(list(itertools.combinations_with_replacement(types,2)))
    return pairs, lj_eps, lambdamap, sigmamap, fasta, types

def genParamsDH(df,name,prot):
    kT = 8.3145*prot.temp*1e-3
    r = df.copy()
    # Set the charge on HIS based on the pH of the protein solution
    r.loc['H','q'] = 1. / ( 1 + 10**(prot.pH-6) )
    # Calculate the prefactor for the Yukawa potential
    qq = pd.DataFrame(r.q.values*r.q.values.reshape(-1,1),index=r.q.index,columns=r.q.index)
    fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3
    epsw = fepsw(prot.temp)
    lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.022*1000/kT
    yukawa_eps = qq*lB*kT
    # Calculate the inverse of the Debye length
    yukawa_kappa = np.sqrt(8*np.pi*lB*prot.ionic*6.022/10)
    return yukawa_eps, yukawa_kappa
