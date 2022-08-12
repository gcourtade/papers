from analyse import *
import hoomd
import hoomd.md
import time
import os
import sys
import psutil
from argparse import ArgumentParser
from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB

parser = ArgumentParser()
parser.add_argument('--name',nargs='?',const='', type=str)
args = parser.parse_args()

def simulate(residues,name,prot):
    Nsteps = 7e7
    hoomd.context.initialize("--mode=cpu --nthreads=1");
    hoomd.option.set_notice_level(0)
    hoomd.util.quiet_status()
    pairs, lj_eps, lj_lambda, lj_sigma, fasta, types = genParamsLJ(residues,name,prot)
    yukawa_eps, yukawa_kappa = genParamsDH(residues,name,prot)
    N = len(fasta)
    snapshot = hoomd.data.make_snapshot(N=N,
                                    box=hoomd.data.boxdim(Lx=200, Ly=200, Lz=200),
                                    particle_types=types,
                                    bond_types=['polymer']);

    geo = Geometry.geometry(prot.fasta[0])
    geo.phi = -120
    geo.psi_im1 = 150
    structure = PeptideBuilder.initialize_res(geo)
    for residue in prot.fasta[1:]:
        structure = PeptideBuilder.add_residue(structure, residue, geo.phi, geo.psi_im1)
    out = Bio.PDB.PDBIO()
    out.set_structure(structure)
    xyz = []
    for atom in out.structure.get_atoms():
        if atom.name == 'CA':
            xyz.append([atom.coord[0]/10.,atom.coord[1]/10.,atom.coord[2]/10.])
    snapshot.particles.position[:] = xyz;

    fst = ''.join(prot.fasta)
    ids = [types.index(a) for a in fasta]
    snapshot.particles.typeid[:] = ids
    snapshot.particles.mass[:] = [residues.loc[a].MW for a in prot.fasta]
    snapshot.particles.mass[0] += 2
    snapshot.particles.mass[-1] += 16
    snapshot.bonds.resize(N-1);
    snapshot.bonds.group[:] = [[i,i+1] for i in range(N-1)];
    snapshot.bonds.typeid[:] = [0] * (N-1)
    hoomd.init.read_snapshot(snapshot);

    kT = 8.3145*prot.temp*1e-3
    hb = hoomd.md.bond.harmonic();
    hb.bond_coeff.set('polymer', k=8033.0, r0=0.38);

    nl = hoomd.md.nlist.cell();

    lj1 = hoomd.md.pair.lj(r_cut=4.0, nlist=nl, name="1")
    lj2 = hoomd.md.pair.lj(r_cut=4.0, nlist=nl, name="2")
    yukawa = hoomd.md.pair.yukawa(r_cut=4.0, nlist=nl)
    for a,b in pairs:
        lj1.pair_coeff.set(a, b, epsilon=lj_eps*(1-lj_lambda.loc[a,b]), sigma=lj_sigma.loc[a,b],
                    r_cut=np.power(2.,1./6)*lj_sigma.loc[a,b])
        lj2.pair_coeff.set(a, b, epsilon=lj_eps*lj_lambda.loc[a,b], sigma=lj_sigma.loc[a,b], r_cut=4.)
        yukawa.pair_coeff.set(a, b, epsilon=yukawa_eps.loc[a,b], kappa=yukawa_kappa, r_cut=4.);
    lj1.set_params(mode='shift')
    yukawa.set_params(mode='shift')
    nl.reset_exclusions(exclusions = ['bond'])

    integrator_mode = hoomd.md.integrate.mode_standard(dt=0.005);
    integrator = hoomd.md.integrate.langevin(group=hoomd.group.all(),kT=kT,seed=np.random.randint(100));
    for a in types:
        integrator.set_gamma(a, residues.loc[a].MW/100)
    nl.tune()
    hoomd.dump.gsd("{:s}/{:s}.gsd".format(name,name), period=2e3, group=hoomd.group.all(), overwrite=True);
    hoomd.dump.gsd("{:s}/restart.gsd".format(name), group=hoomd.group.all(), truncate=True, period=1e6, phase=0)
    hoomd.run(Nsteps)

def genDCD(residues,name,prot):
    """
    Generates coordinate and trajectory
    in convenient formats
    """
    top = md.Topology()
    chain = top.add_chain()
    for resname in prot.fasta:
        residue = top.add_residue(residues.loc[resname,'three'], chain)
        top.add_atom(residues.loc[resname,'three'], element=md.element.carbon, residue=residue)
    for i in range(len(prot.fasta)-1):
        top.add_bond(top.atom(i),top.atom(i+1))
    traj = md.load("{:s}/{:s}.gsd".format(name,name), top)
    traj.top = top
    traj = traj.image_molecules(inplace=False, anchor_molecules=[set(traj.top.chain(0).atoms)], make_whole=True)
    traj.center_coordinates()
    traj.xyz += traj.unitcell_lengths[0,0]/2
    traj.save_dcd("{:s}/{:s}.dcd".format(name,name))
    traj[0].save_pdb("{:s}/{:s}.pdb".format(name,name))
    os.remove("{:s}/{:s}.gsd".format(name,name))

proteins = pd.read_pickle('proteins.pkl')
residues = pd.read_csv('residues.csv',index_col='one')

t0 = time.time()
simulate(residues,args.name,proteins.loc[args.name])
genDCD(residues,args.name,proteins.loc[args.name])
print('Timing Simulation {:.3f}'.format(time.time()-t0))
