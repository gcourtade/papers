import mdtraj as md
t = md.load('md_prot.xtc',top='all_LPMO_theta1.10.top')
print(t)
