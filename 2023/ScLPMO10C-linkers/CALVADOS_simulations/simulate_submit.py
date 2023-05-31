from analyse import *
import os
import subprocess
from jinja2 import Template

proteins = initProteinsDF()
proteins.to_pickle('proteins.pkl')

submission = Template("""#!/bin/bash
#SBATCH --job-name={{name}}
#SBATCH --nodes=1           
#SBATCH --partition=sbinlab
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH -t 16:00:00
#SBATCH -o {{name}}.out
#SBATCH -e {{name}}.err
#SBATCH --mem-per-cpu=2000

source /groups/sbinlab/giulio/.bashrc
conda activate hoomd

python ./simulate.py --name {{name}}""")

for name in proteins.index:
    if not os.path.isdir(name):
        os.mkdir(name)
    with open('{:s}.sh'.format(name), 'w') as submit:
        submit.write(submission.render(name=name))
    subprocess.run(['sbatch','--partition','sbinlab','{:s}.sh'.format(name)])
