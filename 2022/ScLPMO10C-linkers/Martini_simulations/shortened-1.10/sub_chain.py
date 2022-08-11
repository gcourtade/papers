#!/usr/bin/python

import sys
import os


partition='sbinlab'

slurm=sys.argv[1] #run.slurm
try:
	ch_len=int(sys.argv[2])
	print 'Requested to submit %s jobs to be chained!'%str(ch_len)
except IndexError:
	ch_len=5
	print 'Chain length not give - default 5!'

try:
	jid=sys.argv[3]
	print 'Submitted jobs will be dependent in jid %s'%jid
except IndexError:
	jid=0
	print 'No jid supplied, will assume first submission!'
	
ans=raw_input('Last chance to check submission info... \nIs %s the correct partition (y/n)?'%partition)

if ans=='y':
	print "Let's go then!"
elif ans is not 'y':
	partition=raw_input('Please set correct partition for seetings in %s: '%slurm,)



jid_list=[]
def sub_j(slurm,jid):
        if jid==0:
                os.system('sbatch --partition '+partition+' '+slurm+' > jid')
                with open('jid','r') as jobid:
                        for line in jobid:
                                jid=line.split()[3]
                        print 'Submitted jid '+jid+' as first job!'

                #sub first job of series
                #save jid (+ print)
        else:
                os.system('sbatch --partition '+partition+' --dependency=afterany:'+jid+' '+slurm+'>jid')
                jid_d=jid
                with open('jid','r') as jobid:
                        for line in jobid:
                                jid=line.split()[3]
                     
                        print 'Submitted jid '+jid+' as -d on '+jid_d
	return jid

def jid_list_out(jid_list):
	with open('jid.list','w') as out:
		for job in jid_list:
			out.write('%s \n'% (job))

if __name__=="__main__":
	for i in range(ch_len):
		jid=sub_j(slurm,jid)
		jid_list.append(jid)
	jid_list_out(jid_list)
