import global_stuff

import wc
import param
import objects
import sys

which_job = int(sys.argv[1])
total_jobs = int(sys.argv[2])

which_obj = objects.easy_to_read_msa_format

evalue = 1e-10

f = open(global_stuff.protein_list_file, 'r')

which_completed_file = '../dump/which_completeds/job' + str(which_job) + '_' + str(total_jobs)
g = open(which_completed_file, 'a')


print 'starting', which_job, total_jobs

i = 0

for line in f:
    if i % total_jobs == which_job:
        protein_name = line.strip()
        print which_job, total_jobs, 'running ', i
        wc.get_stuff(which_obj, param.param({'uniprot_id':protein_name, 'ev':evalue}), True, True, False)
        g.write(protein_name + '\n')
    i += 1

print 'ending', which_job, total_jobs
        
