import global_stuff

import wc
import param
import objects
import sys

which_job = int(sys.argv[1])
total_jobs = int(sys.argv[2])

which_object = objects.dW

f = open(global_stuff.protein_list_file, 'r')

i = 0
for line in f:
    if i % total_jobs == which_job:
        protein_name = line.strip()
        wc.get_stuff(which_obj, param.param({'uniprot_id':protein_name}), True, True, False)


