import global_stuff

import wc
import param
import objects
import sys

import pdb

which_job = int(sys.argv[1])
total_jobs = int(sys.argv[2])


try:
    to_skip = []
    skip_file = sys.argv[3]
    h = open(skip_file, 'r')
    for line in h:
        to_skip.append(line.strip())
except:
    skip_file = None



which_obj = objects.pairwise_dist

evalue = 1e-10

f = open(global_stuff.protein_list_file, 'r')

which_completed_file = '../dump/which_completeds/job' + str(which_job) + '_' + str(total_jobs)
g = open(which_completed_file, 'a')


print 'starting', which_job, total_jobs

i = 0

import datetime


#p = param.param({'ev':1e-10, 'protein_list_file':'asdfdone', 'uniprot_id':'Q9NVL1', 'avg_deg':1, 'n_cutoff':0, 'f_cutoff':15, 'which_msa':1, 'which_weight':1, 'which_dist':1, 'pseudo_c':1})

p = global_stuff.get_param()

past = datetime.datetime.now()

def get(obj, p):
    global past
    print 'starting: ', p.get_param('uniprot_id'), obj, which_job, total_jobs
    ans = wc.get_stuff(obj, p, False, False, False)
    print 'took: ', datetime.datetime.now() - past
    past = datetime.datetime.now()
    return ans

for line in f:


    
    if i % total_jobs == which_job:
        protein_name = line.strip()
        
        if protein_name not in to_skip:


            p.set_param('uniprot_id',protein_name)
            seq = get(objects.dW, p)


            for which_weight in range(2):
                p.set_param('which_weight',which_weight)
            
                if len(seq) < 500:
                    import pdb

                    try:
                        get(objects.general_msa,p)
                        get(objects.general_seq_weights,p)
                        get(objects.neighbors_w_weight_w,p)
                        get(objects.edge_to_rank,p)
                    except:
                        print 'fail', protein_name
                else:
                    print 'skipping', protein_name
        else:
            print 'skipping ', protein_name
    i += 1

print 'ending', which_job, total_jobs
        
