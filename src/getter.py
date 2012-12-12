import global_stuff

import wc
import param
import objects
import sys

import ssh

import pdb

print  sys.argv

which_job = int(sys.argv[1])
total_jobs = int(sys.argv[2])

protein_list = sys.argv[3]

try:
    to_skip = []
    skip_file = sys.argv[4]
    h = open(skip_file, 'r')
    for line in h:
        to_skip.append(line.strip())
except:
    skip_file = None

username = 'fultonw'
password = 'lc7140$$'
hostname = 'ent.csail.mit.edu'
port = 22

#remote_folder = '/mnt/work/fw27/deleterious/deleterious/data/proteins/humvar/'
remote_base_folder = '/home/fultonw/scratch/'

#f = open(global_stuff.protein_list_file, 'r')
f = open(protein_list, 'r')

print 'starting', which_job, total_jobs

to_send = True
whether_to_delete = True

i = 0

import datetime
import helper
import os
#p = param.param({'ev':1e-10, 'protein_list_file':'asdfdone', 'uniprot_id':'Q9NVL1', 'avg_deg':1, 'n_cutoff':0, 'f_cutoff':15, 'which_msa':1, 'which_weight':1, 'which_dist':1, 'pseudo_c':1})

p = global_stuff.get_param()

past = datetime.datetime.now()

def get(obj, p, gotten_stuff, used_ps, check = True):
    global past
    used_ps.add(p.get_copy())
    print 'starting: ', p.get_param('uniprot_id'), obj, which_job, total_jobs
    to_get = True
    if check:
        gotten_stuff.append([obj, p.get_copy()])
        if wc.get_wrapper_instance(obj).has(p, False, True):
            to_get = False
    if to_get:
        ans = wc.get_stuff(obj, p, False, False, False)
        print 'took: ', datetime.datetime.now() - past
        past = datetime.datetime.now()
        gotten_stuff.append([obj, p.get_copy()])
        return ans
    else:
        gotten_stuff.append([obj, p.get_copy()])
        print 'already have: ', p.get_param('uniprot_id'), obj

import pdb


#specify the files that should be remote.  specify the files that should be deleted

used_ps = set()

#this is the stuff to send over.  delete these when u send them over
to_gets = set([objects.general_msa, objects.general_seq_weights, objects.neighbors_w_weight_w, objects.edge_to_rank, objects.general_distance])

#this is the stuff to delete right after one protein is processed
to_deletes = set([objects.general_msa, objects.general_seq_weights, objects.neighbors_w_weight_w, objects.edge_to_rank, objects.dW, objects.adW, objects.afW, objects.agW, objects.their_agW, objects.pairwise_dist, objects.general_distance, objects.mf_distance, objects.general_msa, objects.div_weights, objects.general_seq_weights]) - to_gets

sender = helper.file_sender(global_stuff.lock_folder + str(which_job % 5), 20)

for line in f:

    
    
    if i % total_jobs == which_job:


        gotten_stuff = []

        
        protein_name = line.strip()
        
        if protein_name not in to_skip:




            i += 1
            print '                        TRAVERSED:', i

            p.set_param('uniprot_id',protein_name)

            seq = get(objects.dW, p, gotten_stuff, used_ps, False)



            for which_weight in range(2):
                p.set_param('which_weight',which_weight)

                if len(seq) < 10000:
                    import pdb

                    try:
                        for to_get in to_gets:
                            get(to_get, p, gotten_stuff, used_ps)
              
                    except Exception, err:
                        print 'fail', protein_name, err
                else:
                    print 'skipping', protein_name

            if to_send:

                # move stuff to ent
                # first try to create the folder
                there_folder = remote_base_folder + protein_name + '/'

                for gotten in gotten_stuff:
                    try:

                        obj = gotten[0]
                        p_used = gotten[1]
                        instance = wc.get_wrapper_instance(obj)
                        here_file = instance.get_file_location(p_used)
                        there_folder = remote_base_folder + p_used.get_param('uniprot_id') + '/'
                        file_name = instance.get_file_name(p_used)
                        there_file = there_folder + file_name
                        sender.send(here_file, there_file, hostname, there_folder, username, password, port, instance, p_used)
                    except Exception, err:
                        pass


            if whether_to_delete:
                for to_delete in to_deletes:
                    for used_p in used_ps:
                        the_file = wc.get_wrapper_instance(to_delete).get_file_location(used_p)
                        try:
                            import subprocess

                            if os.path.isfile(the_file):
                                print '              blind removing:', the_file
                                subprocess.call(['rm', the_file])

                        except Exception, err:
                            #print err
                            pass
        else:
            print 'skipping ', protein_name
    i += 1

print 'ending', which_job, total_jobs
        
