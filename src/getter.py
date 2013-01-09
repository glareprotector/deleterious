"""
which_job
total_jobs
protein_list absolute path
whether_to_send T or F
whether_to_delete T or F
whether_to_get_anything T or F
whether_to_temp T or F
whether_to_check_remote T or F
optional: params followed by their values and type(i,f,s)
"""


import global_stuff

import wc
import param
import objects
import sys

import ssh
import os
import pdb
import shutil
import subprocess
print  sys.argv

which_job = int(sys.argv[1])
total_jobs = int(sys.argv[2])

protein_list = sys.argv[3]

if sys.argv[4] == 'T':
    whether_to_send = True
else:
    whether_to_send = False

if sys.argv[5] == 'T':
    whether_to_delete = True
else:
    whether_to_delete = False

if sys.argv[6] == 'T':
    whether_to_get_anything = True
else:
    whether_to_get_anything = False

if sys.argv[7] == 'T':
    whether_to_temp = True
else:
    whether_to_temp = False

if sys.argv[8] == 'T':
    whether_to_check_remote = True
else:
    whether_to_check_remote = False

skip_file = None

username = 'fultonw'
password = 'lc7140$$$$'
hostname = 'ent.csail.mit.edu'
port = 22

#remote_folder = '/mnt/work/fw27/deleterious/deleterious/data/proteins/humvar/'
#remote_base_folder = '/mnt/work/fultonw/scratch_cosmic/'
#remote_base_folder = '/mnt/work/fultonw/scratch/'

f = open(protein_list, 'r')

print >> sys.stderr, 'starting', which_job, total_jobs



i = 0

import datetime
import helper
import os
#p = param.param({'ev':1e-10, 'protein_list_file':'asdfdone', 'uniprot_id':'Q9NVL1', 'avg_deg':1, 'n_cutoff':0, 'f_cutoff':15, 'which_msa':1, 'which_weight':1, 'which_dist':0, 'pseudo_c':1})

p = global_stuff.get_param()

to_skip = []


# get param values.
print >> sys.stderr, sys.argv
assert (len(sys.argv)-9)%3 == 0
helper.parse_p_input(p, sys.argv[9:])



past = datetime.datetime.now()

def get(obj, p, gotten_stuff, used_ps, check = True):

    global past
    used_ps.add(p.get_copy())
    print >> sys.stderr, 'starting: ', p.get_param('uniprot_id'), obj, which_job, total_jobs
    to_get = True

    gotten_stuff.append([obj, p.get_copy()])
    if check:

        if wc.get_wrapper_instance(obj).has(p, False, True):
            to_get = False
    if to_get:
        global whether_to_get_anything
        if whether_to_get_anything:
            if not wc.get_wrapper_instance(obj).has(p, False, check):
                ans = wc.get_stuff(obj, p, False, False, False)

                print >> sys.stderr, 'took: ', datetime.datetime.now() - past
                past = datetime.datetime.now()
            
                return ans


    print >> sys.stderr, 'already have: ', p.get_param('uniprot_id'), obj


import pdb


#specify the files that should be remote.  specify the files that should be deleted

used_ps = set()

#this is the stuff to send over.  delete these when u send them over
#to_gets = set([objects.general_msa, objects.general_distance, objects.general_seq_weights, objects.edge_to_rank, objects.neighbors_w_weight_w])
import wrapper
to_gets = set([objects.general_msa, wrapper.my_msa_obj_wrapper, objects.general_distance, objects.neighbors_w_weight_w])

to_blind_sends = set([objects.general_msa, objects.general_distance, objects.neighbors_w_weight_w])


#to_gets = set()

#this is the stuff to delete right after one protein is processed
to_deletes = set([objects.general_msa, objects.general_seq_weights, objects.neighbors_w_weight_w, objects.edge_to_rank, objects.dW, objects.adW, objects.aeW, objects.afW, objects.agW, objects.their_agW, objects.pairwise_dist, objects.general_distance, objects.mf_distance, objects.div_weights, objects.general_seq_weights, objects.MIP_input_msa, objects.MIP_input_msa_file, objects.MIP_distance_file, objects.MIP_distance, objects.bW])# - to_gets


sender = helper.file_sender(global_stuff.lock_folder + str(which_job % 5), 0)

protein_list_file = protein_list.split('/')[-1]




log_file = global_stuff.process_folder + str(which_job) + '_' + str(total_jobs) + '_' + protein_list_file

import datetime
past = (datetime.datetime.now())
past2 = (datetime.datetime.now())
print >> sys.stderr, 'ggggggggggg', past
g = open(log_file,'w')

f2 = open(protein_list, 'r')
num_proteins = len(f2.readlines())

for line in f:

    
    
    if i % total_jobs == which_job:


        gotten_stuff = []

        
        protein_name = line.strip()
        p.set_param('uniprot_id',protein_name)

        
        
        import wc
        import pdb

        seq = wc.get_stuff(objects.dW,p)

        print >> sys.stderr, "currently getting: ", protein_name, len(seq)
        
        if len(seq) < 1000:




            if whether_to_temp:
                global_stuff.home = global_stuff.temp_home
                assert global_stuff.base_folder == global_stuff.real_base_folder

                try:
                    os.makedirs(global_stuff.temp_base_folder)
                except:
                    pass
                try:
                    os.makedirs(global_stuff.get_holding_folder())
                except:
                    pass
                real_uniprot_folder = wc.get_wrapper_instance(objects.dW).get_folder(p)
                global_stuff.base_folder = global_stuff.temp_base_folder
                temp_uniprot_folder = wc.get_wrapper_instance(objects.dW).get_folder(p)
                if os.path.isdir(temp_uniprot_folder):
                    shutil.rmtree(temp_uniprot_folder)
                mv_cmd = 'cp -r ' + real_uniprot_folder + ' ' + temp_uniprot_folder
                subprocess.call(mv_cmd, shell=True, executable='/bin/bash')




            print >> sys.stderr, 'hhhhhhhhhhhhhh', past
            g.write('started: ' + protein_name + ' ' + str(i) + ' out of ' + str(num_proteins) + ' by ' + str(total_jobs) + ' ' +  str(datetime.datetime.now()) + ' ' + str(datetime.datetime.now()-past2) + '\n')
            past2  = datetime.datetime.now()

            g.flush()
            

            print >> sys.stderr, '                        TRAVERSED:', i

            

            #print >> sys.stderr, 'seq length: ', len(seq)

            for which_filter_co in [0.2]:
                for avg_deg in [1,2,3,4,5,6,7,8]:

                    p.set_param('avg_deg', avg_deg)
                    p.set_param('filter_co',which_filter_co)
                    
                    import pdb

                    try:
                        for to_get in to_gets:

                            get(to_get, p, gotten_stuff, used_ps)

                    except Exception, err:
                        print >> sys.stderr, 'fail', protein_name, err


            g.write('finished: ' + protein_name + ' ' + str(i) + ' out of ' + str(num_proteins) + ' by ' + str(total_jobs) + ' ' +  str(datetime.datetime.now()) + ' ' + str(datetime.datetime.now()-past2) + '\n')
            past2 = datetime.datetime.now()
            g.flush()

            if whether_to_send:

                # move stuff to ent
                # first try to create the folder
                there_folder = global_stuff.remote_base_folder + protein_name + '/'


                # also try sending other stuff if it's there
                copy = []
                for gotten in gotten_stuff:
                    copy.append(gotten)
                    for to_blind_send in to_blind_sends:
                        copy.append([to_blind_send, gotten[1]])


                for gotten in copy:
                    try:

                        obj = gotten[0]
                        p_used = gotten[1]
                        instance = wc.get_wrapper_instance(obj)
                        here_file = instance.get_file_location(p_used)
                        there_folder = global_stuff.remote_base_folder + p_used.get_param('uniprot_id') + '/'
                        file_name = instance.get_file_name(p_used)
                        there_file = there_folder + file_name
                        import pdb
                        
                        sender.send(here_file, there_file, hostname, there_folder, username, password, port, instance, p_used, whether_to_delete)
                    except Exception, err:
                        pass


            if whether_to_delete:
                for to_delete in to_deletes:
                    for used_p in used_ps:
                        the_file = wc.get_wrapper_instance(to_delete).get_file_location(used_p)
                        try:
                            import subprocess
                            #print >> sys.stderr, the_file
                            if os.path.isfile(the_file):
                                print >> sys.stderr, '              blind removing:', the_file
                                subprocess.call(['rm', the_file])

                        except Exception, err:
                            #print >> sys.stderr, err
                            pass
            if whether_to_temp:
                #make copy of old folder, delete old folder, move temp folder, delete copy
                real_uniprot_folder_copy = real_uniprot_folder[:-1] + '.backup'
                if os.path.isdir(real_uniprot_folder_copy):
                    subprocess.call('rm -r ' + real_uniprot_folder_copy, shell=True, executable='/bin/bash')
                try:
                    subprocess.call('mv ' + real_uniprot_folder + ' ' + real_uniprot_folder_copy, shell=True, executable='/bin/bash')
                    shutil.move(temp_uniprot_folder, real_uniprot_folder[:-1])
                    subprocess.call('rm -r ' + real_uniprot_folder_copy, shell=True, executable='/bin/ash')
                except Exception, err:
                    print >> sys.stderr, Exception, err
                    
                global_stuff.base_folder = global_stuff.real_base_folder
                global_stuff.home = global_stuff.real_home
        else:
            print >> sys.stderr, 'skipping ', protein_name
    i += 1

print >> sys.stderr, 'ending', which_job, total_jobs
        
