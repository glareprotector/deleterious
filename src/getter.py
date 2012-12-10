import global_stuff

import wc
import param
import objects
import sys

import ssh

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

username = 'fultonw'
password = 'lc7140$$'
hostname = 'dragon.csail.mit.edu'
port = 20

#remote_folder = '/mnt/work/fw27/deleterious/deleterious/data/proteins/humvar/'
remote_folder = '~/scratch/'

f = open(global_stuff.protein_list_file, 'r')

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


            if i > 2:
                break

            i += 1
            

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

            pdb.set_trace()

            # move stuff to ent
            # first try to create the folder
            folder = remote_folder + protein_name
            client = ssh.SSHClient()
            client.load_system_host_keys()
            client.set_missing_host_key_policy(ssh.AutoAddPolicy())
            client.connect(hostname, port, username, password)
            client.exec_command('mkdir ' + folder)


            cur_folder = wc.get_wrapper_instance(objects.by_uniprot_id_wrapper).get_folder(p)
            
            scp = ssh.SCPSclient(ssh.get_transport())

            all_files = os.listdir(cur_folder)
            for a_file in all_files:
                scp.put(cur_folder + a_file, folder + a_file)
                

            
        else:
            print 'skipping ', protein_name
    i += 1

print 'ending', which_job, total_jobs
        
