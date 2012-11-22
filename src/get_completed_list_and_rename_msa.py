import global_stuff
import os
import subprocess
import wc
import objects
import param

f = open(global_stuff.protein_list_file, 'r')

completed = []

evalue = 1e-10

for line in f:
    name = line.strip()
    folder = global_stuff.base_folder + name + '/'
    files = os.listdir(folder)
    has_easy = False
    has_dist = False
    enough_rows = False
    for a_file in files:
        if 'easy' in a_file:
            has_easy = True
            subprocess.call(['cp', folder+a_file, folder+'msa'])
        if 'pairwise' in a_file:
            has_dist = True
            subprocess.call(['cp', folder+a_file, folder+'dists'])

            # copy to better file_name
    msa = wc.get_stuff(objects.agW, param.param({'uniprot_id':name, 'ev':evalue}), False, False, False)
    if len(msa) > 50:
        enough_rows = True

    if has_easy and has_dist and enough_rows:
        completed.append(name)

g = open(global_stuff.completed_list_file, 'w')
for name in completed:
    g.write(name + '\n')

f.close
g.close()
