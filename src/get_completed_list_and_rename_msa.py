import global_stuff
import os
import subprocess

f = open(global_stuff.protein_list_file, 'r')

completed = []

for line in f:
    name = line.strip()
    folder = global_stuff.base_folder + name + '/'
    files = os.listdir(folder)
    for a_file in files:
        if 'easy' in a_file:
            completed.append(name)
            # copy to better file_name
            subprocess.call(['cp', folder+a_file, folder+'msa'])
            break
    

g = open(global_stuff.completed_list_file, 'w')
for name in completed:
    g.write(name + '\n')

f.close
g.close()
