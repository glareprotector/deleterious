import global_stuff
import os
import pdb
import sys

in_list = sys.argv[1]
out_list = sys.argv[2]

searching_for = sys.argv[3:]

f = open(global_stuff.data_folder + in_list, 'r')

completed = []

for line in f:
    name = line.strip()
    folder = global_stuff.base_folder + name + '/'
    files = os.listdir(folder)
    num_present = 0
    for a_file in files:
        for it in searching_for:
            if it in a_file:
                num_present += 1

    if num_present == len(searching_for):
        completed.append(name)
    

g = open(global_stuff.data_folder + out_list, 'w')
for name in completed:
    g.write(name + '\n')

f.close
g.close()
