import global_stuff
import os
import pdb
import sys
import subprocess

in_list = sys.argv[1]

searching_for = sys.argv[2:]

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
                the = folder + a_file
                subprocess.call('rm \'' + the + '\'', shell=True, executable='/bin/bash')

                break



