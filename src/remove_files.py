import global_stuff
import os
import pdb
import sys
import subprocess

in_list = sys.argv[1]

searching_for = sys.argv[2:]

f = open(in_list, 'r')

completed = []
i = 0
bad = 0
for line in f:
    try:
        name = line.strip()

        folder = global_stuff.base_folder + name + '/'
        files = os.listdir(folder)
        num_present = 0
        for a_file in files:
            for it in searching_for:
                if it in a_file:
                    the = folder + a_file
                    print the
                    subprocess.call('rm \'' + the + '\'', shell=True, executable='/bin/bash')
                    print i
                    i += 1
                    break
    except:
        bad += 1
        print bad, 'bad'



