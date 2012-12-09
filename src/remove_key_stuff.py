import global_stuff
import os
import sys
import subprocess

files = os.listdir(global_stuff.BIN_FOLDER)

to_remove = sys.argv[1:]

for name in files:
    for it in to_remove:
        if it in name:
            the = global_stuff.BIN_FOLDER + name
            subprocess.call('rm \'' + the + '\'', shell=True, executable='/bin/bash')
    
