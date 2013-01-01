import sys
import subprocess

"""
total jobs
protein_list absolute path
mode
whether_to_send T or F
whether_to_delete T or F
whether_to_get_anything T or F
whether to operate in temp folder T or F
queue(only for bsub)
mem_limit(only for bsub)
hours:minutes(only for bsub)
tmp space
then possibly param arguments to pass to getter
"""
import os
import global_stuff


getter_script = 'getter.py'

# make the folder, or delete its contents.  same with process folder
if os.path.isdir(global_stuff.log_folder):
    subprocess.call('rm ' + global_stuff.log_folder + '*', shell=True, executable='/bin/bash')
else:
    os.path.makedirs(global_stuff.log_folder)

if os.path.isdir(global_stuff.process_folder):
    subprocess.call('rm ' + global_stuff.process_folder + '*', shell=True, executable='/bin/bash')
else:
    os.makedirs(global_stuff.process_folder)

subprocess.call('rm ' + global_stuff.lock_folder + '*', shell=True, executable='/bin/bash')


total_jobs = int(sys.argv[1])

protein_list = sys.argv[2]

mode = sys.argv[3]

protein_list_file = protein_list.split('/')[-1]

whether_to_send = sys.argv[4]
whether_to_delete = sys.argv[5]
whether_to_get_anything = sys.argv[6]
whether_to_temp = sys.argv[7]

if mode == 'b':

    queue = sys.argv[8]
    mem_limit = sys.argv[9]
    time_limit = sys.argv[10]
    temp_mem_limit = sys.argv[11]
import pdb

if mode == 'b':
    arg_start = 12
elif mode == 'r':
    arg_start = 8

args = sys.argv[arg_start:]
import string
arg_string = string.join(args,sep=' ')

skip_file = None

import pdb



for i in range(total_jobs):
    print i

    cmd = 'python ' + getter_script + ' ' + str(i) + ' ' + str(total_jobs) + ' ' + protein_list + ' ' + whether_to_send + ' ' + whether_to_delete + ' ' + whether_to_get_anything + ' ' + whether_to_temp + ' ' + arg_string
    print cmd
    if mode == 'b':
        log_file = global_stuff.log_folder + 'b'+ '_' + str(i) + '_' + str(total_jobs) + '_' + protein_list_file + '_' + 'log'
        error_file = global_stuff.log_folder + 'b' + '_' + str(i) + '_' + str(total_jobs) + '_' + protein_list_file + '_' + 'err'
    if mode == 'r':
        log_file = global_stuff.log_folder + 'r'+ '_' + str(i) + '_' + str(total_jobs) + '_' + protein_list_file + '_' + 'log'
        error_file = global_stuff.log_folder + 'r' + '_' + str(i) + '_' + str(total_jobs) + '_' + protein_list_file + '_' + 'err'
        

    
    if skip_file != None:
        subprocess.Popen(['python', getter_script, str(i), str(total_jobs), skip_file], stdout = open(log_file,'w'), stderr = open(error_file,'w'))
    else:
        if mode == 'r':
            #subprocess.Popen(cmd.split())
            subprocess.Popen(cmd.split(), stdout = open(log_file,'w'), stderr = open(error_file,'w'))
        if mode == 'b':
            print mem_limit, 'fdsa', time_limit, temp_mem_limit, log_file, error_file
            import string
            bsub_cmd = string.join(['bsub', '-q', queue, '-o', log_file, '-e', error_file,'-W', time_limit, '-R', '\"rusage[mem='+mem_limit+']' + ' && ' + 'rusage[tmp='+temp_mem_limit+']\"', cmd],sep = ' ')
            subprocess.call(bsub_cmd, shell=True, executable='/bin/bash')
