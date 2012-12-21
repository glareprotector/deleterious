import sys
import subprocess

"""
total jobs
protein_list absolute path
mode
whether_to_send T or F
whether_to_delete T or F
queue(only for bsub)
mem_limit(only for bsub)
hours:minutes(only for bsub)
then possibly param arguments to pass to getter
"""

getter_script = 'getter.py'

log_dir = '../dump/'

total_jobs = int(sys.argv[1])

protein_list = sys.argv[2]

mode = sys.argv[3]

protein_list_file = protein_list.split('/')[-1]

whether_to_send = sys.argv[4]
whether_to_delete = sys.argv[5]

if mode == 'b':

    queue = sys.argv[6]
    mem_limit = sys.argv[7]
    time_limit = sys.argv[8]
import pdb

if mode == 'b':
    arg_start = 9
elif mode == 'r':
    arg_start = 6

args = sys.argv[arg_start:]
import string
arg_string = string.join(args,sep=' ')

skip_file = None

import pdb



for i in range(total_jobs):
    print i

    cmd = 'python ' + getter_script + ' ' + str(i) + ' ' + str(total_jobs) + ' ' + protein_list + ' ' + whether_to_send + ' ' + whether_to_delete + ' ' + arg_string
    print cmd
    if mode == 'b':
        log_file = log_dir + 'b'+ '_' + str(i) + '_' + str(total_jobs) + '_' + protein_list_file + '_' + 'log'
        error_file = log_dir + 'b' + '_' + str(i) + '_' + str(total_jobs) + '_' + protein_list_file + '_' + 'err'
    if mode == 'r':
        log_file = log_dir + 'r'+ '_' + str(i) + '_' + str(total_jobs) + protein_list_file + '_' + '_log'
        error_file = log_dir + 'r' + '_' + str(i) + '_' + str(total_jobs) + protein_list_file + '_' + '_err'
        

    
    if skip_file != None:
        subprocess.Popen(['python', getter_script, str(i), str(total_jobs), skip_file], stdout = open(log_file,'w'), stderr = open(error_file,'w'))
    else:
        if mode == 'r':
            #subprocess.Popen(cmd.split())
            subprocess.Popen(cmd.split(), stdout = open(log_file,'w'), stderr = open(error_file,'w'))
        if mode == 'b':
            subprocess.call(['bsub', '-R', 'rusage[mem='+mem_limit+']', '-o', log_file, '-e', error_file,'-W', time_limit, '-q', queue, cmd])
            
