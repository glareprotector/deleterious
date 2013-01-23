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
whether_to_check_remote T or F
uniprot_or_pdb_chain U or P
fqueue(only for bsub)
mem_limit(only for bsub)
hours:minutes(only for bsub)
tmp space
then possibly param arguments to pass to getter.  param name, value, param type
"""
import os
import global_stuff
import pdb


print 'remember to source profile'
pdb.set_trace()

#a=subprocess.call('source ~/.profile', shell=False, executable='/bin/bash')

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
whether_to_check_remote = sys.argv[8]
uniprot_or_pdb_chain = sys.argv[9]

if mode == 'b':

    queue = sys.argv[10]
    mem_limit = sys.argv[11]
    time_limit = sys.argv[12]
    temp_mem_limit = sys.argv[13]
import pdb

if mode == 'b':
    arg_start = 14
elif mode == 'r':
    arg_start = 10

args = sys.argv[arg_start:]
import string
arg_string = string.join(args,sep=' ')

skip_file = None

import pdb



for i in range(total_jobs):
    print >> sys.stderr, i

    cmd = 'python ' + getter_script + ' ' + str(i) + ' ' + str(total_jobs) + ' ' + protein_list + ' ' + whether_to_send + ' ' + whether_to_delete + ' ' + whether_to_get_anything + ' ' + whether_to_temp + ' ' + whether_to_check_remote + ' ' + uniprot_or_pdb_chain + ' ' + arg_string
    print >> sys.stderr, cmd
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
            print >> sys.stderr, mem_limit, 'fdsa', time_limit, temp_mem_limit, log_file, error_file
            import string
            bsub_cmd = string.join(['bsub', '-q', queue, '-o', log_file, '-e', error_file,'-W', time_limit, '-R', '\"rusage[mem='+mem_limit+']' + ' && ' + 'rusage[tmp='+temp_mem_limit+']\"', cmd],sep = ' ')
            subprocess.call(bsub_cmd, shell=True, executable='/bin/bash')
