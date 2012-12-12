import sys
import subprocess

getter_script = 'getter.py'

log_dir = '../dump/'

total_jobs = int(sys.argv[1])

protein_list = sys.argv[2]

mode = sys.argv[3]

try:
    queue = sys.argv[4]
except:
    pass

try:
    mem_limit = sys.argv[5]
except:
    pass

try:
    time_limit = sys.argv[6]
except:
    pass

skip_file = None

import pdb



for i in range(total_jobs):
    print i

    cmd = 'python ' + getter_script + ' ' + str(i) + ' ' + str(total_jobs) + ' ' + protein_list
    print cmd
    if mode == 'b':
        log_file = log_dir + 'b'+ '_' + str(i) + '_' + str(total_jobs) + '_log'
        error_file = log_dir + 'b' + '_' + str(i) + '_' + str(total_jobs) + '_err'
    if mode == 'r':
        log_file = log_dir + 'r'+ '_' + str(i) + '_' + str(total_jobs) + '_log'
        error_file = log_dir + 'r' + '_' + str(i) + '_' + str(total_jobs) + '_err'
        

    
    if skip_file != None:
        subprocess.Popen(['python', getter_script, str(i), str(total_jobs), skip_file], stdout = open(log_file,'w'), stderr = open(error_file,'w'))
    else:
        if mode == 'r':
            #subprocess.Popen(cmd.split())
            subprocess.Popen(cmd.split(), stdout = open(log_file,'w'), stderr = open(error_file,'w'))
        if mode == 'b':
            subprocess.call(['bsub', '-R', 'rusage[mem='+mem_limit+']', '-o', log_file, '-e', error_file,'-W', time_limit, '-q', queue, cmd])
            
