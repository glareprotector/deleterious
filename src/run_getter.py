import sys
import subprocess

getter_script = 'getter.py'

log_dir = '../dump/'

total_jobs = int(sys.argv[1])

mode = sys.argv[2]

try:
    queue = sys.argv[3]
except:
    pass

try:
    mem_limit = sys.argv[4]
except:
    pass

try:
    skip_file = sys.argv[2]
except:
    skip_file = None

for i in range(total_jobs):
    print i

    cmd = 'python ' + getter_script + ' ' + str(i) + ' ' + str(total_jobs)

    log_file = log_dir + str(i) + '_' + str(total_jobs) + '_log'
    error_file = log_dir + str(i) + '_' + str(total_jobs) + '_err'
    
    if skip_file != None:
        subprocess.Popen(['python', getter_script, str(i), str(total_jobs), skip_file], stdout = open('../dump/'+str(i)+'_out','w'), stderr = open('../dump/'+str(i)+'_err','w'))
    else:
        if mode == 'r':
            subprocess.Popen(['python', getter_script, str(i), str(total_jobs)], stdout = open('../dump/'+str(i)+'_out','w'), stderr = open('../dump/'+str(i)+'_err','w'))
        if mode == 'b':
            subprocess.call(['bsub', '-R', 'rusage[mem='+mem_limit+']', '-o', log_file, '-e', error_file, '-q', queue, cmd])
            
