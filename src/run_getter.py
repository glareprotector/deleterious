import sys
import subprocess

getter_script = 'getter.py'

total_jobs = int(sys.argv[1])

try:
    skip_file = sys.argv[2]
except:
    skip_file = None

for i in range(total_jobs):
    print i
    if skip_file != None:
        subprocess.Popen(['python', getter_script, str(i), str(total_jobs), skip_file], stdout = open('../dump/'+str(i)+'_out','w'), stderr = open('../dump/'+str(i)+'_err','w'))
    else:
        subprocess.Popen(['python', getter_script, str(i), str(total_jobs)], stdout = open('../dump/'+str(i)+'_out','w'), stderr = open('../dump/'+str(i)+'_err','w'))
