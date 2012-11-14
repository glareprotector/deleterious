import sys
import subprocess

getter_script = 'getter.py'

total_jobs = int(sys.argv[1])

for i in range(total_jobs):
    print i
    subprocess.Popen(['python', getter_script, str(i), str(total_jobs)], stdout = open('../dump/'+str(i)+'_out','w'), stderr = open('../dump/'+str(i)+'_err','w'))
