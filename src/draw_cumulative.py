import sys
from matplotlib.backends.backend_pdf import PdfPages, PdfFile
from pylab import *
import matplotlib.pyplot as plt
import pdb

out_file = sys.argv[1]

def read_score(in_file):
    f=open(in_file,'r')
    # first line holds possible labels
    scores = {}
    for line in f:
        s = line.strip().split(',')
        label = s[1]
        score = float(s[0])
        try:
            scores[label].append(score)
        except:
            scores[label] = [score]
    return scores


in_files = sys.argv[2:]

pdf = PdfFile(out_file)

for in_file in in_files:
    scores = read_score(in_file)
    pdb.set_trace()
    for label in scores:
        plt.hist(scores[label], cumulative=True, histtype='step', bins=1000, label=in_file + '_' + label)

plt.legend()
savefig(out_file)
close()
