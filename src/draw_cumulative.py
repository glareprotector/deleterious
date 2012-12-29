import sys
from matplotlib.backends.backend_pdf import PdfPages, PdfFile
from pylab import *
import matplotlib.pyplot as plt

out_file = sys.argv[1]

def read_score(in_file):
    f=open(in_file,'r')
    scores = []
    for line in f:
        scores.append(float(line.strip()))
    return scores

in_files = sys.argv[2:]

pdf = PdfFile(out_file)

for in_file in in_files:
    scores = read_score(in_file)
    plt.hist(scores, cumulative=True, histtype='step', bins=1000, label=in_file)

plt.legend()
savefig(infile)
close()
