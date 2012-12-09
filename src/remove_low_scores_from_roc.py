import sys

num = len(sys.argv)-1
files = sys.argv[1:]
fs = [open(files[i], 'r') for i in range(len(files))]

lines = [f.readlines() for f in fs]

outfiles = [x+'new' for x in files]
gs = [open(outfiles[i],'w') for i in range(len(outfiles))]

for i in range(len(lines[0])):
    s = lines[0][i].strip().split(',')
    score = float(s[0])
    if abs(score) > .0001:
        for j in range(len(fs)):
            gs[j].write(lines[j][i])

"""

first = sys.argv[1]
second = sys.argv[2]

f = open(first, 'r')
g = open(second, 'r')

first_out = first + 'new'
second_out = second + 'new'

firstlines = f.readlines()
secondlines = g.readlines()

assert len(firstlines) == len(secondlines)

fo = open(first_out, 'w')
go = open(second_out, 'w')

for i in range(len(firstlines)):

    s = firstlines[i].strip().split(',')
    score = float(s[0])
    print score
    if score > .00001:
        fo.write(firstlines[i])
        go.write(secondlines[i])
"""
