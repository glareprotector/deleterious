import sys

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

