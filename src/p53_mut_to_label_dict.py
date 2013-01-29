import global_stuff, pdb

f = open(global_stuff.P53_MESSY_MUTATIONS, 'r')
f.next()

m = {}
seen = set()
for line in f:

    b = line.strip().split('\t')
    mut = b[2]
    loss = b[10]
    gain = b[12]
    seen.add(mut)

    if loss != 'NA':
        try:
            m[mut].append(loss)
        except:
            m[mut] = [loss]
    
    if gain != 'N1':
        try:
            m[mut].append(gain)
        except:
            m[mut] = [gain]

pdb.set_trace()
