#!/public/home/fgl/public/miniconda3/bin/python3

import os
import re
import sys
from xml.etree.ElementTree import parse
 
def convert2gb(size):
    size = re.sub('\s+', '', size).replace(',', '').lower().rstrip('b') + 'b'
    regex = re.compile('^([0-9.]+)([kmgt]*)i*[b]$')
    assert regex.match(size), 'Invalid memory size: ' + size
    num, unit = regex.sub('\g<1> \g<2>b', size).split(' ')
    mags = ['b', 'kb', 'mb', 'gb', 'tb']
    mag = mags.index(unit) - 3
    num = round(float(num) * 1000 ** mag, 3)
    return '{:.3f}gb'.format(num)

def getNodesStatus():
    qs = [ re.split(' +', i) for i in os.popen('bqueues -w').readlines() ]
    qs = [ [i[16].split(','),i[1]] for i in qs[1:]]
    dq = {}
    for i in qs:
        for j in i[0]:
            try:
                dq[j] = dq[j] + ',' + i[1]
            except KeyError:
                dq[j] = i[1]
    dNodes = {}
    nodes = [ i.split(' ')[0] for i in os.popen('qnodes -l all').readlines() ]
    for node in nodes:
        attribs = list(list(parse(os.popen('qnodes {:s} -x'.format(node))).getroot())[0])
        njobs = -1
        for i in attribs:
            if i.tag == 'np':
                ncpus = int(i.text)
            elif i.tag == 'jobs':
                njobs = len(i.text.split(','))
            elif i.tag == 'status':
                avmem = [ j.split('availmem=')[1] for j in i.text.split(',') if j.startswith('availmem=')][0]
                avmem = convert2gb(avmem)
        if njobs == -1:
            njobs = 0
        nfree = ncpus - njobs
        dNodes[node] = [nfree, avmem, dq[node]]
    dNodes = dict(sorted(dNodes.items(), key=lambda x:int(re.sub('[A-Za-z]*', '', x[0]))))
    return dNodes

if __name__ == '__main__':
    d = getNodesStatus()
    s = '{:s}\t{:s}\t{:s}\t{:s}'
    print(s.format('Node', 'FreeCpus', 'AvailMem', 'Queue'))
    for i in d:
        print(s.format(i, str(d[i][0]), d[i][1], d[i][2]))


