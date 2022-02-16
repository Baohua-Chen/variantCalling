#!/home/bin/python3

import os
from subprocess import Popen

def makeDirs(wd, ref):
    for i in ['ref', 'reads', 'mapping', 'logs', 'vcfs', 'tmp']:
        try:
            os.mkdir('{:s}/{:s}'.format(wd, i))
        except FileExistsError:
            pass
    refCopy = '{:s}/ref/ref.fa'.format(wd)
    os.system('cp -f {:s} {:s}'.format(ref, refCopy))


def bwaRefIndex(bwa, wd):
    cmd = '{:s} index {:s}/ref/ref.fa && '.format(bwa, wd)
    return cmd

def bitmapRefIndex(bitmapper, wd):
    cmd = '{bitmapper:s} --index {wd:s}/ref/ref.fa'.format(bitmapper, wd)
    return cmd

def gatkRefIndex(gatk, wd):
    refDict = '{:s}/ref/ref.fa.dict'.format(wd)
    if os.path.isfile(refDict):
        os.remove(refDict)
    cmd = '{:s} CreateSequenceDictionary -R {:s} &&'.format(gatk, refDict.removesuffix('.dict'))
    return cmd

def makeRefIndex(wd, mappingTool, varcallingTool):
    cmd = ''
    if mappingTool == 'bwa':
        cmd += bwaRefIndex('bwa', wd)
    if varcallingTool == 'gatk':
        cmd += gatkRefIndex('gatk', wd)
    cmd = cmd.removesuffix('&&')
    with open('{:s}/logs/refIndex.log'.format(wd), 'w') as log:
        print(cmd, file=log)
        p = Popen(cmd, stdout=log, stderr=log, shell=True, executable='/bin/bash')
        p.wait()

def main(wd, ref, mappingTool, varcallingTool):
    makeDirs(wd, ref)
    makeRefIndex(wd, mappingTool, varcallingTool)

if __name__ == '__main__':
    wd = 'reseq'
    ref = 'ref.fa'
    main(wd, ref, 'bwa', 'gatk')
