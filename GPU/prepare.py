#!/home/bin/python3

import os
from subprocess import Popen

def makeDirs(wd, ref, samples, mapper):
    for i in ['ref', 'reads', 'mapping', 'logs', 'gvcfs', 'vcfs', 'tmp']:
        try:
            os.mkdir('{:s}/{:s}'.format(wd, i))
        except FileExistsError:
            pass
    refCopy = '{:s}/ref/ref.fa'.format(wd)
    os.system('cp -f {:s} {:s}'.format(ref, refCopy))
    for indv, r1, r2 in samples:
        lt =  r1.split('.')
        suf = lt.pop()
        if suf in ['gz', 'bgzip', 'zip', 'gzip', 'bz2', 'bz', 'xz']:
            if mapper in ['Arioc', 'arioc']:
                raise Exception('Arioc can not read compressed reads.', r1)
            suf = suf + '.' + lt.pop()
#        link1 = '{:s}/reads/{:s}_1.{:s}'.format(wd, indv, suf)
#        link2 = '{:s}/reads/{:s}_2.{:s}'.format(wd, indv, suf)
#        refreshLink(r1, link1)
#        refreshLink(r2, link2)

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

def makeRefIndex(wd, mappingTool, varcallingTool, mapper, varcaller):
    cmd = ''
    if mappingTool == 'bwa':
        cmd += bwaRefIndex(mapper, wd)
    if varcallingTool == 'gatk':
        cmd += gatkRefIndex(varcaller, wd)
    cmd = cmd.removesuffix('&&')
    with open('{:s}/logs/refIndex.log'.format(wd), 'w') as log:
        print(cmd, file=log)
        p = Popen(cmd, stdout=log, stderr=log, shell=True, executable='/bin/bash')
        p.wait()

def main(wd, ref, samples, mappingTool, varcallingTool, mapper, varcaller):
    makeDirs(wd, ref, samples, mapper)
#    makeRefIndex(wd, mappingTool, varcallingTool, mapper, varcaller)
