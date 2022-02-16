#!/home/bin/python3

import os
from time import time
from subprocess import Popen

def main(pbrun, samples, wd):
    for indv, r1, r2 in samples:
        cmd = '{:s} fq2bam --ref {:s}/ref/ref.fa --in-fq {:s} {:s} "@RG\\tID:{:s}\\tLB:{:s}\\tPL:Illumina\\tSM:{:s}\\tPU:unit1" --out-bam {:s}/mapping/{:s}.bam --tmp-dir {:s}/tmp'.format(pbrun, wd, r1, r2, indv, indv, indv, wd, indv, wd)
        with open('{:s}/logs/mapping.{:s}.log'.format(wd, indv), 'w') as f:
            print(cmd, file=f)
            t1 = time()
            p = Popen(cmd, stdout=f, stderr=f, shell=True)
            p.wait()
            t = round(time() - t1, 2)
            print('Elapsed Time: {:f}s.'.format(t), file=f)
