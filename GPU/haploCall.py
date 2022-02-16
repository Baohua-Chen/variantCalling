#!/home/bin/python3

import os
from time import time
from subprocess import Popen

def run(cmd, log):
    with open(log, 'w') as f:
        print(cmd, file=f)
        t1 = time()
        p = Popen(cmd, stdout=f, stderr=f, shell=True)
        p.wait()
        t = round(time() - t1, 2)
        print('Elapsed Time: {:f}s.'.format(t), file=f)

def hcGPU(pbrun, samples, wd):
    for indv, r1, r2 in samples:
        cmd = '{pbrun:s} haplotypecaller --ref {wd:s}/ref/ref.fa --in-bam {wd:s}/mapping/{indv:s}.bam --tmp-dir {wd:s}/tmp --gvcf --out-variants {wd:s}/gvcfs/{indv:s}.g.vcf.gz'.format(pbrun=pbrun, wd=wd, indv=indv)
        log = '{:s}/logs/hc.gpu.{:s}.log'.format(wd, indv)
        run(cmd, log)

def hcCPU(gatk, samples, wd):
    for indv, r1, r2 in samples:
        cmd = '{gatk:s} HaplotypeCaller --java-options "-Xmx32G" -R {wd:s}/ref/ref.fa --emit-ref-confidence GVCF -I {wd:s}/mapping/{indv:s}.bam -O {wd:s}/gvcfs/{indv:s}.g.vcf.gz'.format(gatk=gatk, wd=wd, indv=indv)
        cmd += ' && tabix {wd:s}/gvcfs/{indv:s}.g.vcf.gz'
        log = '{:s}/logs/hc.cpu.{:s}.log'.format(wd, indv)
        run(cmd, log)

def deepvariant(pbrun, samples, wd):
    for indv, r1, r2 in samples:
        cmd = '{pbrun:s} deepvariant --ref {wd:s}/ref/ref.fa --in-bam {wd:s}/mapping/{indv:s}.bam --tmp-dir {wd:s}/tmp --out-variants {wd}/gvcfs/{indv:s}.g.vcf'.format(pbrun=pbrun, wd=wd, indv=indv)
        cmd += ' && bgzip {wd:s}/gvcfs/{indv:s}.g.vcf && tabix {wd}/gvcfs/{indv:s}.g.vcf.gz'.format(wd=wd, indv=indv)
        log = '{:s}/logs/dp.{:s}.log'.format(wd, indv)
        run(cmd, log)
