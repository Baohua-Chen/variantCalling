#!/home/bin/python3

import os
from time import time
from subprocess import Popen
from multiprocessing import Process

def run(cmd, log):
    with open(log, 'w') as f:
        print(cmd, file=f)
        t1 = time()
        p = Popen(cmd, stdout=f, stderr=f, shell=True)
        p.wait()
        t = round(time() - t1, 2)
        print('Elapsed Time: {:f}s.'.format(t), file=f)

def filter(vcftools, bcftools, tabix, wd, outprefix):
    cmd1 = '{vcftools:s} --gzvcf {wd:s}/vcfs/{outprefix:s}.raw.vcf.gz --max-missing 0.25 --recode --recode-INFO-all --stdout | {bcftools:s} view -v snps --max-alleles 2 --min-af 0.05 -i "QD>2.5 & QUAL>30 & SOR<3 & FS<20 & MQ>40.0 & MQRankSum>-2.5 & ReadPosRankSum<2 & ReadPosRankSum>-2" -O z -o {wd:s}/vcfs/{outprefix:s}.snp.vcf.gz'.format(wd=wd, outprefix=outprefix, vcftools=vcftools, bcftools=bcftools)
    cmd1 += ' && {tabix:s} {wd:s}/vcfs/{outprefix:s}.snp.vcf.gz'.format(wd=wd, outprefix=outprefix, tabix=tabix)
    log1 = '{:s}/logs/snp.filter.log'.format(wd)
    mp1 = Process(target=run, args=(cmd1, log1))
    mp1.start()
    cmd2 = '{vcftools:s} --gzvcf {wd:s}/vcfs/{outprefix:s}.raw.vcf.gz --max-missing 0.25 --recode --recode-INFO-all --stdout | {bcftools:s} view -v indels --min-alleles 2 --min-af 0.05 -i "QD>2 & QUAL>30 & FS<200 & ReadPosRankSum>-20" -O z -o {wd:s}/vcfs/{outprefix:s}.indel.vcf.gz'.format(wd=wd, outprefix=outprefix, vcftools=vcftools, bcftools=bcftools)
    cmd2 += ' && {tabix:s} {wd:s}/vcfs/{outprefix:s}.indel.vcf.gz'.format(wd=wd, outprefix=outprefix, tabix=tabix)
    log2 = '{:s}/logs/indel.filter.log'.format(wd)
    mp2 = Process(target=run, args=(cmd2, log2))
    mp2.start()
    for mp in [mp1, mp2]:
        mp.join()
        mp.close()

