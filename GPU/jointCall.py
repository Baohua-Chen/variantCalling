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

def glnexus(pbrun, samples, wd, outprefix):
    cmd = '{:s} glnexus --tmp-dir {:s}/tmp --out-bcf {:s}/vcfs/{:s}.bcf '.format(pbrun, wd, wd, outprefix)
    cmd += ' '.join([ '--in-gvcf {:s}/gvcfs/{:s}.g.vcf'.format(wd, i[0]) for i in samples])
    log = '{:s}/logs/glnexus.log'.format(wd)
    run(cmd, log)

def genotypegvcf_discarded(pbrun, samples, wd, outprefix):
#    cmd = '{:s} triocombinegvcf --ref {:s}/ref/ref.fa --out-variants {:s}/vcfs/{:s}.combined.vcf '.format(pbrun, wd, wd, outprefix)
#    cmd += ' '.join([ '--in-gvcf {:s}/gvcfs/{:s}.g.vcf'.format(wd, i[0]) for i in samples])
#    cmd += ' && {:s} genotypegvcf --tmp-dir {:s}/tmp --ref {:s}/ref/ref.fa --in-gvcf {:s}/vcfs/{:s}.combined.vcf --out-vcf {:s}/vcfs/{:s}.raw.vcf.gz '.format(pbrun, wd, wd, wd, outprefix, wd, outprefix)
    cmd = '{:s} genotypegvcf --tmp-dir {:s}/tmp --ref {:s}/ref/ref.fa --out-vcf {:s}/vcfs/{:s}.raw.vcf.gz --keep-tmp '.format(pbrun, wd, wd, wd, outprefix)
    cmd += ' --in-gvcf {:s}/vcfs/{:s}.combined.g.vcf.gz'.format(wd, outprefix)
#    cmd += '--in-selectvariants-dir {:s}/gvcfs'.format(wd)
    log = '{:s}/logs/genotypegvcf.log'.format(wd)
    run(cmd, log)

def genotypegvcf(gatk, samples, wd, outprefix):
    cmd = '{:s} CombineGVCFs --java-options "-Djava.io.tmpdir={:s}/tmp -Xmx32G" -R {:s}/ref/ref.fa -O {:s}/vcfs/{:s}.combined.vcf.gz '.format(gatk, wd, wd, wd, outprfix)
    cmd += ' '.join([ '-V {:s}/gvcfs/{:s}.g.vcf.gz'.format(wd, i[0]) for i in samples ])
    cmd += ' && {:s} GenotypeGVCFs --java-options "-Djava.io.tmpdir={:s}/tmp -Xmx32G" -R {:s}/ref/ref.fa -G StandardAnnotation -V {:s}/vcfs/{:s}.combined.vcf.gz -O {:s}/vcfs/{:s}.raw.vcf.gz '.format(gatk, wd, wd, wd, outprefix, wd, outprefix)
    log = '{:s}/logs/genotypegvcf.log'.format(wd)
    run(cmd, log)

