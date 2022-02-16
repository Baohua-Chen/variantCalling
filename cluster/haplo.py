#!/public/home/fgl/public/miniconda3/bin/python3

import pandas as pd
from subprocess import Popen, PIPE

pre = '/public/home/fgl/chenbh/lc_reseq/gt/'
ref = pre + 'ref/3d.fa'

def readList():
    lt = pd.read_csv(pre+'/scripts/test.list', sep='\t', header=None).to_numpy()
#    queue = ['low']*11 + ['middle']*8 + ['high']*6 + ['defaultApp']*2
    queue = ['high'] * 6
    for i in lt:
        q = queue.pop(0)
        run(i, q)

def bwaMem(lt, nt):
    indv,pop,r1,r2 = lt
    cmd = 'time bwa mem -M -R $"RG\\tID:{:s}PL:XTen\\tSM:{:s} -t {:d} -K 100000000 -o {:s}/mapping/raw/{:s}.{:d}.sam {:s} {:s} {:s}'.format(indv, indv, nt, pre, indv, nt, ref, r1, r2)
    return cmd

def sortBam(lt):
    indv,pop,r1,r2 = lt
    cmd = 'echo Sorting Starts at $(date). && gatk SortSam --java-options "-Djava.io.tmpdir=' + pre + 'java_tmp -Xmx4G" --INPUT '+ pre +'/mapping/raw/' + indv + '.bam'
    cmd = cmd + ' --OUTPUT ' + pre + 'mapping/sorted/' + indv + '.sorted.bam --SORT_ORDER coordinate && echo Sorting Ends at $(date). &&\\'
    return cmd

def markDup(lt):
    indv,pop,r1,r2 = lt
    cmd = 'echo Deduplicating Starts at $(date). && gatk MarkDuplicates --java-options "-Djava.io.tmpdir=' + pre + 'java_tmp -Xmx4G" --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000'
    cmd = cmd + ' --INPUT ' + pre + 'mapping/sorted/' + indv + '.sorted.bam --OUTPUT ' + pre + 'mapping/dd/' + indv + '.dd.bam'
    cmd = cmd + ' --METRICS_FILE ' + pre + 'mapping/dd/' + indv + '.dd.bam.metrics && echo Deduplicating Ends at $(date) &&\\'
    cmd = cmd + 'samtools index -@ 4 ' + pre + 'mapping/dd/' + indv + '.dd.bam &&\\'
    return cmd

def hapCall(lt):
    indv,pop,r1,r2 = lt
    cmd = 'echo Haplotyping Starts at $(date). && time gatk HaplotypeCaller --java-options "-Djava.io.tmpdir=' + pre + 'java_tmp -Xmx10G" -R ' + ref + ' --emit-ref-confidence GVCF'
    cmd = cmd + ' -I ' + pre + 'mapping/raw/' + indv + '.bam -O ' + pre + 'mapping/test_gvcf/' + indv + '.g.vcf '
    cmd = cmd + ' && echo Haplotyping Ends at $(date) &&\\'
    return cmd

def run(lt, q):
    pbs = "#!/bin/bash\n"
    indv,pop,r1,r2 = lt
    steps = [bwaMem, sortBam, markDup, hapCall]
    for i in steps:
        pbs = pbs + '\n' + i(lt)
    pbs = pbs.strip('&&\\')
    fpbs = '{:s}/logs/test.call.{:s}.pbs'.format(pre, indv)
    print(fpbs)
    with open(fpbs, 'w') as f:
        print(pbs, file=f)
    flog = open(fpbs.replace('.pbs', '.log'), 'w')
    cmd = 'qsub -N {:s} -j oe -l nodes=1:ppn=4 -q {:s} {:s}'.format(indv, q, fpbs)
    p = Popen(cmd, shell=True, executable='/bin/bash')
    p.wait()

d = readList()
