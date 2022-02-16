#!/home/bin/python3

import pandas as pd
from subprocess import Popen, PIPE
from multiprocessing import Pool

def mapping(lt):
    indv,pop,r1,r2 = lt
    cmd = f'bwa mem -M -R $"@RG\\tID:{indv}PL:XTen\\tSM:{indv}" -t {ncp} -K 100000000 {wd}/ref/ref.fa {r1} {r2}| samtools view -bS - -F 4 -o {wd}/mapping/raw/{indv}.bam'
    return cmd

def sortBam(lt):
    indv,pop,r1,r2 = lt
    cmd = f'echo Sorting Starts at $(date). && gatk SortSam --java-options "-Djava.io.tmpdir={wd}/java_tmp -Xmx4G" --INPUT {wd}/mapping/raw/{indv}.bam'
    cmd += f' --OUTPUT {wd}/mapping/sort/{indv}.sorted.bam --SORT_ORDER coordinate && echo Sorting Ends at $(date).'
    return cmd

def markDup(lt):
    indv,pop,r1,r2 = lt
    cmd = f'echo Deduplicating Starts at $(date). && gatk MarkDuplicates --java-options "-Djava.io.tmpdir={wd}/java_tmp -Xmx4G" --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000'
    cmd += f' --INPUT {wd}/mapping/sort/{indv}.sorted.bam --OUTPUT {wd}/mapping/dd/{indv}.dd.bam'
    cmd += f' --METRICS_FILE {wd}/mapping/dd/{indv}.dd.bam.metrics && echo Deduplicating Ends at $(date) &&'
    cmd += f'samtools index -@ {ncp} {wd}/mapping/dd/{indv}.dd.bam'
    return cmd

def hapCall(lt):
    indv,pop,r1,r2 = lt
    cmd = f'samtools index -@ {ncp} {wd}/mapping/dd/{indv}.dd.bam && echo Haplotyping Starts at $(date). && time gatk HaplotypeCaller --java-options "-Djava.io.tmpdir={wd}/java_tmp -Xmx10G" -R {wd}/ref/ref.fa --emit-ref-confidence GVCF'
    cmd += f' -I {wd}/mapping/dd/{indv}.dd.bam -O {wd}/mapping/gvcf/{indv}.g.vcf && bgzip -@ {ncp} {wd}/mapping/gvcf/{indv}.g.vcf && echo Haplotyping Ends at $(date)'
    return cmd

def runStep(lt, step):
    indv,pop,r1,r2 = lt
    cmd = step(lt).strip(' &')
    with open(f'{wd}/logs/{step.__name__}/{indv}.{step.__name__}.log', 'w') as f:
        p = Popen(cmd, shell=True, executable='/bin/bash', stdout=f, stderr=f)
        p.wait()

def process(lt):
    steps = [hapCall]
    for step in steps:
        runStep(lt, step)

def main(lst):
    d = pd.read_csv(lst, sep=' |\t', header=None, engine='python').values
    pool = Pool(np)
    pool.map(process, d)
    pool.close()
    pool.join()

if __name__ == '__main__':
    wd = 'reseq'
    lst = f'{wd}/list'
    ncp = 6
    np = 14
    main(lst)
