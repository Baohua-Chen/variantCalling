#!/public/home/fgl/public/miniconda3/bin/python3

import pandas as pd
import numpy as np
from subprocess import Popen, PIPE

pre = '/public/home/fgl/chenbh/lc_reseq/gt/'
ref = pre + 'ref/3d.fa'

def grouping(lc):
    lcCopy = lc[:]
    lcCopy = lcCopy.set_index(1).sort_values(2, ascending=False)
    ml = lcCopy[2].sum()/24 # Number of groups is 24. ml: mean length of chromosomes in a group
    lcCopy = list(lcCopy[2].to_dict().items())
    nl,res,v = [],[],[]
    for i in range(24):
        nl.append([lcCopy.pop(0)])
    for i in nl:
        gl = i[0][1]
        i = [i[0][0]]
        while gl < ml and len(lcCopy)>0:
            for j in lcCopy:
                if gl + j[1] <= ml*1.0001:
                    i.append(j[0])
                    lcCopy.remove(j)
                    gl += j[1]
        res.append(i)
        v.append(gl)
        std = np.std(v)
    return res

def readList():
    li = pd.read_csv(pre+'/scripts/list', sep='\t', header=None)[0].to_numpy() # List of individuals
    lc = pd.read_csv(ref.replace('fa','dict'), sep='\t', skiprows=1, header=None)[[1,2]] # List of chromosomes
    lc[1] = lc[1].str.replace('SN:','')
    lc[2] = lc[2].str.replace('LN:','').astype('int')
    gc = grouping(lc) # grouped chromosomes
    lc = lc[1].to_list()
    jids = processPerChrom(li, gc)
#    [ mergeVCFs(lc, i, jids) for i in ['SNP', 'INDEL'] ]
#    [ mergeVCFs(lc, i, None) for i in ['SNP', 'INDEL'] ]

def processPerChrom(li, gc):
    n = 0
#    steps = [ importGenomicsDB, genotype, filterVariants ,variantsToTable]
    steps = [ filterVariants, variantsToTable ]
    qs = ['high']*100
    jids = []
    for i in gc:
        pbs = '#!/bin/bash\n'
        for step in steps:
            pbs = pbs + step(li,i)
        pbs = pbs.strip('&&\\\n')
        fpbs = pre + 'logs/processPerChrom.' + str(n) + '.pbs'
        with open(fpbs, 'w') as f:
            print(pbs, file=f)
        flog = pre + 'logs/processPerChrom.' + str(n) + '.log'
        q = qs.pop(0)
        p = Popen('qsub -N pPerChr_' + str(n) + ' -j oe -o ' + flog + ' -l nodes=1:ppn=3 -q ' + q + ' ' + fpbs, shell=True, stdout=PIPE, executable='/bin/bash')
        p.wait()
        jids.append(p.stdout.readline().decode('utf-8').strip())
        n += 1
    return jids
    
def importGenomicsDB(li, lc):
    cmd = ''
    for c in lc:        
        cmd = cmd + 'echo Start to import ' + c + '. && gatk GenomicsDBImport --java-options "-Djava.io.tmpdir=' + pre + 'java_tmp -Xmx4G"'
        cmd = cmd +  ''.join([ ' -V ' + pre + 'mapping/gvcf/' + i + '.g.vcf.gz' for i in li ])
        cmd = cmd + ' --genomicsdb-workspace-path ' + pre + 'gdb/' + c + ' --intervals ' + c + ' &&\\\n'
    return cmd

def genotype(li, lc):
    cmd = ''
    for c in lc:
        cmd = cmd + 'echo Start to genotype ' + c + '. && gatk GenotypeGVCFs --java-options "-Djava.io.tmpdir=' + pre + 'java_tmp -Xmx4G"'
        cmd = cmd + ' -R ' + ref  + ' -V gendb://' + pre + 'gdb/' + c + ' -G StandardAnnotation -O /dev/stdout | bgzip -@ 3 -c > ' + pre + 'vcfs/premerge/' + c + '.premerge.vcf.gz'
        cmd = cmd + ' && tabix ' + pre + 'vcfs/premerge/' + c + '.premerge.vcf.gz &&\\\n'
    return cmd

def filterVariants(li, lc):
    cmd = ''
    for c in lc:
        cmd = cmd + 'echo Start filter SNPs on ' + c + '. &&\\\n'
#        cmd = cmd + 'gatk --java-options "-Djava.io.tmpdir=' + pre + 'java_tmp -Xmx4G"'
#        cmd = cmd + ' SelectVariants -select-type SNP -R ' + ref + ' -V ' + pre + 'vcfs/premerge/' + c + '.premerge.vcf.gz'
#        cmd = cmd + ' -O ' + pre + 'vcfs/rawSNPs/' + c + '.rawSNP.vcf && \\\n'
        cmd = cmd + 'gatk --java-options "-Djava.io.tmpdir=' + pre + 'java_tmp -Xmx4G"'
        cmd = cmd + ' VariantFiltration -R ' + ref + ' -V ' + pre + 'vcfs/rawSNPs/' + c + '.rawSNP.vcf'
        cmd = cmd + ' -filter "QD < 2.0" --filter-name "QD2"'
        cmd = cmd + ' -filter "QUAL < 30.0" --filter-name "QUAL30"'
        cmd = cmd + ' -filter "SOR > 3.0" --filter-name "SOR3"'
        cmd = cmd + ' -filter "FS > 60.0" --filter-name "FS60"'
        cmd = cmd + ' -filter "MQ < 40.0" --filter-name "MQ40"'
        cmd = cmd + ' -filter "MQRankSum < -12.5" --filter-name "MQRS-12.5"'
        cmd = cmd + ' -filter "ReadPosRankSum < -8.0" --filter-name "RPRS-8"'
        cmd = cmd + ' -O /dev/stdout | bgzip -@ 3 -c > ' + pre + 'vcfs/markedSNPs/' + c + '.markedSNP.vcf.gz &&'
        cmd = cmd + ' tabix ' + pre + 'vcfs/markedSNPs/' + c + '.markedSNP.vcf.gz &&'
        cmd = cmd + ' vcftools --gzvcf ' + pre + 'vcfs/markedSNPs/' + c + '.markedSNP.vcf.gz --maf 0.05 --min-alleles 2 --max-alleles 2 --max-missing 0.25 --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -@ 3 -c > ' + pre + 'vcfs/cleanSNPs/' + c + '.cleanSNP.vcf.gz &&'
        cmd = cmd + ' tabix ' + pre + 'vcfs/cleanSNPs/' + c + '.cleanSNP.vcf.gz &&\\\n'
        cmd = cmd + 'echo Start filter INDELs on ' + c + '. &&\\\n'
#        cmd = cmd + 'gatk --java-options "-Djava.io.tmpdir=' + pre + 'java_tmp -Xmx4G"'
#        cmd = cmd + ' SelectVariants -select-type INDEL -R ' + ref + ' -V ' + pre + 'vcfs/premerge/' + c + '.premerge.vcf.gz'
#        cmd = cmd + ' -O ' + pre + 'vcfs/rawINDELs/' + c + '.rawINDEL.vcf && '
        cmd = cmd + 'gatk --java-options "-Djava.io.tmpdir=' + pre + 'java_tmp -Xmx4G"'
        cmd = cmd + ' VariantFiltration -R ' + ref + ' -V ' + pre + 'vcfs/rawINDELs/' + c + '.rawINDEL.vcf'
        cmd = cmd + ' -filter "QD < 2.0" --filter-name "QD2"'
        cmd = cmd + ' -filter "QUAL < 30.0" --filter-name "QUAL30"'
        cmd = cmd + ' -filter "FS > 200.0" --filter-name "FS200"'
        cmd = cmd + ' -filter "ReadPosRankSum < -20.0" --filter-name "RPRS-20"'
        cmd = cmd + ' -O /dev/stdout | bgzip -@ 3 -c > ' + pre + 'vcfs/markedINDELs/' + c + '.markedINDEL.vcf.gz &&'
        cmd = cmd + ' tabix '+ pre + 'vcfs/markedINDELs/' + c + '.markedINDEL.vcf.gz &&'
        cmd = cmd + ' vcftools --gzvcf ' + pre + 'vcfs/markedINDELs/' + c + '.markedINDEL.vcf.gz --maf 0.05 --min-alleles 2 --max-missing 0.25 --remove-filtered-all --recode --recode-INFO-all --stdout | bgzip -@ 3 -c > ' + pre + 'vcfs/cleanINDELs/' + c + '.cleanINDEL.vcf.gz &&'
        cmd = cmd + ' tabix ' + pre + 'vcfs/cleanINDELs/' + c + '.cleanINDEL.vcf.gz &&\\\n'
    return cmd

def variantsToTable(li, lc):
    cmd = ''
    for c in lc:
        cmd = cmd + 'echo Starts to convert VCF to table on ' + c + '. && gatk --java-options "-Djava.io.tmpdir=' + pre + 'java_tmp -Xmx4G"'
        cmd = cmd + ' VariantsToTable -R ' + ref + ' -V ' + pre + 'vcfs/markedSNPs/' + c + '.markedSNP.vcf'
        cmd = cmd + ' -O ' + pre + 'vcfs/tables/snp/' + c + '.snp.txt'
        cmd = cmd + ' -F CRHOM -F FILTER -F POS -F QUAL -F QD -F FS -F MQ -F SOR -F MQRankSum -F ReadPosRankSum ‐‐show-filtered true &&\\\n'
        cmd = cmd + ' gatk --java-options "-Djava.io.tmpdir=' + pre + 'java_tmp -Xmx4G"'
        cmd = cmd + ' VariantsToTable -R ' + ref + ' -V ' + pre + 'vcfs/markedINDELs/' + c + '.markedINDEL.vcf'
        cmd = cmd + ' -O ' + pre + 'vcfs/tables/indel/' + c + '.indel.txt'
        cmd = cmd + ' -F CRHOM -F FILTER -F POS -F QUAL -F QD -F FS -F MQ -F SOR -F MQRankSum -F ReadPosRankSum ‐‐show-filtered true &&\\\n'
    return cmd

def mergeVCFs(lc, typ, w):
    fpbs = pre + '/logs/merge' + typ + '.pbs'
    pbs = '#!/bin/bash\n'
    pbs = pbs + 'gatk --java-options "-Djava.io.tmpdir=' + pre + 'java_tmp -Xmx4G"'
    pbs = pbs + ' GatherVcfs -I ' + ' -I '.join([ pre + 'vcfs/clean' + typ + 's/' + c + '.clean' + typ + '.vcf.gz' for c in lc ])
    pbs = pbs + ' -O /dev/stdout | bgzip -@ 3 -c > ' + pre + 'vcfs/' + typ.lower() + '.vcf.gz && tabix ' + pre + 'vcfs/' + typ.lower() + '.vcf.gz'
    with open(fpbs, 'w') as f:
        print(pbs, file=f)
    flog = pre + 'logs/merge' + typ + '.log'
    q = 'low'
    sub = 'qsub -N merge' + typ + ' -j oe -o ' + flog
    if w != None:
        sub = sub + ' -W depend=afterok:' + ':'.join(w)
    sub = sub + ' -l nodes=1:ppn=3 -q ' + q + ' ' + fpbs
    p = Popen(sub, shell=True, executable='/bin/bash', stdout=PIPE)
    p.wait()
        
readList()
