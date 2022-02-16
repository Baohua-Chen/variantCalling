#!/home/bin/python3

import os
import re
import sys
from configparser import ConfigParser
sys.path.append(os.path.split(os.path.realpath(__file__))[0])
import prepare
import mapping
import haploCall
import jointCall
import varFilter

class NewConfigParser(ConfigParser):
    def optionxform(self, optionstr):
        return optionstr

def assignToolBin(tool):
    try:
        binary = dict(cfg.items(tool))['bin']
    except KeyError:
        binary = tool
    return binary

cfg = NewConfigParser()
cfg.read('./my.cnf')
cfg.base = dict(cfg.items('base'))

mappingTool = cfg.base['mappingTool']
if mappingTool == 'bwa':
    mapper = assignToolBin('bwa')
elif mappingTool in ['ariocp', 'AriocP']:
    mapper = assignToolBin('AriocP')

varcallingTool = cfg.base['varcallingTool']
if varcallingTool == 'gatk':
    varcaller = assignToolBin('gatk')
elif varcallingTool == 'deepvariant':
    pass

jointcallingTool = cfg.base['jointcallingTool']
vcftools = assignToolBin('vcftools')
bcftools = assignToolBin('bcftools')
tabix = assignToolBin('tabix')

pbrun = cfg.base['pbrun']
wd = cfg.base['workdir']
ref = cfg.base['ref']
sampleList = cfg.base['samplelist']
outprefix = cfg.base['outprefix']


with open(sampleList, 'r') as f:
    samples = [ re.split('[ \t]', line.strip()) for line in f.readlines() ]
prepare.main(wd, ref, samples, mappingTool, varcallingTool, mapper, varcaller)
mapping.main(pbrun, samples, wd)
haploCall.hcCPU(varcaller, samples, wd)
if jointcallingTool == 'glnexus':
    jointCall.glnexus(pbrun, samples, wd, outprefix)
elif jointcallingTool == 'genotypegvcf':
    jointCall.genotypegvcf(varcaller, samples, wd)
else:
    raise Exception('Unrecognized jointcalling tool.', jointcallingTool)
varFilter.filter(vcftools, bcftools, tabix, wd, outprefix)
