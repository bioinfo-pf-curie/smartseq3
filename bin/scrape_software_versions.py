#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re
from os import path

regexes = {
    'Pipeline': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'seqkit': ['v_seqkit.txt', r"Version: (\S+)"],
    'cutadapt': ['v_cutadapt.txt', r"(\S+)"],
    'umi_tools': ['v_umi_tools.txt', r"UMI-tools version: (\S+)"],
    'STAR': ['v_star.txt', r"(\S+)"],
    'samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'deeptools': ['v_deeptools.txt', r"bamCoverage (\S+)"],
    'rseqc': ['v_rseqc', r"geneBody_coverage.py (\S+)"],
    'R': ['v_R.txt', r"R version (\S+)"],
    'preseq' : ['v_preseq.txt', r"Version: (\S+)"]
}

results = OrderedDict()
results['Pipeline'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['cutadapt'] = '<span style="color:#999999;\">N/A</span>'
results['seqkit'] = '<span style="color:#999999;\">N/A</span>'
results['STAR'] = '<span style="color:#999999;\">N/A</span>'
results['umi_tools'] = '<span style="color:#999999;\">N/A</span>'
results['samtools'] = '<span style="color:#999999;\">N/A</span>'
results['deeptools'] = '<span style="color:#999999;\">N/A</span>'
results['rseqc'] = '<span style="color:#999999;\">N/A</span>'
results['R'] = '<span style="color:#999999;\">N/A</span>'
results['preseq'] = '<span style="color:#999999;\">N/A</span>'


# Search each file using its regex
for k, v in regexes.items():
    if path.isfile(v[0]):
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
<<<<<<< HEAD
                
=======
>>>>>>> 64a61a80d2d39f3f015a1eff7b2e3fe3519f200a
# Dump to YAML
print ('''
id: 'software versions'
section_name: 'Software Versions'
section_href: 'https://gitlab.curie.fr/sc-platform/smartseq3'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in results.items():
    print("        <dt>{}</dt><dd>{}</dd>".format(k,v))
print ("    </dl>")
