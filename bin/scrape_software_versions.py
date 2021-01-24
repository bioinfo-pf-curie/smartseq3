#!/usr/bin/env python
from __future__ import print_function
from collections import OrderedDict
import re

# TODO Add additional regexes for new tools in process get_software_versions
regexes = {
    'Pipeline': ['v_pipeline.txt', r"(\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'seqkit': ['v_seqkit.txt', r"Version: (\S+)"],
	'umi_tools': ['v_umi_tools.txt', r"UMI-tools version: (\S+)"],
	'STAR': ['v_star.txt', r"(\S+)"],
	'samtools': ['v_samtools.txt', r"samtools (\S+)"],
	'deeptools': ['v_deeptools.txt', r"bamCoverage (\S+)"],
    'rseqc': ['v_rseqc', r"geneBody_coverage.py (\S+)"],
	'R': ['v_R.txt', r"R version (\S+)"],
}

results = OrderedDict()
results['Pipeline'] = '<span style="color:#999999;\">N/A</span>'
results['Nextflow'] = '<span style="color:#999999;\">N/A</span>'
results['STAR'] = '<span style="color:#999999;\">N/A</span>'
results['samtools'] = '<span style="color:#999999;\">N/A</span>'
results['deeptools'] = '<span style="color:#999999;\">N/A</span>'
results['rseqc'] = '<span style="color:#999999;\">N/A</span>'
results['R'] = '<span style="color:#999999;\">N/A</span>'


# Search each file using its regex
for k, v in regexes.items():
    with open(v[0]) as x:
        versions = x.read()
        match = re.search(v[1], versions)
        if match:
            results[k] = "v{}".format(match.group(1))

# Dump to YAML
print ('''
id: 'scRNA-SmartSeq3 pipeline software versions'
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
