#TISU regex

import regex

f = open('ATG_segmented.tsv', 'r')
lines = f.readlines()
f.close

# TISU is SAASATGGCGGC
# regex_basic_pattern = '[GC]AA[GC]ATG'
# test_list = ['EIFGGAAGATGXX', 'EIFGGAAXATGXX', 'EIFGXAAXATGXX']

# errors = {0<=e<=2} # not used
# Allows for up to 2 deletions, insertions or substitutions
# PMID: 28584194 allows for 2 mismatches.
'''
o = 'full_mm_tisu_segmented.tsv'
o = open(o, 'w+')
for line in lines:
    # full
    # obj = regex.search(r'[GC]AA[GC]ATGGCGGC', line)
    # partial (potential)
    # obj = regex.search(r'[GC]AA[GC]ATG', line)
    # full with mismatch
    obj = regex.search(r'[GC]AA[GC]ATGGCGGC{0<=e<=2}', line)    
    if obj is not None:
        o.write(line)
'''

#COUNT TISU SEGMENTS
from collections import Counter
confirmed = 'full_mm_tisu_segmented.tsv'
segments = []

with open(confirmed, 'r') as inF:
    for line in inF:
        linea = line.split('\t')
        if 'ATG' in linea[3]:
            segments.append(linea[2])
print len(segments)
counted = Counter(segments)
print counted
