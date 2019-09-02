#!/usr/bin/python
"""ATG analysis.
"""

######################################################################Generic header
#import regex
import re
import os
import sys
import collections
from collections import Counter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches
import multiprocessing as mp
from multiprocessing import Pool
#import numpy as np
#import scipy
#from scipy import stats
#from multiprocessing import Pool
#import random

colors = ['#FF0000', '#FF8000', '#FFFF00', '#80FF00', '#00FFBF', '#00BFFF', '#0040FF', '#BF00FF', '#808080', 'k']
colors10 = ['k', '#FF0000', '#FF8000', '#FFFF00', '#80FF00', '#00FFBF', '#00BFFF', '#0040FF', '#BF00FF', '#808080', '#151B54'] 
greys = ['#d3d3d3', '#a8a8a8', '#7e7e7e', '#545454' ]
plt.style.use('ggplot')

def GetAtg(size_lim, atg_in=True):
    """Gets leader sequences as fasta.

    Choose True or False and define a length of leader.
    """
    atg = []

    with open(segmented, 'r') as inF:
        for line in inF:
            if 'ATG' in (line.split('\t'))[3]:
                linea = line.split('\t')
                atg.append(linea[3])

    if atg_in == False:
        outfile = (('no_atg_out_%s.txt') % (size_lim))
        o = open(outfile, 'w')
        i = 0
        for entry in atg:
            if len(entry) >= size_lim:
                if 'ATG' not in entry[0:size_lim]:
                    i = i + 1
                    o.write('>' + str(i) + '\n' + entry[0:size_lim] + '\n')
    else:
        outfile = (('atg_out_%s.txt')%(size_lim))
        o = open(outfile, 'w')
        i = 0
        for entry in atg:
            if len(entry) >= size_lim:
                if 'ATG' in entry[0:size_lim]:
                    i = i + 1
                    o.write('>' + str(i) + '\n' + entry[0:size_lim] + '\n')


def ATGinBam(bampath):
    #
    files = []

    dirs = os.listdir(bampath)
    for infile in dirs:
        totes = 0
        if 'collated' in infile:
            files.append(infile)
    #now you have a list of bamfiles in order with their total lines (reads)
    #print files
    #reorganise
    bamatg = []
    i = 0

    for infile in files:
        i = 0
        with open(infile, 'r') as inF:
            for line in inF:
                
                if '>' not in line:
                    i = i + 1
                    if 'ATG' in line:
                        bamatg.append(line.strip())



    o = open('bam_atg_seqs.txt', 'ab+')
    o.write(('Total bam reads %s\nTotal atgbam reads %s\n')%(str(i), str(len(bamatg))))
    for b in bamatg:
        o.write(str(b) + '\n')
    o.close()
    bamstarts = []
    for b in bamatg:
        i = 0
        while i < len(b) - 3:
            if b[i] == 'A':
                if b[i:i+3] == 'ATG':
                    bamstarts.append(i)
            i = i + 1

    countedbam = Counter(bamstarts)
    countedbam = collections.OrderedDict(sorted(countedbam.items()))

    x = []
    y = []
    y1 = []
    for key in countedbam:
        x.append(int(key))
        y1.append(int(countedbam[key]))
        
        
    sumy = sum(y1)

    for i in y1:
        iy = ((float(i)/float(sumy)))
        y.append(iy)

    print x
    print y

    w = 1
    plt.bar(x,y,w, align = 'center', color = 'r', label = 'ATG start', alpha = 0.5)
    plt.savefig('bam_atg.pdf', format = 'pdf', dpi = 1200)
    plt.show()

    o = open(('bam_atg_seqs.txt'), 'ab+')
    o.write(str(countedbam) + '\n')
    o.close()


def GraphATGStart(atgleaders):
    
    leaders = []

    with open(atgleaders, 'r') as inF:
        for line in inF:
            if '>' in line:
                pass
            else:
                line = line.strip()
                leaders.append(line)



    atgstarts = []
    for b in leaders:
        i = 0
        while i < len(b) - 3:
            if b[i] == 'A':
                if b[i:i+3] == 'ATG':
                    atgstarts.append(i)
            i = i + 1



    countedlens = Counter(atgstarts)

    countedlens = collections.OrderedDict(sorted(countedlens.items()))

    x = []
    y = []
    y1 = []
    for key in countedlens:
        x.append(int(key))
        y1.append(int(countedlens[key]))
        
        
    sumy = sum(y1)

    z = files.index(atgleaders)

    for i in y1:
        iy = ((float(i)/float(sumy)))
        y.append(iy)

    w = 1
    plt.bar(x,y,w, align = 'center', color = colors[z], alpha = 0.5)

    print atgleaders
    print x
    print y

    plt.ylim(ymin = 0, ymax = 0.25)
    plt.xlim(xmin = 0, xmax = 20)
    plt.xticks(range(0, 15, 5), fontsize = 10)

    plt.title('Length of Leader Sequences', fontsize = 10)
    plt.ylabel('Frequency')
    plt.xlabel('Length(nts)')
    plt.savefig(('start_position_%s.pdf')%(str(z)), format='pdf', dpi=1200)
    plt.show()


def SegmentATG(atgleaders):

    outfile = ('ATG_start_0h.txt')
    o = open(outfile, 'w')


    segment_names = ['Segment 8',
                     'Segment 7',
                     'Segment 6',
                     'Segment 5',
                     'Segment 4',
                     'Segment 3',
                     'Segment 2',
                     'Segment 1',
                     'Orphan',
                     'Short']

    f = open(atgleaders, 'r')
    lines = f.readlines()
    f.close()

    len_lines = len(lines)

    for segment in segment_names:

        print segment

        leaders = []

        i = 0
        while i < len(lines):
            if '>' in lines[i]:
                if segment in lines[i]:
                    seq = lines[i + 1].strip()
                    seq = seq.replace('N', '')
                    leaders.append(seq)
            i = i + 1

        print len(leaders)

        atgstarts = []
        for b in leaders:
            print b
            i = 0
            while i < len(b) - 3:
                if b[i] == 'A':
                    if b[i:i + 3] == 'ATG':
                        atgstarts.append(i)
                i = i + 1

        print len(atgstarts)

        countedlens = Counter(atgstarts)
        countedlens = collections.OrderedDict(sorted(countedlens.items()))

        x = []
        y = []
        y1 = []
        for key in countedlens:
            x.append(int(key))
            y1.append(int(countedlens[key]))



        perc = float(len_lines) / 2

        perc2 = sum(y1)

        for i in y1:
            iy = ((float(i) / perc2))
            y.append(iy)

        print y

        if segment == 'Segment 8':
            o.write('Segment\t')
            for entry in x:
                o.write(str(entry) + '\t')
        o.write(segment + '\t')
        for entry in y:
            o.write(str(entry) + '\t')
        o.write('\n')

        w = 1

        n = ('n = %s')%(str(perc2))

        #plt.clf()
        #plt.bar(x,y,w, align = 'center', color = 'k', alpha = 0.5)
        #plt.ylim(ymin = 0, ymax = 0.25)
        #plt.xlim(xmin = 0, xmax = 20)
        #plt.xticks(range(-1, 15, 5), fontsize = 10)
        #plt.text(14, 0.24, n, color = 'k', fontsize=12 )
        #plt.title(('ATG start site for %s at 0h'%(str(segment))), fontsize = 12)
        #plt.ylabel('Frequency')
        #plt.xlabel('Start site')
        #plt.savefig(('start_position_%s_0h_inseg.pdf')%(str(segment)), format='pdf', dpi=1200)
        #plt.show()

'''
ATGinBam(bampath)
files = []
z = 0
dirs = os.listdir(atgleaders)
for infile in dirs:
    if 'collated' in infile:
        files.append(infile)


pool = mp.Pool(processes=4)
results = pool.map(GraphATGStart, files)
'''
atgleaders = "atg_0h_30.fa"
SegmentATG(atgleaders)

######################################
segmented = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/3_segment/publish/segment/20170718_GCAAAGCAGG_segmented.txt'
size_lim = 12
atg_in = False
######################################
GetAtg(size_lim, atg_in=False)
######################################

######################################################################Files required
#segmented = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/3_segment/publish/segment/20170718_GCAAAGCAGG_segmented.txt'
bampath = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/abandoned/4_atg/seqlogo'
#outpath = '/mnt/ris-fas1a/linux_groups2/fantom5/capsnatch/abandoned/4_atg/publish'
######################################################################Functions

