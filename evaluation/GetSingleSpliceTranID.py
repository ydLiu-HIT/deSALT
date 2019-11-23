#! /usr/bin/python
import os
import datetime
import sys
import random

GFF_STRANDFW = '+'
GFF_STRANDRV = '-'
GFF_FRAME = [0, 1, 2]
DEFAULT_ALLOWED_INACCURACY = 10

class GFFLine:
    def __init__(self):
        self.seqname = ''
        self.source = ''
        self.feature = ''
        self.start = 0
        self.end = 0
        self.score = 0.0
        self.strand = GFF_STRANDFW
        self.frame = 0
        self.attribute = {}

class ExonItem:
    def __init__(self):
        self.start = 0
        self.end = 0

    def getLength(self):
        return self.end - self.start + 1


class TranscripitItem:
    def __init__(self):
        self.seqname = ''
        self.source = ''
        self.genename = ''
        self.transcriptname = ''
        self.strand = GFF_STRANDFW
        self.start = -1
        self.end = -1
        self.exonitems = []
    
    # recallculate gene start and end position from exons
    def calcBoundsFromItems(self):
        if len(self.exonitems) == 0:
            pass
        else:
            self.start = self.exonitems[0].start
            self.end = self.exonitems[0].end

            start = self.exonitems[-1].start
            end = self.exonitems[-1].end

            if self.start > start:
                self.start = start
            elif self.end < end:
                self.end = end

def Annotation_TranscriptItem(gffline):
    transcript = TranscripitItem()
    transcript.seqname = gffline.seqname
    transcript.source = gffline.source
    transcript.start = gffline.start
    transcript.end = gffline.end
    transcript.strand = gffline.strand

    transcript.transcriptname = gffline.attribute['transcript_id']
    transcript.genename = gffline.attribute['gene_id']

    exonitem = ExonItem()
    exonitem.start = gffline.start
    exonitem.end = gffline.end

    transcript.exonitems.append(exonitem)

    return transcript

def Get_Transcript_Annotation(gff_lines):
    transcripts = []
    old_annt_name = ''
    curr_annt = None
    for gffline in gff_lines:
        new_annt = Annotation_TranscriptItem(gffline)
        new_annt_name = new_annt.transcriptname
        if old_annt_name != new_annt_name:
            if old_annt_name != '':
                curr_annt.calcBoundsFromItems()
                transcripts.append(curr_annt)
            curr_annt = new_annt
        else:
            if new_annt.seqname != curr_annt.seqname or \
                new_annt.source != curr_annt.source or \
                new_annt.strand != curr_annt.strand or \
                new_annt.genename != curr_annt.genename or \
                new_annt.transcriptname != curr_annt.transcriptname:
                raise Exception('Invalid GFF/GTF line for transcript %s' % new_annt_name)
            assert len(new_annt.exonitems) == 1
            curr_annt.exonitems.append(new_annt.exonitems[0])
        
        old_annt_name = new_annt_name

    # Add the last collected annotation
    if old_annt_name != '':
        transcripts.append(curr_annt)
    
    
    return transcripts

def Get_Gene_Annotation(transcripts):
    gene_dict = {}
    gene_name = []
    
    for trans in transcripts:
        gene_name.append(trans.genename)
    gene_name = list(set(gene_name))

    for name in gene_name:
        gene_dict[name] = []

    for trans in transcripts:
        gene_dict[trans.genename].append(trans)
    
    return gene_dict

def Load_Annotation_From_GTF(filename, check_duplicates = True):
    
    fname, ftext = os.path.splitext(filename)
    if ftext not in ['.gtf', '.gff']:
        raise Exception('Invalid annotation file type: %s' %ftext)
    gtffile = open(filename, 'r')

    gene_dict = {}
    transcripts = []
    gff_lines = []

    for line in gtffile:
        elements = line.strip().split('\t')
        gffline = GFFLine()

        if elements[2] == '.':
            gffline.feature = ''
        else:
            gffline.feature = elements[2]

        if gffline.feature != 'exon':
            continue
        
        if elements[0] == '.':
            gffline.seqname = ''
        else:
            gffline.seqname = elements[0]

        if elements[1] == '.':
            gffline.source = ''
        else:
            gffline.source = elements[1]

        if elements[3] == '.':
            gffline.start = 0
        else:
            gffline.start =  int(elements[3])

        if elements[4] == '.':
            gffline.end = 0
        else:
            gffline.end = int(elements[4])

        if elements[5] == '.':
            gffline.score = 0.0
        else:
            gffline.score = float(elements[5])

        if elements[6] not in [GFF_STRANDFW, GFF_STRANDRV]:
            gffline.strand = GFF_STRANDFW
        else:
            gffline.strand = elements[6]

        if elements[7] not in GFF_FRAME:
            gffline.frame = 0
        else:
            gffline.frame = int(elements[7])

        if elements[8] == '.':
            gffline.attribute = {}
        else:
            att_line = elements[8].split(';')
            for i in xrange(len(att_line)):
                att = att_line[i].split()
                if len(att) == 2:
                    gffline.attribute[att[0]] = att[1].strip('"') #dict for attribute
        gff_lines.append(gffline)

    #get transcript
    transcripts = Get_Transcript_Annotation(gff_lines)
    #get gene dict and gene with short exons
    gene_dict = Get_Gene_Annotation(transcripts)

    #print all single splicing isoform, alternative splcing isoforms to file, for the use of evaluation
    All_SS_iso = "All_SS_iso.txt"
    with open(All_SS_iso, 'w') as f_ss:
        for key in gene_dict.keys():
            if len(gene_dict[key]) == 1:
                f_ss.write(gene_dict[key][0].transcriptname + '\n')

    gtffile.close()

    return transcripts, gene_dict


if __name__ == '__main__':
    anno_path = sys.argv[1] #the route of GTF file

    starttime = datetime.datetime.now()

    transcripts, gene_dict = Load_Annotation_From_GTF(anno_path)

    endtime = datetime.datetime.now()

    print "run time: ", (endtime - starttime).seconds
