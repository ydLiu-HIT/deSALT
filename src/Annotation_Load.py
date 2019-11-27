#! /usr/bin/python
import os
import datetime
import sys
import random

GFF_STRANDFW = '+'
GFF_STRANDRV = '-'
GFF_FRAME = [0, 1, 2]

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
        self.genename = ''
        self.source = ''
        self.transcriptname = ''
        self.strand = GFF_STRANDFW
        self.start = -1
        self.end = -1
        self.exonitems = []

    # Recallculate gen start and end position from exons
    def calcBoundsFromItems(self):
        if len(self.exonitems) == 0:
            pass
        else:
            self.start = self.exonitems[0].start
            self.end = self.exonitems[0].end

            start = self.exonitems[-1].start
            end = self.exonitems[-1].end

            if self.start > start:  #decrease order
                self.start = start
            elif self.end < end:   #increase order
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

def outPutTrans2bed(transcripts, fpath):
    TC_exon = 0
    TC_intron = 0
    #out as exon pos
    for trans in transcripts:
        if trans.strand == GFF_STRANDRV:
            trans.exonitems = trans.exonitems[::-1]

        TC_exon += len(trans.exonitems)
        TC_intron += len(trans.exonitems) - 1

    with open(fpath, 'w') as fw:
        fw.write(str(TC_exon) + '\t' + str(TC_intron) + "\n")
        for trans in transcripts:
            header = "%s\t%s|%s\t%s\t%s\t" %(trans.seqname, trans.transcriptname, trans.genename, trans.strand, str(len(trans.exonitems)))
            fw.write(header)
            
            for  I in trans.exonitems:
                fw.write(str(I.start) + "," + str(I.end) + ",")
            fw.write('\n')

def Load_Annotation_From_GTF(filename, fpath, check_duplicates = True):
    
    fname, ftext = os.path.splitext(filename)
    if ftext not in ['.gtf', '.gff']:
        raise Exception('Invalid annotation file type: %s' %ftext)
    gtffile = open(filename, 'r')

    gene_dict = {}
    transcripts = []
    gff_lines = []

    for line in gtffile:
        if line[0] == "#":
            continue
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
            for i in range(len(att_line)):
                att = att_line[i].split()
                if len(att) == 2:
                    gffline.attribute[att[0]] = att[1].strip('"') #dict for attribute
        gff_lines.append(gffline)

    #get transcript
    transcripts = Get_Transcript_Annotation(gff_lines)

    outPutTrans2bed(transcripts, fpath)

    gtffile.close()

    return transcripts


if __name__ == '__main__':
    gffpath = sys.argv[1]
    bedpath = sys.argv[2]

    starttime = datetime.datetime.now()

    transcripts = Load_Annotation_From_GTF(gffpath, bedpath)

    endtime = datetime.datetime.now()
