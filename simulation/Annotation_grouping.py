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

def paras_gene_dict(gene_dict, short_gene_dict, annotation_path):
    single_gene_list = list()
    multip_gene_list = list()
    short_gene_list = list()
    limits = [700, 3000, 300] #AS, SS, short, corresponding gene numbers, can be changed by user

    for key in gene_dict.keys():
        if len(gene_dict[key]) > 1:
            multip_gene_list.append(key)
        elif len(gene_dict[key]) == 1:
            single_gene_list.append(key)

    for key in short_gene_dict.keys():
        short_gene_list.append(key)

    #random
    AS_len = len(multip_gene_list)
    SS_len = len(single_gene_list)
    short_len = len(short_gene_list)
    list_AS = list()
    list_SS = list()
    list_short = list()
    for i in xrange(AS_len):
        rnum = random.randint(0, AS_len)
        if rnum < limits[0]:
            list_AS.append(multip_gene_list[i])

    for i in xrange(SS_len):
        rnum = random.randint(0, SS_len)
        if rnum < limits[1]:
            list_SS.append(single_gene_list[i])

    for i in xrange(short_len):
        rnum = random.randint(0, short_len)
        if rnum < limits[2]:
            list_short.append(short_gene_list[i])
    
    print "number of alternative splicing genes: ", len(list_AS)
    print "number of single splicing genes: ", len(list_SS)
    print "number of genes with short exons: ", len(list_short)


    
    #write coppesponding gene_id to file suffix by "SS.gtf", "AS.gtf", "short.gtf"
    filename, file_ext = os.path.splitext(annotation_path)
    file_AS = filename + '_AS' + file_ext
    file_SS = filename + '_SS' + file_ext
    file_short = filename + "_short" + file_ext
    #file_AS = "AS.gtf"
    #file_SS = "SS.gtf"
    #file_short = "short.gtf"

    with open(file_AS, 'w') as f_AS, open(file_SS, 'w') as f_SS, open(file_short, 'w') as f_short, open(annotation_path, 'r') as afile:
        for line in afile:
            genename = "*"
            elements = line.split('\t')
            att_line = elements[8]
            att_list = att_line.split(';')
            for i in xrange(len(att_list)):
                ele = att_list[i].split()
                if len(ele) > 1 and ele[0] == "gene_id":
                    genename = ele[1][1:-1]
            if genename in list_SS:
                f_SS.write(line)
            if genename in list_AS:
                f_AS.write(line)
            if genename in list_short:
                f_short.write(line)


def Get_Gene_Annotation(transcripts):
    gene_dict = {}
    short_gene_dict = {}

    gene_name = []
    short_gene_name = []
    
    for trans in transcripts:
        sig = 1
        for i in xrange(len(trans.exonitems)):
            if trans.exonitems[i].getLength() < 31:
                short_gene_name.append(trans.genename)
                sig = 0
                break
        if sig == 1:
            gene_name.append(trans.genename)
    gene_name = list(set(gene_name))
    short_gene_name = list(set(short_gene_name))

    for name in gene_name:
        gene_dict[name] = []
    for name in short_gene_name:
        short_gene_dict[name] = []

    for trans in transcripts:
        if trans.genename in short_gene_name:
            short_gene_dict[trans.genename].append(trans)
        else:
            gene_dict[trans.genename].append(trans)

    #print gene and corresponding transcripts
    # for key in gene_dict.keys():
    #     print "gene_name:", key
    #     print "transcriptname:",
    #     for trans in gene_dict[key]:
    #         print trans.transcriptname,
    #     print

    return gene_dict, short_gene_dict

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
    gene_dict, short_gene_dict = Get_Gene_Annotation(transcripts)
    
    #seperate genes into three groups: single splicing, alternative splicing and genes with short exons
    paras_gene_dict(gene_dict, short_gene_dict, filename)


    #print all single splicing isoform, alternative splcing isoforms to file, for the use of evaluation
    All_SS_iso = "All_SS_iso.txt"
    All_AS_iso = "All_AS_iso.txt"
    with open(All_SS_iso, 'w') as f_ss, open(All_AS_iso, 'w') as f_as:
        for key in gene_dict.keys():
            if len(gene_dict[key]) == 1:
                f_ss.write(gene_dict[key][0].transcriptname + '\n')
            else:
                for i in xrange(len(gene_dict[key])):
                    f_as.write(gene_dict[key][i].transcriptname + '\n')
        for key in short_gene_dict.keys():
            if len(short_gene_dict[key]) == 1:
                f_ss.write(short_gene_dict[key][0].transcriptname + '\n')
            else:
                for i in xrange(len(short_gene_dict[key])):
                    f_as.write(short_gene_dict[key][i].transcriptname + '\n')


    gtffile.close()

    return transcripts, gene_dict

def print_lines(gfflines, path):
    f = open(path, 'w')
    for line in gfflines:
        f.write(line.seqname+'\t'+line.source+'\t'+line.feature+'\t'+str(line.start)+'\t'+str(line.end)+'\t'+str(line.score)+'\t'+line.strand+'\t'+str(line.frame)+'\t')
        for key in line.attribute.keys():
            f.write(key+' '+line.attribute[key]+'; ')
        f.write('\n')
       


if __name__ == '__main__':
    path = sys.argv[1] #the route of GTF file

    starttime = datetime.datetime.now()

    transcripts, gene_dict = Load_Annotation_From_GTF(path)

    endtime = datetime.datetime.now()

    print "run time: ", (endtime - starttime).seconds
