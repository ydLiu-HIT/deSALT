#! /usr/bin/python
import os
import datetime
import sys

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

    def overlapsItem(self, startpos, endpos, allowed_inacc = DEFAULT_ALLOWED_INACCURACY):
        if (endpos <= self.start + allowed_inacc) or (startpos >= self.end - allowed_inacc):
            return False
        else:
            return True

    def overlapsItemBases(self, startpos, endpos, cigars):
        if startpos < self.start:
            minpos = self.start
        else:
            minpos = startpos

        if endpos < self.end:
            maxpos = endpos + 1
        else:
            maxpos = self.end + 1

        aligned_bases = 0
        toff = startpos
        
        #find starting cigar operator
        pre_op = -1
        start_idx = 0
        if toff < minpos:
            for i in xrange(len(cigars)):
                op1 = cigars[i][0]
                op2 = cigars[i][1]
                start_idx = i
                if toff < minpos:
                    if op1 == 0 or op1 == 2:
                        toff += op2
                else:
                    if pre_op == 0:
                        aligned_bases += (toff - minpos)
                    break
                pre_pos = op1
        for i in xrange(len(cigars)):
            if i >= start_idx:
                op1 = cigars[i][0]
                op2 = cigars[i][1]
                #print toff, maxpos
                if toff <= maxpos:
                    if op1 == 0:
                        toff += op2
                        aligned_bases += op2
                    elif op1 == 1:
                        aligned_bases += op2
                    elif op1 == 2:
                        toff += op2
                else:
                    if pre_op == 0:
                        aligned_bases -= (toff - maxpos)
                    break
                pre_op = op1


        #for op in cigars:
        #    op1 = op[0]
        #    op2 = op[1]
        #    if toff < minpos:
        #        if op1 == 0 or op1 == 2:
        #            toff += op2
        #    elif toff >= minpos and toff <= maxpos:
        #        if op1 == 0:
        #            toff += op2
        #            aligned_bases += op2
        #        elif op1 == 2:
        #            toff += op2
        #    else:
        #        break


        return aligned_bases


    # def equalsItem(self, startpos, endpos, allowed_inacc = DEFAULT_ALLOWED_INACCURACY):
    #     if startpos < self.start - allowed_inacc:
    #         return False
    #     if startpos > self.start + allowed_inacc:
    #         return False
    #     if endpos < self.end - allowed_inacc:
    #         return False
    #     if endpos > self.end + allowed_inacc:
    #         return False

    #     return True

    #exact match exon
    def exactEqualsItem(self, startpos, endpos):
        if startpos == self.start and endpos == self.end:
            return True
        else:
            return False

    #approximate match exon 
    def approxEqualsItem(self, startpos, endpos, allowed_inacc = DEFAULT_ALLOWED_INACCURACY):
        if startpos < self.start - allowed_inacc:
            return False
        if startpos > self.start + allowed_inacc:
            return False
        if endpos < self.end - allowed_inacc:
            return False
        if endpos > self.end + allowed_inacc:
            return False

        return True

    def startsItem(self, startpos, endpos, allowed_inacc = DEFAULT_ALLOWED_INACCURACY):
        if startpos < self.start - allowed_inacc:
            return False
        if startpos > self.start + allowed_inacc:
            return False

        return True

    # Returns true if a given interval (startpos, endpos) and a GeneItem (exon)
    # end at the same position (within allowed_inacc)
    # NOTE: Consider if it might be usefull to also look at the start of the interval
    def endsItem(self, startpos, endpos, allowed_inacc = DEFAULT_ALLOWED_INACCURACY):
        if endpos < self.end - allowed_inacc:
            return False
        if endpos > self.end + allowed_inacc:
            return False

        return True

    # Returns a number of bases by which a given interval and a GeneItem ovelap
    def basesInside(self, startpos, endpos):
        count = 0

        if startpos > self.start:
            maxstart = startpos
        else:
            maxstart = self.start

        if endpos < self.end:
            minend = endpos
        else:
            minend = self.end

        count = minend - maxstart
        if count < 0:
            count = 0

        return count

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

    # this implementation of overlap, will include inside
    def overlapsTranscript(self, startpos, endpos):
        if endpos <= self.start or startpos >= self.end:
            return False
        else:
            return True

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


# class GeneItem(self):
#     self.seqname = ''
#     self.source = ''
#     self.genename = ''
#     self.strand = GFF_STRANDFW
#     self.start = -1
#     self.end = -1
#     self.transcriptitems = []


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

# def Annotation_GeneItem(transcript):
#     gene = GeneItem()
#     gene.seqname = transcript.seqname
#     gene.source = transcript.source
#     gene.start = transcript.start
#     gene.end = transcript.end
#     gene.strand = transcript.strand
#     gene.genename = transcript.genename

#     gene.transcriptitems.append(transcript)

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
    
    
    #cal the average and max exons numbers of one transcript
    max_exon_cnt = 0
    ave_exon_cnt = 0
    for trans in transcripts:
        l = len(trans.exonitems)
        ave_exon_cnt += l
        if max_exon_cnt < l:
            max_exon_cnt = l
    #print "ave = %d, max = %d" %(ave_exon_cnt/len(transcripts), max_exon_cnt)

    return transcripts

def Get_Gene_Annotation(transcripts):
    gene_dict = {}

    gene_name = []

    #method 1, slowly
    # for trans in transcripts:
    #     name = trans.genename
    #     if name not in gene_name:
    #         gene_name.append(name)
    #         gene_dict[name] = []
    #         gene_dict[name].append(trans)
    #     else:
    #         gene_dict[name].append(trans)

    for trans in transcripts:
        gene_name.append(trans.genename)
    gene_name = list(set(gene_name))

    for name in gene_name:
        gene_dict[name] = []

    for trans in transcripts:
        gene_dict[trans.genename].append(trans)

    #print gene and corresponding transcripts
    # for key in gene_dict.keys():
    #     print "gene_name:", key
    #     print "transcriptname:",
    #     for trans in gene_dict[key]:
    #         print trans.transcriptname,
    #     print

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
            for i in xrange(len(att_line)):
                att = att_line[i].split()
                if len(att) == 2:
                    gffline.attribute[att[0]] = att[1].strip('"') #dict for attribute
        gff_lines.append(gffline)

    #get transcript
    transcripts = Get_Transcript_Annotation(gff_lines)
    #get gene
    gene_dict = Get_Gene_Annotation(transcripts)

    return transcripts, gene_dict

def print_lines(gfflines, path):
    f = open(path, 'w')
    for line in gfflines:
        f.write(line.seqname+'\t'+line.source+'\t'+line.feature+'\t'+str(line.start)+'\t'+str(line.end)+'\t'+str(line.score)+'\t'+line.strand+'\t'+str(line.frame)+'\t')
        for key in line.attribute.keys():
            f.write(key+' '+line.attribute[key]+'; ')
        f.write('\n')
       


if __name__ == '__main__':
    # path = '/data/ydliu/human_ONT_2D_simulate/generator/Homo_sapiens.GRCh38.92.chr_gfread.gtf'
    # path1 = '/home/ydliu/out1.gtf'

    path = sys.argv[1]
    #path1 = sys.argv[2]

    starttime = datetime.datetime.now()

    transcripts, gene_dict = Load_Annotation_From_GTF(path)

    for tran in transcripts:
        print tran.transcriptname, tran.genename, len(tran.exonitems)

    #f = open(path1, 'w')

    #for trans in transcripts:
    #    f.write(trans.transcriptname+'\t'+str(trans.start)+'\t'+str(trans.end)+'\n')

    endtime = datetime.datetime.now()

    print "run time: ", (endtime - starttime).seconds
    # print_lines(gfflines, path1)
