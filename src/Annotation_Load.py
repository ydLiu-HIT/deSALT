#! /usr/bin/python
import os
import datetime
import sys
import operator

GFF_STRANDFW = '+'
GFF_STRANDRV = '-'
GFF_FRAME = [0, 1, 2]
DEFAULT_ALLOWED_INACCURACY = 5

class GFFLine:
    def __init__(self):
        self.seqname = ''
        self.feature = ''
        self.start = 0
        self.end = 0
        self.strand = GFF_STRANDFW
        self.attribute = {}


def Get_Annotation_by_Gene(gff_lines):
    gene_dict = {}

    gene_name = []

    for line in gff_lines:
        # print line.attribute['gene_id']
        gene_name.append(line.attribute['gene_id'])
    gene_name = list(set(gene_name))

    for name in gene_name:
        gene_dict[name] = []

    for line in gff_lines:
        gene_dict[line.attribute['gene_id']].append(line)


    return gene_dict
    

def parse_exon_by_gene(gene_dict, path_out):

    file_out = open(path_out, "w")
    
    pre_processed_exons = []
    total_exons = 0
    for key in gene_dict.keys():
        exon_by_gene = []
        for line in gene_dict[key]:
            exon_by_gene.append((line.start, line.end, line.seqname, line.strand))
        #merge and reduplication
        exon_by_gene = list(set(exon_by_gene))
        exon_by_gene = sorted(exon_by_gene, key = operator.itemgetter(0,1))

        pre_processed_exons.extend(exon_by_gene)

    pre_processed_exons = list(set(pre_processed_exons))
    pre_processed_exons = sorted(pre_processed_exons, key = operator.itemgetter(0,1))
    
    total_exons = len(pre_processed_exons)
    file_out.write(str(total_exons)+'\n')
    for i in xrange(len(pre_processed_exons)):
        file_out.write(str(pre_processed_exons[i][2]) + '\t' + str(pre_processed_exons[i][3]) + '\t' + str(pre_processed_exons[i][0]) + '\t' + str(pre_processed_exons[i][1]) + '\n')

    file_out.close()


def Load_Annotation_From_GTF(filename, path_out, check_duplicates = True):
    
    fname, ftext = os.path.splitext(filename)
    if ftext not in ['.gtf', '.gff']:
        raise Exception('Invalid annotation file type: %s' %ftext)
    gtffile = open(filename, 'r')

    gene_dict = {}
    gff_lines = []
    strand_dict = {"+":0, "-":1}

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

        if elements[3] == '.':
            gffline.start = 0
        else:
            gffline.start =  int(elements[3]) - 1

        if elements[4] == '.':
            gffline.end = 0
        else:
            gffline.end = int(elements[4]) - 1


        if elements[6] not in [GFF_STRANDFW, GFF_STRANDRV]:
            gffline.strand = 0
        else:
            gffline.strand = strand_dict[elements[6]]

        if elements[8] == '.':
            gffline.attribute = {}
        else:
            att_line = elements[8].split(';')
            for i in xrange(len(att_line)):
                att = att_line[i].split()
                if len(att) == 2:
                    gffline.attribute[att[0]] = att[1].strip('"') #dict for attribute
        gff_lines.append(gffline)

    #get annotation by gene
    gene_dict = Get_Annotation_by_Gene(gff_lines)
    parse_exon_by_gene(gene_dict, path_out)

def print_lines(gfflines, path):
    f = open(path, 'w')
    for line in gfflines:
        f.write(line.seqname+'\t'+line.source+'\t'+line.feature+'\t'+str(line.start)+'\t'+str(line.end)+'\t'+str(line.score)+'\t'+line.strand+'\t'+str(line.frame)+'\t')
        for key in line.attribute.keys():
            f.write(key+' '+line.attribute[key]+'; ')
        f.write('\n')
       


if __name__ == '__main__':

    path = sys.argv[1]
    path_out = sys.argv[2]

    starttime = datetime.datetime.now()

    Load_Annotation_From_GTF(path, path_out)

    #f = open(path1, 'w')

    #for trans in transcripts:
    #    f.write(trans.transcriptname+'\t'+str(trans.start)+'\t'+str(trans.end)+'\n')

    endtime = datetime.datetime.now()

    print "run time: ", (endtime - starttime).seconds
    # print_lines(gfflines, path1)
