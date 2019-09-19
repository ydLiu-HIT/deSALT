#! /usr/bin/python

import sys, os
import Annotation_formats
# To enable importing from samscripts submodulew
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
SCRIPT_PATH = os.path.abspath(os.path.join(SCRIPT_PATH, ".."))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))

from fastqparser import read_fastq


class Report:
    def __init__(self):
        self.Total_expected_exons = 0
        self.Total_reads = 0
        self.Total_bases = 0
        self.Total_level2_reads = 0
        self.Total_level2_r_reads = 0
        self.Total_level3_SS_reads = 0
        self.Total_level3_AS_reads = 0
        self.Total_level4_2_5_reads = 0
        self.Total_level4_6_9_reads = 0
        self.Total_level4_10_reads = 0

        self.Total_level2_bases = 0
        self.Total_level2_r_bases = 0
        self.Total_level3_SS_bases = 0
        self.Total_level3_AS_bases = 0
        self.Total_level4_2_5_bases = 0
        self.Total_level4_6_9_bases = 0
        self.Total_level4_10_bases = 0

        self.Total_level2_expected_exons = 0
        self.Total_level2_r_expected_exons = 0
        self.Total_level3_SS_expected_exons = 0
        self.Total_level3_AS_expected_exons = 0
        self.Total_level4_2_5_expected_exons = 0
        self.Total_level4_6_9_expected_exons = 0
        self.Total_level4_10_expected_exons = 0


def processData(read_fastq, annotationfile, ss_list):
    report = Report()
    #load annotation:
    annotations = Annotation_formats.Load_Annotation_From_File(annotationfile)
    annotation_dict = {}
    for annotation in annotations:
        if annotation.transcriptname in annotation_dict:
            pass
        else:
            annotation_dict[annotation.transcriptname] = annotation

    SS_list = list()
    with open(ss_list, 'r') as f:
        for s in f:
            SS_list.append(s.strip())

    fread = open(read_fastq, 'r')

    unaligned = False
    for line in fread:
        if line.startswith('>'):
            pos = line.split('_')
            simGeneName = pos[0].strip('>')
            if "transcript" in simGeneName:
                simGeneName = simGeneName.split(':')[1]

            maf_startpos = int(pos[1])
            aln_sig = pos[2]
            read_idx = int(pos[3])
            maf_strand = pos[4]
            l_clip = int(pos[5])
            maf_length = int(pos[6])
            r_clip = int(pos[7])

            if aln_sig == "unaligned":
                unaligned = True
                continue


            annotation = annotation_dict[simGeneName]       # Getting the correct annotation

            maf_reflen = 0
            for i in range(len(annotation.items)):
                maf_reflen += annotation.items[i].getLength()  # get the reference length from exons itemso

            # IMPORTANT: If the reads were generated from an annotation on reverse strand
            #            expected partial alignments must be reversed
            if annotation.strand == Annotation_formats.GFF_STRANDRV:
                maf_startpos = maf_reflen - maf_length - maf_startpos
            
            # Calculating expected partial alignmetns from MAF and annotations
            sigA = False
            sigB = True
            sigC = False
            sigD = False
            sigE = False
            sigF = False
            sigG = False

            # 1. Calculating the index of the first exon
            # i - the index of exon currently being considered
            i = 0
            while annotation.items[i].getLength() <= maf_startpos:
                maf_startpos -= annotation.items[i].getLength()
                i += 1

            # Calculating expected partial alignments by filling up exons using maf_length
            expected_partial_alignments = []
            while maf_length > 0:
                start = annotation.items[i].start + maf_startpos
                end = annotation.items[i].end
                assert start <= end

                # OLD: length = end-start+1
                # KK: End is already indicating position after the last base, so adding one when callculating length is not correct
                length = end - start
                if length <= maf_length:
                    expected_partial_alignments.append((start, end))
                    maf_length -= length
                    i += 1
                else:
                    expected_partial_alignments.append((start, start + maf_length))
                    maf_length = 0
                    i += 1

                # Start position should only be considered for the first exon
                maf_startpos = 0
            #*****************************************
            #*****************************************

            report.Total_reads += 1 
            num = len(expected_partial_alignments)
            report.Total_expected_exons += len(expected_partial_alignments)
            #level2
            for ele in expected_partial_alignments[1:-1]:
                if ele[1] - ele[0] < 30:
                    report.Total_level2_reads += 1
                    report.Total_level2_expected_exons += num
                    sigA = True
                    break
            if sigA == False:
                report.Total_level2_r_reads += 1
                report.Total_level2_r_expected_exons += num
            #level4
            if num < 6:
                report.Total_level4_2_5_reads += 1
                report.Total_level4_2_5_expected_exons += num
                sigE = True
            elif num > 5 and num < 10:
                report.Total_level4_6_9_reads += 1
                report.Total_level4_6_9_expected_exons += num
                sigF = True
            else:
                report.Total_level4_10_reads += 1
                report.Total_level4_10_expected_exons += num

            #level3
            #print simGeneName
            if simGeneName in SS_list:
                report.Total_level3_SS_reads += 1
                report.Total_level3_SS_expected_exons += num
                sigC = True
            else:
                report.Total_level3_AS_reads += 1
                report.Total_level3_AS_expected_exons += num
                sigD = True
        else:
            if unaligned == True:
                unaligned = False
                continue
            sim_bases = int(len(line))
            report.Total_bases += sim_bases
            #level2
            if sigA == True:
                report.Total_level2_bases += sim_bases
            else:
                report.Total_level2_r_bases += sim_bases

            #level3
            if sigC == True:
                report.Total_level3_SS_bases += sim_bases
            else:
                report.Total_level3_AS_bases += sim_bases
            #level4
            if sigE == True:
                report.Total_level4_2_5_bases += sim_bases
            elif sigF == True:
                report.Total_level4_6_9_bases += sim_bases
            else:
                report.Total_level4_10_bases += sim_bases
 
    fread.close()
    return report

        

if __name__ == '__main__':

    readfile = sys.argv[1]
    annotationfile = sys.argv[2]
    ss_list = sys.argv[3]

    report = processData(readfile, annotationfile, ss_list)
    
