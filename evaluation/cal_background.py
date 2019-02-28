#! /usr/bin/python

import sys, os
import commands
import Annotation_formats

# To enable importing from samscripts submodulew
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
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


def processData(datafolder, annotationfile, ss_list, as_list):

    #load annotation:
    annotations = Annotation_formats.Load_Annotation_From_File(annotationfile)
    annotation_dict = {}
    for annotation in annotations:
        if annotation.transcriptname in annotation_dict:
            pass
        else:
            annotation_dict[annotation.transcriptname] = annotation

    #cal file count
    command = "ls -l %s*.maf | grep '^-'|wc -l" % datafolder
    #print command
    file_count = commands.getstatusoutput(command)[1]
    file_count = int(file_count)
    #print "total %d reads in %s" %(file_count, datafolder)
    AS_list = list()
    SS_list = list()
    with open(as_list, 'r') as f_as, open(ss_list, 'r') as f_ss:
        for line in f_as:
            AS_list.append(line.strip())
        for line in f_ss:
            SS_list.append(line.strip())

    report = Report()

    simFileSuffix = 'sd'

    for i in xrange(file_count):
        simFileName = simFileSuffix + '_%04d' % (i + 1)
        simRefFileName = simFileName + '.ref'
        simMafFileName = simFileName + '.maf'

        simFilePath = datafolder
        simRefFilePath = os.path.join(simFilePath, simRefFileName)
        simMafFilePath = os.path.join(simFilePath, simMafFileName)

        if not os.path.exists(simRefFilePath):
            raise Exception('Reference file for simulated read %s does not exist!' % qname)
        if not os.path.exists(simMafFilePath):
            raise Exception('Sequence alignment (MAF) for simulated read %s does not exist!' % qname)

        # Reading reference file
        [headers, seqs, quals] = read_fastq(simRefFilePath)
        simGeneName = headers[0]
        annotation = annotation_dict[simGeneName]       # Getting the correct annotation

        maf_startpos = maf_length = 0
        i = 0
        l_c = 0
        sigA = False
        sigE = False
        sigF = False
        total_sim_bases = 0
        total_sim_exons = 0
        with open(simMafFilePath, 'rU') as maffile:
            for line in maffile:
                if line[0] == 's':
                    if line.split()[1] == 'ref': # sim ref
                        l_c += 1
                        elements = line.split()
                        maf_startpos = int(elements[2])
                        maf_length = int(elements[3])

                        # Calculating expected partial alignmetns from MAF and annotations

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

                        report.Total_expected_exons += len(expected_partial_alignments)
                        num = len(expected_partial_alignments)
                        total_sim_exons += num
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
                    else: #sim read
                        sim_bases = int(line.split()[3])
                        report.Total_bases += sim_bases
                        total_sim_bases += sim_bases
                        #level2
                        if sigA == True:
                            report.Total_level2_bases += sim_bases
                            sigA = False
                        else:
                            report.Total_level2_r_bases += sim_bases
                        #level4
                        if sigE == True:
                            report.Total_level4_2_5_bases += sim_bases
                            sigE = False
                        elif sigF == True:
                            report.Total_level4_6_9_bases += sim_bases
                            sigF = False
                        else:
                            report.Total_level4_10_bases += sim_bases


            #level3
            #print simGeneName
            if simGeneName in SS_list:
                report.Total_level3_SS_reads += l_c
                report.Total_level3_SS_bases += total_sim_bases
                report.Total_level3_SS_expected_exons += total_sim_exons
            else:
                report.Total_level3_AS_reads += l_c
                report.Total_level3_AS_bases += total_sim_bases
                report.Total_level3_AS_expected_exons += total_sim_exons

    report.Total_reads = report.Total_level3_SS_reads + report.Total_level3_AS_reads
    # print report.Total_reads, report.Total_bases, report.Total_expected_exons
    # print report.Total_level2_reads, report.Total_level2_bases, report.Total_level2_expected_exons
    # print report.Total_level2_r_reads, report.Total_level2_r_bases, report.Total_level2_r_expected_exons
    # print report.Total_level3_AS_reads, report.Total_level3_AS_bases, report.Total_level3_AS_expected_exons
    # print report.Total_level3_SS_reads, report.Total_level3_SS_bases, report.Total_level3_SS_expected_exons
    # print report.Total_level4_2_5_reads, report.Total_level4_2_5_bases, report.Total_level4_2_5_expected_exons
    # print report.Total_level4_6_9_reads, report.Total_level4_6_9_bases, report.Total_level4_6_9_expected_exons
    # print report.Total_level4_10_reads, report.Total_level4_10_bases, report.Total_level4_10_expected_exons
    return report

def merge_static(report1, report2):
    report1.Total_reads += report2.Total_reads
    report1.Total_bases += report2.Total_bases
    report1.Total_expected_exons += report2.Total_expected_exons

    report1.Total_level2_reads += report2.Total_level2_reads
    report1.Total_level2_bases += report2.Total_level2_bases
    report1.Total_level2_expected_exons += report2.Total_level2_expected_exons

    report1.Total_level2_r_reads += report2.Total_level2_r_reads
    report1.Total_level2_r_bases += report2.Total_level2_r_bases
    report1.Total_level2_r_expected_exons += report2.Total_level2_r_expected_exons

    report1.Total_level3_AS_reads += report2.Total_level3_AS_reads
    report1.Total_level3_AS_bases += report2.Total_level3_AS_bases
    report1.Total_level3_AS_expected_exons += report2.Total_level3_AS_expected_exons
    
    report1.Total_level3_SS_reads += report2.Total_level3_SS_reads
    report1.Total_level3_SS_bases += report2.Total_level3_SS_bases
    report1.Total_level3_SS_expected_exons += report2.Total_level3_SS_expected_exons

    report1.Total_level4_2_5_reads += report2.Total_level4_2_5_reads
    report1.Total_level4_2_5_bases += report2.Total_level4_2_5_bases
    report1.Total_level4_2_5_expected_exons += report2.Total_level4_2_5_expected_exons

    report1.Total_level4_6_9_reads += report2.Total_level4_6_9_reads
    report1.Total_level4_6_9_bases += report2.Total_level4_6_9_bases
    report1.Total_level4_6_9_expected_exons += report2.Total_level4_6_9_expected_exons

    report1.Total_level4_10_reads += report2.Total_level4_10_reads
    report1.Total_level4_10_bases += report2.Total_level4_10_bases
    report1.Total_level4_10_expected_exons += report2.Total_level4_10_expected_exons


def process(data_folder, group_list, annotationfile, ss_list, as_list):
    g_list = []
    with open(group_list, 'r') as g:
        lines = g.readlines()
        for line in lines:
            g_list.append(data_folder + line.strip())

    report_total = Report()
    for group_path in g_list:
        report = processData(group_path, annotationfile, ss_list, as_list)
        #merge
        merge_static(report_total, report)

    return report_total


if __name__ == '__main__':
    data_folder = sys.argv[1]
    group_list = sys.argv[2]
    annotationfile = sys.argv[3]
    ss_list = sys.argv[4]
    as_list = sys.argv[5]

    g_list = []
    with open(group_list, 'r') as g:
        lines = g.readlines()
        for line in lines:
            g_list.append(data_folder + group_folder)

    report_total = Report()
    for group_path in g_list:
        report = processData(group_path, annotationfile, ss_list, as_list)
        #merge
        merge_static(report_total, report)
