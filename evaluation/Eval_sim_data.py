#! /usr/bin/python

import sys, os
import paramsparser
import csv

from datetime import datetime

# To enable importing from samscripts submodulew
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))
import Annotation_formats
import cal_background
import part_cal
import RNAseqEval
from report import EvalReport, ReportType
from RNAseq_benchmark import benchmark_params
from fastqparser import read_fastq

# Determines whether to check the strand whene analyzing data
# Due to complications in generating simulated RNA reads, this is False
P_CHECK_STRAND = False

# OLD: Predefined dictionaries for analyzing different datasets
# simFolderDict_d1 = {'SimG1' : 'group1'
#                   , 'SimG2' : 'group2'
#                   , 'SimG3' : 'group3'}
#
# simFolderDict_all = {'SimG1' : 'group1'
#                    , 'SimG2' : 'group2'
#                    , 'SimG3' : 'group3'
#                    , 'SimG1AS' : 'group1_AS'
#                    , 'SimG1SS' : 'group1_SS'
#                    , 'SimG2AS' : 'group2_AS'
#                    , 'SimG2SS' : 'group2_SS'
#                    , 'SimG3AS' : 'group3_AS'
#                    , 'SimG3SS' : 'group3_SS'}


# A dictionary connecting fasta/fastq header prefix with the folder with pbsim generated data
# Containing information for reads with each prefix
# This is used because data is simulated using several pbsim runs to get different
# coverages for different sets of references (in this case transcripts)
# NOTE: this should be changed for different simulations
simFolderDict = benchmark_params.simFolderDict

DEBUG = False
#DEBUG = True

paramdefs = {'--version' : 0,
             '-v' : 0,
             '--split-qnames' : 1,
             '-sqn' : 0,
             '--save_query_names' : 0,
             '--debug' : 0,
             '--print_mapping' : 1}

# Obsolete
def interval_equals(interval1, interval2, allowed_inacc = Annotation_formats.DEFAULT_ALLOWED_INACCURACY):
    if interval1[0] < interval2[0] - allowed_inacc:
        return False
    if interval1[0] > interval2[0] + allowed_inacc:
        return False
    if interval1[1] < interval2[1] - allowed_inacc:
        return False
    if interval1[1] > interval2[1] + allowed_inacc:
        return False

    return True

# Obsolete
def interval_overlaps(interval1, interval2, allowed_inacc = Annotation_formats.DEFAULT_ALLOWED_INACCURACY):

    if (interval1[1] <= interval2[0] + allowed_inacc) or (interval1[0] >= interval2[1] - allowed_inacc):
        return False
    else:
        return True

def basesInside(startpos, endpos, startpos1, endpos1):
        count = 0

        if startpos > startpos1:
            maxstart = startpos
        else:
            maxstart = startpos1

        if endpos < endpos1:
            minend = endpos
        else:
            minend = endpos1

        count = minend - maxstart
        if count < 0:
            count = 0

        return count


class Static:
    def __init__(self):
        self.Total_reads = 0
        self.Total_bases = 0
        self.Total_expected_exons = 0
        self.Total_aligned_reads = 0
        self.Total_aligned_exons = 0
        self.Total_aligned_bases = 0
        self.ExR80 = 0
        self.ExR90 = 0
        self.ExR100 = 0
        self.ExA80 = 0
        self.ExA90 = 0
        self.ExA100 = 0
        self.Hit100 = 0
        self.Hit80 = 0

def processData(datafolder, resultfile, annotationfile, paramdict, Array, SS_list, csv_path):

    split_qnames = False
    filename = ''
    if '--split-qnames' in paramdict:
        split_qnames = True
        filename = paramdict['--split-qnames'][0]

    filename_correct = filename + '_correct.names'
    filename_hitall = filename + '_hitall.names'
    filename_hitone = filename + '_hitone.names'
    filename_bad = filename + '_incorrect.names'
    filename_unmapped = filename + '_unmapped.names'

    printMap = False
    filename_mapping = ''
    if '--print_mapping' in paramdict:
        filename_mapping = paramdict['--print_mapping'][0]
        printMap = True

    file_correct = None
    file_hitall = None
    file_hitone = None
    file_bad = None
    file_unmapped = None
    folder = os.getcwd()

    # If splittng qnames into files, have to open files first
    if split_qnames:
        file_correct = open(os.path.join(folder, filename_correct), 'w+')
        file_hitall = open(os.path.join(folder, filename_hitall), 'w+')
        file_hitone = open(os.path.join(folder, filename_hitone), 'w+')
        file_bad = open(os.path.join(folder, filename_bad), 'w+')

    # Loading results SAM file
    report = EvalReport(ReportType.FASTA_REPORT)    # not really needed, used for unmapped query names
    # Have to preserve the paramdict
    # paramdict = {}

    sys.stderr.write('\n(%s) Loading and processing SAM file with mappings ... ' % datetime.now().time().isoformat())
    all_sam_lines = RNAseqEval.load_and_process_SAM(resultfile, paramdict, report, BBMapFormat = True)


    # Reading annotation file
    annotations = Annotation_formats.Load_Annotation_From_File(annotationfile)

    mapfile = None
    if printMap:
        mapfile = open(filename_mapping, 'w+')

    # Hashing annotations according to name
    annotation_dict = {}
    for annotation in annotations:
        if annotation.transcriptname in annotation_dict:
            pass
            #sys.stderr.write('\nWARNING: anotation with name %s already in the dictionary!' % annotation.genename)
        else:
            #annotation_dict[annotation.genename] = annotation
            annotation_dict[annotation.transcriptname] = annotation

    #***********************************
    #***********************************
    static_dict = {}
    #"A": with exon < 30 "B": exon > 30
    #"C": single splicing "D": alternative splicing
    #"E": 2-5 exons "F": 6-9 exons "G": >10 exons
    key = ["All", "A", "B", "C", "D", "E", "F", "G"]
    for i in xrange(len(key)):
        static_dict[key[i]] = Static()

    ss_array = list()
    with open(SS_list, 'r') as f_ss:
        for line in f_ss:
            ss_array.append(line.strip())
    #**********************************

    allowed_inacc = Annotation_formats.DEFAULT_ALLOWED_INACCURACY       # Allowing some shift in positions
    # Setting allowed inaccuracy
    # allowed_inacc = 25

    # All samlines in a list should have the same query name
    for samline_list in all_sam_lines:
        qname = samline_list[0].qname

        # Checking the SAM file if all samlines in a list have the same qname
        for samline in samline_list[1:]:
            if samline.qname != qname:
                sys.stderr.write('\nWARNING: two samlines in the same list with different query names (%s/%s)' % (qname, samline.qname))

        # Look for the first underscore in query name
        # Everything before that is the simulation folder name
        # Everything after that is simulated query name
        pos = qname.find('_')
        if pos < 0:
            raise Exception('Invalid query name in results file (%s)!' % qname)

        simFolderKey = qname[:pos]
        if simFolderKey not in simFolderDict:
            # import pdb
            # pdb.set_trace()
            raise Exception('Bad simulation folder short name (%s)!' % simFolderKey)
        simFolder = simFolderDict[simFolderKey]
        simQName = qname[pos+1:]

        simFileSuffix = 'sd'

        pos = simQName.find('_')
        pos2 = simQName.find('_part')
        if pos < 0:
            raise Exception('Invalid simulated query name in results file (%s)!' % simQName)

        # BBMap separates a query into smaller parts he can manage
        # Extends query with '_part_#', which has to be ignored
        if pos2 <> -1:
            simQName = simQName[:pos2]

        simRefNumber = int(simQName[1:pos])
        simFileName = simFileSuffix + '_%04d' % simRefNumber
        simRefFileName = simFileName + '.ref'
        simSeqFileName = simFileName + '.fastq'
        simMafFileName = simFileName + '.maf'

        simFilePath = os.path.join(datafolder, simFolder)
        simRefFilePath = os.path.join(simFilePath, simRefFileName)
        # simSeqFilePath = os.path.join(simFilePath, simSeqFileName)
        simMafFilePath = os.path.join(simFilePath, simMafFileName)

        if not os.path.exists(simRefFilePath):
            # import pdb
            # pdb.set_trace()
            raise Exception('Reference file for simulated read %s does not exist!' % qname)
        #if not os.path.exists(simSeqFilePath):
        #    raise Exception('Sequence file for simulated read %s does not exist!' % qname)
        if not os.path.exists(simMafFilePath):
            # import pdb
            # pdb.set_trace()
            raise Exception('Sequence alignment (MAF) for simulated read %s does not exist!' % qname)

        # Reading reference file
        [headers, seqs, quals] = read_fastq(simRefFilePath)
        simGeneName = headers[0]
        annotation = annotation_dict[simGeneName]       # Getting the correct annotation

        #---------------------
        #for i in xrange(len(annotation.items)):
        #    print "(%d,%d)" %(annotation.items[i].start, annotation.items[i].end)

        # Reading MAF file to get original position and length of the simulated read
        # Query name should be a second item
        maf_startpos = maf_length = 0
        maf_reflen = 0
        i = 0
        with open(simMafFilePath, 'rU') as maffile:
            i += 1
            for line in maffile:
                if line[0] == 's':
                    elements = line.split()
                    maf_qname = elements[1]
                    if maf_qname == 'ref':              # Have to remember data for the last reference before the actual read
                        maf_startpos = int(elements[2])
                        maf_length = int(elements[3])
                        maf_strand = elements[4]
                        maf_reflen = int(elements[5])
                    if maf_qname == simQName:
                        # maf_startpos = int(elements[2])
                        # maf_length = int(elements[3])
                        break

        if maf_qname != simQName:
            # import pdb
            # pdb.set_trace()
            raise Exception('ERROR: could not find query %s in maf file %s' % (qname, simMafFileName))

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

        #    print "(%d, %d)" %(start, end)
            
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
        #level2
        for ele in expected_partial_alignments[1:-1]:
            if ele[1] - ele[0] < 30:
                sigA = True
                sigB = False
                break

        #level4
        n = len(expected_partial_alignments)
        if n < 6:
            sigE = True
        elif n > 5 and n < 10:
            sigF = True
        else:
            sigG = True

        #level3
        if simGeneName in ss_array:
            sigC = True
        else:
            sigD = True     

        
        if DEBUG:
            print "exon in expected alignment---------------"
            for i in xrange(len(expected_partial_alignments)):
                print "(%d, %d)" %(expected_partial_alignments[i][0], expected_partial_alignments[i][1])
            print "exon in real alignment-------------"

        numparts = len(expected_partial_alignments)
        # For each part of expected partial alignments, these maps will count
        # how many real partial alignments overlap or equal it
        parteqmap = {(i+1):0 for i in xrange(numparts)}
        parthitmap = {(i+1):0 for i in xrange(numparts)}

        if RNAseqEval.getChromName(samline_list[0].rname) != RNAseqEval.getChromName(annotation.seqname):
            pass
        else:
            for samline in samline_list:
                # sl_startpos = samline.pos - 1   # SAM positions are 1-based
                sl_startpos = samline.pos
                reflength = samline.CalcReferenceLengthFromCigar()
                readlength = samline.CalcAlignedBaseFromCigar()
                #************************
                #************************
                sl_endpos = sl_startpos + reflength

                if DEBUG:
                    print "(%d, %d)" %(sl_startpos, sl_endpos)

                # Comparing a samline to all expected partial alignments
                tmp_aln = 0
                for i in xrange(len(expected_partial_alignments)):
                    expected_alignement = expected_partial_alignments[i]
                    maf_startpos = expected_alignement[0]
                    maf_endpos = expected_alignement[1]

                    if numparts > 2 and i == 0 and abs(sl_endpos - maf_endpos) < allowed_inacc:
                        parteqmap[i+1] += 1
                        parthitmap[i+1] += 1
                    elif numparts > 2 and (i == len(expected_partial_alignments) - 1) and abs(sl_startpos - maf_startpos) < allowed_inacc:
                        parteqmap[i+1] += 1
                        parthitmap[i+1] += 1
                    elif interval_equals((sl_startpos, sl_endpos), (maf_startpos, maf_endpos), allowed_inacc) :
                        parteqmap[i+1] += 1
                        parthitmap[i+1] += 1
                    elif interval_overlaps((sl_startpos, sl_endpos), (maf_startpos, maf_endpos), 5):
                        parthitmap[i+1] += 1

                    
                    if interval_overlaps((sl_startpos, sl_endpos), (maf_startpos, maf_endpos), 5):
                        l = basesInside(sl_startpos, sl_endpos, maf_startpos, maf_endpos)
                        if tmp_aln < l:
                            tmp_aln = l
                if tmp_aln > readlength:
                    tmp_aln = readlength
                static_dict["All"].Total_aligned_bases += tmp_aln
                part_cal.cal(static_dict, sigA, sigC, sigE, sigF, "Total_aligned_bases", tmp_aln)

            #*************************************************************************************
            #*************************************************************************************
            num_recover_exons = len([x for x in parteqmap.values() if x == 1])
            num_hit_exons = len([x for x in parthitmap.values() if x == 1])

            if num_hit_exons == numparts:
                static_dict["All"].Hit100 += 1
                part_cal.cal(static_dict,sigA, sigC, sigE, sigF, "Hit100", 1)
            if num_hit_exons >= int(0.8 * numparts):
                static_dict["All"].Hit80 += 1
                part_cal.cal(static_dict,sigA, sigC, sigE, sigF, "Hit80", 1)

            sam_l = len(samline_list)
            if num_recover_exons == numparts:
                static_dict["All"].ExR100 += 1
                part_cal.cal(static_dict,sigA, sigC, sigE, sigF, "ExR100", 1)
                if num_recover_exons == sam_l:
                    static_dict["All"].ExA100 += 1
                    part_cal.cal(static_dict,sigA, sigC, sigE, sigF, "ExA100", 1)
            if num_recover_exons >= int(0.8 * numparts):
                static_dict["All"].ExR80 += 1
                part_cal.cal(static_dict,sigA, sigC, sigE, sigF, "ExR80", 1)
                if num_recover_exons >= int(0.8 * sam_l):
                    static_dict["All"].ExA80 += 1
                    part_cal.cal(static_dict,sigA, sigC, sigE, sigF, "ExA80", 1)
            if num_recover_exons >= int(0.9 * numparts):
                static_dict["All"].ExR90 += 1
                part_cal.cal(static_dict,sigA, sigC, sigE, sigF, "ExR90", 1)
                if num_recover_exons >= int(0.9 * sam_l):
                    static_dict["All"].ExA90 += 1
                    part_cal.cal(static_dict,sigA, sigC, sigE, sigF, "ExA90", 1)
            static_dict["All"].Total_aligned_exons += num_recover_exons
            part_cal.cal(static_dict,sigA, sigC, sigE, sigF, "Total_aligned_exons", num_recover_exons)
            static_dict["All"].Total_aligned_reads += 1
            part_cal.cal(static_dict,sigA, sigC, sigE, sigF, "Total_aligned_reads", 1)
            #*************************************************************************************


    #************************************************
    #******************************************write csv
    static_dict["All"].Total_reads = Array.Total_reads
    static_dict["All"].Total_bases = Array.Total_bases
    static_dict["All"].Total_expected_exons = Array.Total_expected_exons
    static_dict["A"].Total_reads = Array.Total_level2_reads
    static_dict["A"].Total_bases = Array.Total_level2_bases
    static_dict["A"].Total_expected_exons = Array.Total_level2_expected_exons
    static_dict["B"].Total_reads = Array.Total_level2_r_reads
    static_dict["B"].Total_bases = Array.Total_level2_r_bases
    static_dict["B"].Total_expected_exons = Array.Total_level2_r_expected_exons
    static_dict["C"].Total_reads = Array.Total_level3_SS_reads
    static_dict["C"].Total_bases = Array.Total_level3_SS_bases
    static_dict["C"].Total_expected_exons = Array.Total_level3_SS_expected_exons
    static_dict["D"].Total_reads = Array.Total_level3_AS_reads
    static_dict["D"].Total_bases = Array.Total_level3_AS_bases
    static_dict["D"].Total_expected_exons = Array.Total_level3_AS_expected_exons
    static_dict["E"].Total_reads = Array.Total_level4_2_5_reads
    static_dict["E"].Total_bases = Array.Total_level4_2_5_bases
    static_dict["E"].Total_expected_exons = Array.Total_level4_2_5_expected_exons
    static_dict["F"].Total_reads = Array.Total_level4_6_9_reads
    static_dict["F"].Total_bases = Array.Total_level4_6_9_bases
    static_dict["F"].Total_expected_exons = Array.Total_level4_6_9_expected_exons
    static_dict["G"].Total_reads = Array.Total_level4_10_reads
    static_dict["G"].Total_bases = Array.Total_level4_10_bases
    static_dict["G"].Total_expected_exons = Array.Total_level4_10_expected_exons

    #print_static_dict(static_dict)

    with open(csv_path, "w") as fw:
        csv_write = csv.writer(fw, dialect = 'excel')
        header = [" ", resultfile]
        csv_write.writerow(header)
        for item in key:
            level = [item, str(static_dict[item].Total_reads) + ' reads/' + str(static_dict[item].Total_bases) + ' bases/' + str(static_dict[item].Total_expected_exons) + ' exons']
            row1 = ["Aligned", round(100*static_dict[item].Total_aligned_reads/float(static_dict[item].Total_reads), 2)]
            row2 = ["bases%", round(100*static_dict[item].Total_aligned_bases/float(static_dict[item].Total_bases), 2)]
            line = str(round(100*static_dict[item].ExR100/float(static_dict[item].Total_reads), 2)) + '/' + str(round(100*static_dict[item].ExR90/float(static_dict[item].Total_reads), 2)) + '/' + str(round(100*static_dict[item].ExR80/float(static_dict[item].Total_reads), 2))
            row3 = ["ExR100/90/80%", line]
            line = str(round(100*static_dict[item].ExA100/float(static_dict[item].Total_reads), 2)) + '/' + str(round(100*static_dict[item].ExA90/float(static_dict[item].Total_reads), 2)) + '/' + str(round(100*static_dict[item].ExA80/float(static_dict[item].Total_reads), 2))
            row4 = ["ExA100/90/80%", line]
            line = str(round(100*static_dict[item].Hit100/float(static_dict[item].Total_reads), 2)) + '/' + str(round(100*static_dict[item].Hit80/float(static_dict[item].Total_reads), 2))
            row5 = ["Hit100/80%", line]
            row6 = ["Exons%", round(100*static_dict[item].Total_aligned_exons/float(static_dict[item].Total_expected_exons), 2)]
            csv_write.writerow(level)
            csv_write.writerow(row1)
            csv_write.writerow(row2)
            csv_write.writerow(row3)
            csv_write.writerow(row4)
            #csv_write.writerow(row5)
            csv_write.writerow(row6)


    
def print_static_dict(static_dict):
    for item in static_dict.keys():
        print "static item************************", item
        print "Aligned read = %d, total read = %d" %(static_dict[item].Total_aligned_reads, static_dict[item].Total_reads)
        print "Aligned base = %d, total base = %d" %(static_dict[item].Total_aligned_bases, static_dict[item].Total_bases)
        print "ExR100 = %d, ExR90 = %d, ExR80 = %d" %(static_dict[item].ExR100, static_dict[item].ExR90, static_dict[item].ExR80)
        print "ExA100 = %d, ExA90 = %d, ExA80 = %d" %(static_dict[item].ExA100, static_dict[item].ExA90, static_dict[item].ExA80)
        print "Aligned exon = %d, total exon = %d" %(static_dict[item].Total_aligned_exons, static_dict[item].Total_expected_exons)



def verbose_usage_and_exit():
    sys.stderr.write('Process pbsim data - A tool for processing data generated by pbsim.\n')
    sys.stderr.write('                   - Collects data generated for multiple references.\n')
    sys.stderr.write('                   - And adjusts headers to reflect a reference of origin.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [mode]\n' % sys.argv[0])
    sys.stderr.write('\n')
    sys.stderr.write('\tmode:\n')
    sys.stderr.write('\t\tprocess\n')
    sys.stderr.write('\n')
    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        verbose_usage_and_exit()

    mode = sys.argv[1]

    if (mode == 'process'):
        if (len(sys.argv) < 9):
            sys.stderr.write('Processes a folder containing data generated by pbsim.\n')
            sys.stderr.write('Joins all generated reads into a single FASTQ file.\n')
            sys.stderr.write('Expands existing headers with the name of originating reference.\n')
            sys.stderr.write('Usage:\n')
            sys.stderr.write('%s %s <pbsim data folder> <results file> <annotations file> <options>\n'% (sys.argv[0], sys.argv[1]))
            sys.stderr.write('\n')
            sys.stderr.write('\noptions:\n')
            sys.stderr.write('\t\t--split-qnames: while calculating the statistics also sorts query names\n')
            sys.stderr.write('\t\t                into four files - file_correct.names, file_hitall.names\n')
            sys.stderr.write('\t\t                                  file_hitone.names, file_bad.names\n')
            sys.stderr.write('\t\t--print_mapping [filename]: Print information about actual and expected alignments\n')
            sys.stderr.write('\t\t                into a give text file.\n')
            sys.stderr.write('\n')
            exit(1)

        datafolder = sys.argv[2]
        resultfile = sys.argv[3]
        annotationfile = sys.argv[4]
        group_list = sys.argv[5]
        ss_list = sys.argv[6]
        as_list = sys.argv[7]
        csv_path = sys.argv[8]

        #Total_reads = int(sys.argv[5])

        pparser = paramsparser.Parser(paramdefs)
        paramdict = pparser.parseCmdArgs(sys.argv[9:])
        paramdict['command'] = ' '.join(sys.argv)

        Array = cal_background.process(datafolder, group_list, annotationfile, ss_list, as_list)

        print "Total reads: ", Array.Total_reads
        print "Total bases: ", Array.Total_bases
        print "Total exons:", Array.Total_expected_exons
        #print Array.Total_reads, Array.Total_bases, Array.Total_expected_exons
        #print Array.Total_level2_reads, Array.Total_level2_bases, Array.Total_level2_expected_exons
        #print Array.Total_level2_r_reads, Array.Total_level2_r_bases, Array.Total_level2_r_expected_exons
        #print Array.Total_level3_AS_reads, Array.Total_level3_AS_bases, Array.Total_level3_AS_expected_exons
        #print Array.Total_level3_SS_reads, Array.Total_level3_SS_bases, Array.Total_level3_SS_expected_exons
        #print Array.Total_level4_2_5_reads, Array.Total_level4_2_5_bases, Array.Total_level4_2_5_expected_exons
        #print Array.Total_level4_6_9_reads, Array.Total_level4_6_9_bases, Array.Total_level4_6_9_expected_exons
        #print Array.Total_level4_10_reads, Array.Total_level4_10_bases, Array.Total_level4_10_expected_exons

        processData(datafolder, resultfile, annotationfile, paramdict, Array, ss_list, csv_path)

    else:
        print 'Invalid mode!'
