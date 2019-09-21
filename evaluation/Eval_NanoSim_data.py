#! /usr/bin/python

import sys, os
import csv
import re
import copy
from datetime import datetime

import cal_NanoSim_background

# To enable importing from samscripts submodulew
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__))
SCRIPT_PATH = os.path.abspath(os.path.join(SCRIPT_PATH, ".."))
sys.path.append(os.path.join(SCRIPT_PATH, 'samscripts/src'))
import Annotation_formats
import part_cal
import utility_sam
from fastqparser import read_fastq

# Determines whether to check the strand whene analyzing data
# Due to complications in generating simulated RNA reads, this is False
P_CHECK_STRAND = False

DEBUG = False

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


class AllStatic:
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

def load_and_process_SAM(sam_file, BBMapFormat = False):
    # Loading SAM file into hash
    # Keeping only SAM lines with regular CIGAR string, and sorting them according to position
    qnames_with_multiple_alignments = {}
    [sam_hash, sam_hash_num_lines, sam_hash_num_unique_lines] = utility_sam.HashSAMWithFilter(sam_file, qnames_with_multiple_alignments)

    # If variable BBMapFormat is set to true, all samfiles referring to the same query will be collected together
    # Stil have to decide what to do with query names, currently removing '_part'
    if BBMapFormat:
        new_sam_hash = {}
        for (qname, sam_lines) in sam_hash.iteritems():
            pos = qname.find('_part')
            if pos > -1:
                origQname = qname[:pos]
            else:
                origQname = qname
            if origQname not in new_sam_hash:
                new_sam_hash[origQname] = sam_lines
            else:
                # import pdb
                # pdb.set_trace()
                new_sam_lines = sam_lines + new_sam_hash[origQname]
                new_sam_hash[origQname] = sam_lines

        sam_hash = new_sam_hash

    # NOTE: This is a quick and dirty solution
    # Setting this to true so that large deletions are turned into Ns
    # BBMap marks intron RNA alignment gaps with deletions!
    BBMapFormat = True

    # Reorganizing SAM lines, removing unmapped queries, leaving only the first alignment and
    # other alignments that possibly costitute a split alignment together with the first one
    samlines = []
    cnt = 0
    pattern = '(\d+)(.)'
    # for samline_list in sam_hash.itervalues():
    for (samline_key, samline_list) in sam_hash.iteritems():
        cnt += 1
        if samline_list[0].cigar <> '*' and samline_list[0].cigar <> '':            # if the first alignment doesn't have a regular cigar string, skip

            if BBMapFormat:
                # All deletes that are 10 or more bases are replaced with Ns of the same length
                operations = re.findall(pattern, samline_list[0].cigar)
                newcigar = ''
                for op in operations:
                    op1 = op[1]
                    op0 = op[0]
                    if op[1] == 'D' and int(op[0]) >= 10:
                        op1 = 'N'
                    newcigar += op0 + op1
                samline_list[0].cigar = newcigar


            operations = re.findall(pattern, samline_list[0].cigar)
            split = False

            for op in operations[1:-1]:             # Ns cannot appear as the first or the last operation
                if op[1] == 'N':
                    split = True
                    break
            # If the first alignment is split (had Ns in the middle), keep only the first alignment and drop the others
            if split:
                # Transform split alignments containing Ns into multiple alignments with clipping
                temp_samline_list = []
                posread = 0
                posref = 0      # NOTE: I don't seem to be using this, probably should remove it
                newcigar = ''
                readlength = samline_list[0].CalcReadLengthFromCigar()
                new_samline = copy.deepcopy(samline_list[0])
                mapping_pos = new_samline.pos
                clipped_bases = new_samline.pos - new_samline.clipped_pos
                hclip_seq = 0        # Used with hard clipping, how big part of sequence should be removed
                clip_type = 'S'     # Soft_clipping by default
                for op in operations:
                    if op[1] == 'N' and int(op[0]) > 1:        # Create a new alignment with clipping
                        newcigar += '%dS' % (readlength - posread)      # Always use soft clipping at the end
                        new_samline.cigar = newcigar
                        # After some deliberation, I concluded that this samline doesn't have to have its position changed
                        # The next samline does, and by the size of N operation in cigar string + any operations before
                        temp_samline_list.append(new_samline)
                        new_samline = copy.deepcopy(samline_list[0])
                        mapping_pos += int(op[0])
                        new_samline.pos = mapping_pos
                        new_samline.clipped_pos = new_samline.pos - clipped_bases
                        posref += int(op[0])
                        if clip_type == 'H':
                            new_samline.seq = new_samline.seq[hclip_seq:]
                        newcigar = '%d%c' % (posread, clip_type)
                    else:                   # Expand a current alignment
                        newcigar += op[0] + op[1]
                        if op[1] in ('D', 'N'):
                            posref += int(op[0])
                            mapping_pos += int(op[0])
                        elif op[1] == 'I':
                            posread += int(op[0])
                            # Everything besides deletes and Ns will be clipped in the next partial alignment
                            # Therefore have to adjust both pos and clipped pos
                            clipped_bases += int(op[0])
                            hclip_seq += int(op[0])
                        elif op[1] in ('S', 'H'):
                            clip_type = op[1]
                            # Clipped bases can not appear in the middle of the original cigar string
                            # And they have already been added to the position,
                            # so I shouldn't adjust my mapping_pos and clipped_bases again
                            # TODO: I should probably diferentiate between hars and soft clipping
                            posread += int(op[0])
                            posref += int(op[0])
                        else:
                            posref += int(op[0])
                            posread += int(op[0])
                            clipped_bases += int(op[0])
                            mapping_pos += int(op[0])
                            hclip_seq += int(op[0])

                new_samline.cigar = newcigar
                temp_samline_list.append(new_samline)

                samlines.append(temp_samline_list)
            else:
                temp_samline_list = [samline_list[0]]        # add the first alignment to the temp list
                multi_alignment = False
                for samline in samline_list[1:]:            # look through other alignments and see if they could form a split alignment with the current temp_samline_list
                    if BBMapFormat:
                        # All deletes that are 10 or more bases are replaced with Ns of the same length
                        operations = re.findall(pattern, samline.cigar)
                        newcigar = ''
                        for op in operations:
                            op0 = op[0]
                            op1 = op[1]
                            if op[1] == 'D' and int(op[0]) >= 10:
                                op1 = 'N'
                            newcigar += op0 + op1
                        samline.cigar = newcigar
                    if not join_split_alignment(temp_samline_list, samline):
                        multi_alignment = True

                samlines.append(temp_samline_list)
        else:
            pass

    # Sorting SAM lines according to the position of the first alignment
    samlines.sort(key = lambda samline: samline[0].pos)

    
    #for samline_list in samlines:
    #    print samline_list[0].qname, samline_list[0].rname
    return samlines

def getChromName(header):
    chromname = ''
    longre = r'(chromosome )(\w*)'
    shortre = r'(chr)(\w*)'

    if header.find('mitochondrion') > -1 or header.find('chrM') > -1:
        chromname = 'chrM'
    else:
        match1 = re.search(longre, header)
        match2 = re.search(shortre, header)
        if match1:
            designation = match1.group(2)
            chromname = 'chr%s' % designation
        elif match2:
            designation = match2.group(2)
            chromname = 'chr%s' % designation
        else:
            chromname = 'genome'

    return chromname

def processData(resultfile, annotationfile, SS_list, TotalReport, csv_path):
    sys.stderr.write('\n(%s) Loading and processing SAM file with mappings ... ' % datetime.now().time().isoformat())
    all_sam_lines = load_and_process_SAM(resultfile, BBMapFormat = True)

    # Reading annotation file
    annotations = Annotation_formats.Load_Annotation_From_File(annotationfile)

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
    #allowed_inacc = 5

    # All samlines in a list should have the same query name
    for samline_list in all_sam_lines:
        qname = samline_list[0].qname
        seqlen = len(samline_list[0].seq)

        # Checking the SAM file if all samlines in a list have the same qname
        for samline in samline_list[1:]:
            if samline.qname != qname:
                sys.stderr.write('\nWARNING: two samlines in the same list with different query names (%s/%s)' % (qname, samline.qname))
        
        pos = qname.split('_')
        simGeneName = pos[0]
        maf_startpos = int(pos[1])
        aln_sig = pos[2]
        read_idx = int(pos[3])
        maf_strand = pos[4]
        l_clip = int(pos[5])
        maf_length = int(pos[6])
        r_clip = int(pos[7])

        if "transcript" in simGeneName:
            simGeneName = simGeneName.split(':')[1]
        
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
        
        # Total
        num = len(expected_partial_alignments)
        
        #level2
        for ele in expected_partial_alignments[1:-1]:
            if ele[1] - ele[0] < 30:
                sigA = True
                sigB = False
                break

        #level4
        if num < 6:
            sigE = True
        elif num > 5 and num < 10:
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

        if getChromName(samline_list[0].rname) != getChromName(annotation.seqname):
            static_dict["All"].Total_aligned_reads += 1
            part_cal.cal(static_dict,sigA, sigC, sigE, sigF, "Total_aligned_reads", 1)
        else:
            for samline in samline_list:
                # sl_startpos = samline.pos - 1   # SAM positions are 1-based
                sl_startpos = samline.pos
                reflength = samline.CalcReferenceLengthFromCigar()
                readlength = samline.CalcReadLengthFromCigar()
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
                    #file_correct.write(qname + '\n')
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
            #**************************************************************************************


    #************************************************
    #******************************************write csv
    static_dict["All"].Total_reads = TotalReport.Total_reads + 1
    static_dict["All"].Total_bases = TotalReport.Total_bases + 1
    static_dict["All"].Total_expected_exons = TotalReport.Total_expected_exons + 1
    static_dict["A"].Total_reads = TotalReport.Total_level2_reads + 1
    static_dict["A"].Total_bases = TotalReport.Total_level2_bases + 1
    static_dict["A"].Total_expected_exons = TotalReport.Total_level2_expected_exons + 1
    static_dict["B"].Total_reads = TotalReport.Total_level2_r_reads + 1
    static_dict["B"].Total_bases = TotalReport.Total_level2_r_bases + 1
    static_dict["B"].Total_expected_exons = TotalReport.Total_level2_r_expected_exons + 1
    static_dict["C"].Total_reads = TotalReport.Total_level3_SS_reads + 1
    static_dict["C"].Total_bases = TotalReport.Total_level3_SS_bases + 1
    static_dict["C"].Total_expected_exons = TotalReport.Total_level3_SS_expected_exons + 1
    static_dict["D"].Total_reads = TotalReport.Total_level3_AS_reads + 1
    static_dict["D"].Total_bases = TotalReport.Total_level3_AS_bases + 1
    static_dict["D"].Total_expected_exons = TotalReport.Total_level3_AS_expected_exons + 1
    static_dict["E"].Total_reads = TotalReport.Total_level4_2_5_reads + 1
    static_dict["E"].Total_bases = TotalReport.Total_level4_2_5_bases + 1
    static_dict["E"].Total_expected_exons = TotalReport.Total_level4_2_5_expected_exons + 1
    static_dict["F"].Total_reads = TotalReport.Total_level4_6_9_reads + 1
    static_dict["F"].Total_bases = TotalReport.Total_level4_6_9_bases + 1
    static_dict["F"].Total_expected_exons = TotalReport.Total_level4_6_9_expected_exons + 1
    static_dict["G"].Total_reads = TotalReport.Total_level4_10_reads + 1
    static_dict["G"].Total_bases = TotalReport.Total_level4_10_bases + 1
    static_dict["G"].Total_expected_exons = TotalReport.Total_level4_10_expected_exons + 1

    #print_static_dict(static_dict)

    with open(csv_path, "w") as fw:
        csv_write = csv.writer(fw, dialect = 'excel')
        header = [" ", resultfile]
        csv_write.writerow(header)
        for item in key:
            level = [item, str(static_dict[item].Total_reads) + ' reads/' + str(static_dict[item].Total_bases) + ' bases/' + str(static_dict[item].Total_expected_exons) + ' exons']
            row1 = ["Aligned", static_dict[item].Total_aligned_reads, round(100*static_dict[item].Total_aligned_reads/float(static_dict[item].Total_reads), 2)]
            row2 = ["bases%", static_dict[item].Total_aligned_bases,  round(100*static_dict[item].Total_aligned_bases/float(static_dict[item].Total_bases), 2)]
            #indicator for recall
            line = str(round(100*static_dict[item].ExR100/float(static_dict[item].Total_reads), 2)) + '/' + str(round(100*static_dict[item].ExR90/float(static_dict[item].Total_reads), 2)) + '/' + str(round(100*static_dict[item].ExR80/float(static_dict[item].Total_reads), 2))
            row3 = ["ExR100/90/80%", line]
            #indicator for accuracy
            #line = str(round(100*static_dict[item].ExA100/float(static_dict[item].Total_reads), 2)) + '/' + str(round(100*static_dict[item].ExA90/float(static_dict[item].Total_reads), 2)) + '/' + str(round(100*static_dict[item].ExA80/float(static_dict[item].Total_reads), 2))
            #row4 = ["ExA100/90/80%", line]
            line = str(round(100*static_dict[item].ExA100/float(static_dict[item].Total_reads), 2)) + '/' + str(round(100*static_dict[item].ExA80/float(static_dict[item].Total_reads), 2))
            row4 = ["Read100/80%", static_dict[item].ExA100, static_dict[item].ExA80, line]
            line = str(round(100*static_dict[item].Hit100/float(static_dict[item].Total_reads), 2)) + '/' + str(round(100*static_dict[item].Hit80/float(static_dict[item].Total_reads), 2))
            row5 = ["Hit100/80%", line]
            row6 = ["Exons%", static_dict[item].Total_aligned_exons, round(100*static_dict[item].Total_aligned_exons/float(static_dict[item].Total_expected_exons), 2)]
            csv_write.writerow(level)
            csv_write.writerow(row1)
            csv_write.writerow(row2)
            #csv_write.writerow(row3)
            csv_write.writerow(row4)
            #csv_write.writerow(row5)
            csv_write.writerow(row6)


def verbose_usage_and_exit():
    sys.stderr.write('Simulation study evaluation - A script for evaluation of simulation data generated by NanoSim.\n')
    sys.stderr.write('\n')
    sys.stderr.write('Usage:\n')
    sys.stderr.write('\t%s [readfile] [alignmentfile] [annotationfile] [ss_list] [csv_path]\n\n' %sys.argv[0])
    sys.stderr.write('\t[readfile]: the reads generated by NanoSim\n')
    sys.stderr.write('\t[alignmentfile]: the alignment results (SAM) of simulation data by aligners\n')
    sys.stderr.write('\t[annotationfile]: the annotations (GTF) of reference genome\n')
    sys.stderr.write('\t[ss_list]: the list of transcripts id of all single splicing isoforms\n')
    sys.stderr.write('\t[csv_path]: the results of evaluation\n')

    exit(0)

if __name__ == '__main__':
    if (len(sys.argv) < 2):
        verbose_usage_and_exit()
    
    readfile = sys.argv[1]
    resultfile = sys.argv[2]
    annotationfile = sys.argv[3]
    ss_list = sys.argv[4]
    csv_path = sys.argv[5]

    Array = cal_NanoSim_background.processData(readfile, annotationfile, ss_list)

    processData(resultfile, annotationfile, ss_list, Array, csv_path)
