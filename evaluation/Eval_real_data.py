#!/usr/bin/env python
# coding=utf-8

import os
import sys
import Annotation_Load
import pysam
from datetime import datetime
from multiprocessing import Pool

class EvalReport:
    def __init__(self):
        self.Num_Match = 0
        self.Num_MisMatch = 0
        self.Num_Insertion = 0
        self.Num_Deletion = 0
        self.Aligned_readlen = 0
        self.Total_readlen = 0
        self.Percentage_bases_aligned = 0.0
        self.Total_lines = 0
        self.Aligned_reads = 0

        self.Num_genes = 0
        self.Num_transcripts = 0
        self.Num_exons = 0
        self.Max_exon_length = 0
        self.Min_exon_length = 0

        self.num_cover_all_exons = 0  #as for transcript
        self.num_cover_some_exons = 0   #as for transcript, cover at least one exon
        self.num_correct_one_exons = 0  #as for transcript, equal at least one exon
        self.num_cover_no_exons = 0
        self.num_exactequal_exons = 0
        self.num_approxequal_exons = 0

        self.num_good_alignment = 0 # as for alignment,  one alignment may covered all exons but may be a good alignment
        self.num_bad_alignment = 0  # 
        self.num_coverall_good_alignment = 0
        self.num_coverall_bad_alignment = 0
        self.num_bad_split_alignments = 0  #the exons number more than the exons in annotation
        self.num_non_spliced_alignments = 0
        self.num_spliced_alignments = 0

        self.num_aligned_transcripts = 0    #transcripts that be covered by reads
        self.num_aligned_exons = 0
        self.num_exonlength_fiften = 0
        self.num_exonlength_twenty = 0
        self.num_exonlength_thirty = 0
        self.num_exonlength_fourty = 0
        self.num_exonlength_fifty = 0
        self.num_exonlength_sixty = 0
        self.num_exonlength_seventy = 0
        self.num_exonlength_eighty = 0 
        self.num_exonlength_ninety = 0 
        self.num_exonlength_hundred = 0
        self.num_exonlength_other = 0

        self.num_exons_samlines = 0

        self.base_level = 0
        self.intron_bases_overlap = 0
        self.intron_bases_more = 0
        self.intron_bases_out = 0

def report_result(report):
    
    print "\n\nStatic Result"
    print """
    (**)Number of alignments (total/with CIGAR) = %d / %d / %.4f 
    (**)Bases in read (total/aligned/percent) = %d / %d / %.2f%%
    Number of matches / mismatchs /inserts / deletes = %d / %d / %d / %d \n
    """ %(report.Total_lines, report.Aligned_reads, report.Aligned_reads / float(report.Total_lines), \
        report.Total_readlen, report.Aligned_readlen, report.Percentage_bases_aligned, \
        report.Num_Match, report.Num_MisMatch, report.Num_Insertion, report.Num_Deletion)


    print """
    Number of Genes / Transcripts / Exons = %d / %d / %d
    (**)Number of samlines cover_All_exons / cover_Some_exons/ equal_one_exon / cover_No_exons = %d / %d / %d / %d
    (**)Number of samlines cover All good / bad alignments = %d / %d
    (**)Number of good / bad alignment = %d / %d
    Number of spliced/non-spliced alignment = %d / %d
    (**)Num of aligned transcripts / exons / exons_in_samlines = %d / %d / %d
    (**)Num of algined exons length less than 15 / 20 / 30 / 40 / 50 / 60 / 70 / 80 / 90 / 100 / other = %d / %d / %d / %d / %d / %d / %d / %d / %d / %d / %d 
    (**)Num of exact-match / approx-match exons = %d / %d
    (**)exon bases / intron bases (overlap, more, out) / (percent) = %d / (%d, %d, %d) / %.4f / %.4f / %.4f / %.4f
    """ %(report.Num_genes, report.Num_transcripts, report.Num_exons, \
        report.num_cover_all_exons, report.num_cover_some_exons, report.num_correct_one_exons, report.num_cover_no_exons, \
        report.num_coverall_good_alignment, report.num_coverall_bad_alignment, \
        report.num_good_alignment, report.num_bad_alignment, \
        report.num_spliced_alignments, report.num_non_spliced_alignments, \
        report.num_aligned_transcripts, report.num_aligned_exons, report.num_exons_samlines, \
        report.num_exonlength_fiften, report.num_exonlength_twenty, report.num_exonlength_thirty, report.num_exonlength_fourty, report.num_exonlength_fifty, report.num_exonlength_sixty, report.num_exonlength_seventy, report.num_exonlength_eighty, report.num_exonlength_ninety, report.num_exonlength_hundred, report.num_exonlength_other, \
        report.num_exactequal_exons, report.num_approxequal_exons, \
        report.base_level, report.intron_bases_overlap, report.intron_bases_more, report.intron_bases_out, report.base_level / float(report.Total_readlen), report.intron_bases_overlap / float(report.Total_readlen), report.intron_bases_more / float(report.Total_readlen), report.intron_bases_out / float(report.Total_readlen))


def Load_And_Process_SAM(filename):
    samfile = pysam.AlignmentFile(filename, 'r')

    return samfile


def isGoodSplitAlignment(exonhitmap, exoncompletemap, exonstartmap, exonendmap):

    isGood = True
    isSpliced = False
    #if not (len(exonhitmap) == len(exoncompletemap) and len(exonhitmap) == len(exonstartmap) and len(exonhitmap) == len(exonendmap)):
        #raise Exception('ERROR: Exon maps have unequal lengths (%d|%d|%d|%d)!' % (len(exonhitmap), len(exoncompletemap), len(exonstartmap), len(exonendmap)))

    for i in exonhitmap.keys():
        if exonhitmap[i] == 0:
            if exoncompletemap[i] <> 0:
                raise Exception('ERROR: HIT map 0 and COMPLETE map nonzero!')
            if exonstartmap[i] <> 0:
                raise Exception('ERROR: HIT map 0 and START map nonzero!')
            if exonendmap[i] <> 0:
                raise Exception('ERROR: HIT map 0 and END map nonzero!')

    # A list of indices of exons for which a hit map is nonzero
    hitlist = [i for i in exonhitmap.keys() if exonhitmap[i] > 0]

    if len(hitlist) == 0:
        return False, False

    starthit = hitlist[0]
    endhit = hitlist[-1]
    middlelist = hitlist[1:-1]

    

    # For an alignment to be spliced, the hit list has to skip some exons in the middle
    for x in hitlist[:-1]:
        if exonhitmap[x+1] - exonhitmap[x] > 1:
            isSpliced = True
            break           # No need to look further

    # For an alignment to be strictly good, it has to be uninterrupted
    # It has to end the first hit exon, start the last exon, end complete all in the middle
    middleOK = True

    for x in xrange(starthit, endhit + 1):
        if exonhitmap[x] == 0:
            middleOK = False
            break
    if middleOK == True:
        for x in middlelist:
            if exoncompletemap[x] == 0:
                middleOK = False
                break           # No need to look further

    if (not middleOK) or exonstartmap[endhit] == 0 or exonendmap[starthit] == 0:
        isGood = False

    return isGood, isSpliced

part_samlines = {}          # A dictionarry containing a list (or deper hierarchy) of samlines for each part
part_annotations = {}       # A dictionarry containing a list of annotations for each part

def eval_mapping(sam_file, annotation_file, thread):

    report = EvalReport()
    #load annotation file
    sys.stderr.write('\n(%s) Loading and processing annotations file ... ' % datetime.now().time().isoformat())
    transcripts, gene_dict = Annotation_Load.Load_Annotation_From_GTF(annotation_file)

    # Sorting annotations according to position
    # NOTE: Might not be necessary because they are generally already sorted in a file
    transcripts.sort(reverse=False, key=lambda annotation: annotation.start)


    Num_genes = len(gene_dict.keys()) # total number of genes
    Num_transcripts = len(transcripts) #total number of transcripts(isoforms), one gene may have more than one isoforms due to alternate splic
    Num_exons = 0
    for trans in transcripts:
        Num_exons += len(trans.exonitems)

    Max_exon_length = 0
    Min_exon_length = 1000000

    #load SAM
    sys.stderr.write('\n(%s) Loading and processing SAM file with mappings ... ' % datetime.now().time().isoformat())
    samlines = Load_And_Process_SAM(sam_file)

    #get chrososome name
    chromname = []
    header = samlines.text.split('\n')
    for lines in header:
        for line in lines.split('\t'):
            if "SN" in line:
                chromname.append(line.split(':')[1])


    #filte multiple mapped read and remain the best alignment
    read_dict = {}
    read_list = []

    for read in samlines:
        read_list.append(read.query_name)
    read_list2 = {}.fromkeys(read_list).keys()
    for qname in read_list2:
        read_dict[qname] = 0
    samlines.close()
    
    samlines = pysam.AlignmentFile(sam_file, 'r')

    #analyzing cigar strings
    Num_Match = 0
    Num_MisMatch = 0
    Num_Insertion = 0
    Num_Deletion = 0
    Aligned_readlen = 0
    Total_readlen = 0
    Percentage_bases_aligned = 0.0
    Total_lines = 0

    Unmapped_lines = 0
    readlength = 0
    samlines_with_CIGAR = []

    for read in samlines:
        based_aligned = 0
        qname = read.query_name
        if read_dict[qname] == 0:
            Total_lines += 1
            read_dict[qname] += 1
            if read.cigartuples == None:
                Total_readlen += len(read.query_sequence)
                Unmapped_lines += 1
                continue

            cigar = []
            for op in read.cigartuples:
                op1 = op[0] #operation
                op2 = op[1] #operation length
                if op1 == 0:
                    Num_Match += op2
                    based_aligned += op2
                elif op1 == 1: #Insertion
                    Num_Insertion += op2
                    based_aligned += op2
                elif op1 == 2 and op2 <= 20:
                    Num_Deletion += op2
                elif op1 == 2 and op2 > 20:
                    op1 = 3
                cigar.append((op1, op2))
        
            Num_MisMatch += read.get_tag("NM")
            Num_Match -= read.get_tag("NM")

            if read.query_sequence == None:
                readlength = read.infer_read_length()
            else:
                readlength = len(read.query_sequence)

            read.cigartuples = cigar
            samlines_with_CIGAR.append(read)

            if (readlength < based_aligned):
                raise Exception("\nERROR counting aligned end total base!")
            Total_readlen += readlength
            Aligned_readlen += based_aligned


    if Total_readlen == 0:
        Percentage_bases_aligned = -1
    else:
        Percentage_bases_aligned = 100 * float(Aligned_readlen) / Total_readlen

    report.Num_Match = Num_Match
    report.Num_MisMatch = Num_MisMatch
    report.Num_Insertion = Num_Insertion
    report.Num_Deletion = Num_Deletion
    report.Aligned_readlen = Aligned_readlen
    report.Total_readlen = Total_readlen
    report.Percentage_bases_aligned = Percentage_bases_aligned
    report.Total_lines = Total_lines
    report.Aligned_reads = Total_lines - Unmapped_lines

    report.Num_genes = Num_genes
    report.Num_transcripts = Num_transcripts
    report.Num_exons = Num_exons
    report.Max_exon_length = Max_exon_length
    report.Min_exon_length = Min_exon_length

    samlines.close()

    #part_samlines = {}          # A dictionarry containing a list (or deper hierarchy) of samlines for each part
    #part_annotations = {}       # A dictionarry containing a list of annotations for each part
    print chromname
    for chro in chromname:
        part_samlines[chro] = []
        part_annotations[chro] = []

    for read in samlines_with_CIGAR:
        refname = read.reference_name
        #print "refname = ", refname, "readname = ", read.query_name
        if refname is not None:
            pass
        else:
            print "refname = ", refname, "readname = ", read.query_name, read.reference_start 
        part_samlines[refname].append(read)

    for trans in transcripts:
        refname = trans.seqname
        if refname in chromname:
            part_annotations[refname].append(trans)
    
    #start multiple processing
    analysis_pool = Pool(processes=thread)
    Final_result = list()

    #for chro in chromname:
    for i in xrange(len(chromname)):
        chro = chromname[i]
        print "process chrosome %s ......., %d" %(chro, len(part_samlines[chro]))

        Final_result.append(analysis_pool.apply_async(eval_mapping_part, args = (chro, )))
    analysis_pool.close()
    analysis_pool.join()
    
    report_dict = dict()
    for i in xrange(len(chromname)):
        chro = chromname[i]
        report_dict[chro] = Final_result[i].get()
    
    merge_result(report, report_dict)

    report_result(report)
    

def merge_result(report, report_dict):
    for value in report_dict.values():
        report.num_cover_all_exons += value.num_cover_all_exons
        report.num_cover_some_exons += value.num_cover_some_exons
        report.num_correct_one_exons += value.num_correct_one_exons
        report.num_cover_no_exons += value.num_cover_no_exons

        report.num_exactequal_exons += value.num_exactequal_exons
        report.num_approxequal_exons += value.num_approxequal_exons
        report.num_good_alignment += value.num_good_alignment
        report.num_bad_alignment += value.num_bad_alignment

        report.num_coverall_good_alignment += value.num_coverall_good_alignment
        report.num_coverall_bad_alignment += value.num_coverall_bad_alignment
        report.num_bad_split_alignments += value.num_bad_split_alignments
        report.num_non_spliced_alignments += value.num_non_spliced_alignments

        report.num_spliced_alignments += value.num_spliced_alignments
        report.num_aligned_transcripts += value.num_aligned_transcripts
        report.num_aligned_exons += value.num_aligned_exons

        report.num_exonlength_fiften += value.num_exonlength_fiften
        report.num_exonlength_twenty += value.num_exonlength_twenty
        report.num_exonlength_thirty += value.num_exonlength_thirty
        report.num_exonlength_fourty += value.num_exonlength_fourty
        report.num_exonlength_fifty += value.num_exonlength_fifty
        report.num_exonlength_sixty += value.num_exonlength_sixty

        report.num_exonlength_seventy += value.num_exonlength_seventy
        report.num_exonlength_eighty += value.num_exonlength_eighty
        report.num_exonlength_ninety += value.num_exonlength_ninety
        report.num_exonlength_hundred += value.num_exonlength_hundred
        report.num_exonlength_other += value.num_exonlength_other
        report.num_exons_samlines += value.num_exons_samlines
        report.base_level += value.base_level
        report.intron_bases_overlap += value.intron_bases_overlap
        report.intron_bases_more += value.intron_bases_more
        report.intron_bases_out += value.intron_bases_out

## eval for each chrososome
def eval_mapping_part(chro):
    samlines = part_samlines[chro]
    transcripts = part_annotations[chro]

    expressed_genes = {}
    report = EvalReport()
    for trans in transcripts:
        expressed_genes[trans.transcriptname] = [0 for i in xrange(len(trans.exonitems) + 1)]

    for read in samlines:
        qname = read.query_name
        badsplit = False
        isGood = False
        isCoverAll = False
        ref_name = read.reference_name
        num_alignments = 1
        aligned_bases = 0
        intron_bases_overlap = 0
        intron_bases_more = 0
        intron_bases_out = 0
        read_len = read.query_length
        aligned_len = read.query_alignment_length
        #print read_len, aligned_len
        # cal the parts in CIGAR
        for op in read.cigartuples:
            if op[0] == 3 or (op[0] == 2 and op[1] > 20):
                num_alignments += 1

        if num_alignments > 1:
            report.num_spliced_alignments += 1
        else:
            report.num_non_spliced_alignments += 1

        if read.flag & 16 == 0:
            readstrand = Annotation_Load.GFF_STRANDFW
        else:
            readstrand = Annotation_Load.GFF_STRANDRV
        
        startpos = read.reference_start
        endpos = read.reference_end
        # Finding candidate annotations
        candidate_annotations = []
        best_match_annotation = None
        for trans in transcripts:
            #check strand
            # if ref_name == trans.seqname and readstrand == trans.strand and trans.overlapsTranscript(startpos, endpos):
            # if ref_name == trans.seqname and trans.overlapsTranscript(startpos, endpos) and num_alignments <= len(trans.exonitems):
            if ref_name == trans.seqname and trans.overlapsTranscript(startpos, endpos):
                candidate_annotations.append(trans)
                #print candidate_annotations


        #find the best match candidate
        #print "read_name = ", read.query_name
        #print "candidate len = ", len(candidate_annotations)
        if len(candidate_annotations) > 0:
            max_score = 0
            for candidate in candidate_annotations:
                #print candidate.seqname, candidate.start, candidate.end, candidate.transcriptname
                score = 0
                rlen = 0
                start = read.reference_start + 1
                for op in read.cigartuples:
                    op1 = op[0]
                    op2 = op[1]
                    if op1 == 0 or (op1 == 2 and op2 <= 20):
                        rlen += op2
                    elif op1 == 3 or (op1 == 2 and op2 > 20):
                        end = start + rlen - 1
                        #print "(%d, %d), " %(start, end)
                        for item in candidate.exonitems:
                            score += item.basesInside(start, end)
                        start = end + op2 + 1
                        rlen = 0

                #the last 
                end = start + rlen - 1
                #print "(%d, %d), " %(start, end)
                for item in candidate.exonitems:
                    score += item.basesInside(start, end)

                #print "name = %s, score = %d, max_score = %d" %(candidate.transcriptname, score, max_score)
                if score > max_score:
                    max_score = score
                    best_match_annotation = candidate

        
        #have find the best annotation
        if best_match_annotation is not None:
            annotation = best_match_annotation

            #if the annotation contain a exon less than 30bp record
            Have = False
            list_a = []
            cnt = 1
            for item in annotation.exonitems:
                l = item.getLength()
                #if l <= 50 and l > 30:
                if l < 30:
                    Have = True
                    list_a.append(cnt)
                cnt += 1

            #**************print best annotations
            # fw.write(annotation.transcriptname+'--'+read.query_name+'--'+str(Have))
            # for i in xrange(len(list_a)):
            #     fw.write('--'+str(list_a[i]))
            # fw.write('\n')
            # for item in annotation.exonitems:
            #     fw.write("(%d,%d), " %(item.start, item.end))
            # fw.write('\n')

            #**************print alignment parts
            #fw.write(read.query_name+':\n')


            if annotation.transcriptname in expressed_genes.keys():
                expressed_genes[annotation.transcriptname][0] += 1
            else:
                expressed_genes[annotation.transcriptname][0] = 1
                
            if num_alignments > len(annotation.exonitems):
                badsplit = True

            #cal the exons in samline
            report.num_exons_samlines += num_alignments
            
            exonhitmap = {(i+1):0 for i in xrange(len(annotation.exonitems))}
            exoncompletemap = {(i+1):0 for i in xrange(len(annotation.exonitems))}
            exonexactmap = {(i+1):0 for i in xrange(len(annotation.exonitems))}
            exonstartmap = {(i+1):0 for i in xrange(len(annotation.exonitems))}
            exonendmap = {(i+1):0 for i in xrange(len(annotation.exonitems))}
            start = read.reference_start + 1 
            qoff = 0
            rlen = 0
            tmp_cigar = []
            for op in read.cigartuples:
                op1 = op[0]
                op2 = op[1]
                tmp_cigar.append((op1, op2))
                if op1 == 0:
                    qoff += op2
                    rlen += op2
                elif op1 == 1:
                    qoff += op2
                elif op1 == 2 and op2 <= 20:
                    rlen += op2
                elif op1 == 3 or (op1 == 2 and op2 > 20):
                    end = start + rlen - 1
                    # fw.write("(%d, %d), " %(start, end))
                    item_idx = 0
                    overlap = 0
                    for item in annotation.exonitems:
                        item_idx += 1
                        if item.overlapsItem(start, end):
                            overlap = 1
                            exonhitmap[item_idx] += 1
                            #cal bases aligned in annotations
                            exon_base = item.overlapsItemBases(start, end, tmp_cigar)
                            aligned_bases += exon_base
                            intron_bases_overlap += (qoff - exon_base)

                            #print "(%d, %d) - (%d, %d)" %(item.start, item.end, start, end)
                            #print tmp_cigar
                            ##if (qoff - exon_base > 0):
                            #print qoff, exon_base, qoff - exon_base
                            if item.approxEqualsItem(start, end):
                                exoncompletemap[item_idx] = 1
                                exonstartmap[item_idx] = 1
                                exonendmap[item_idx] = 1

                                if item.getLength() < 16:
                                    report.num_exonlength_fiften += 1
                                elif item.getLength() < 21:
                                    report.num_exonlength_twenty += 1
                                elif item.getLength() < 31:
                                    report.num_exonlength_thirty += 1
                                elif item.getLength() < 41:
                                    report.num_exonlength_fourty += 1
                                elif item.getLength() < 51:
                                    report.num_exonlength_fifty += 1
                                elif item.getLength() < 61:
                                    report.num_exonlength_sixty += 1
                                elif item.getLength() < 71:
                                    report.num_exonlength_seventy += 1
                                elif item.getLength() < 81:
                                    report.num_exonlength_eighty += 1
                                elif item.getLength() < 91:
                                    report.num_exonlength_ninety += 1
                                elif item.getLength() < 100:
                                    report.num_exonlength_hundred += 1
                                else:
                                    report.num_exonlength_other += 1

                                if item.exactEqualsItem(start, end):
                                    exonexactmap[item_idx] = 1
                            elif item.startsItem(start, end):
                                exonstartmap[item_idx] = 1
                            elif item.endsItem(start, end):
                                exonendmap[item_idx] = 1

                            expressed_genes[annotation.transcriptname][item_idx] += 1 
                            # break
                    if overlap == 0 and (start < annotation.end and end > annotation.start):
                        if (start < annotation.end and end > annotation.start):
                            intron_bases_more += qoff
                        else:
                            intron_bases_out += qoff

                    qoff = 0
                    start = end + op2 + 1
                    rlen = 0
                    tmp_cigar = []

                    

            #the last
            end = start + rlen - 1
            # fw.write("(%d, %d)\n" %(start, end))
            item_idx = 0
            overlap = 0
            for item in annotation.exonitems:
                item_idx += 1
                if item.overlapsItem(start, end):
                    overlap = 1
                    exonhitmap[item_idx] += 1
                    exon_base = item.overlapsItemBases(start, end, tmp_cigar)
                    aligned_bases += exon_base
                    intron_bases_overlap += (qoff - exon_base)

                    #print "(%d, %d) - (%d, %d)" %(item.start, item.end, start, end)
                    #if (qoff - exon_base > 0):
                    #    print qoff - exon_base

                    if item.approxEqualsItem(start, end):
                        exoncompletemap[item_idx] = 1
                        exonstartmap[item_idx] = 1
                        exonendmap[item_idx] = 1

                        if item.getLength() < 16:
                            report.num_exonlength_fiften += 1
                        elif item.getLength() < 21:
                            report.num_exonlength_twenty += 1
                        elif item.getLength() < 31:
                            report.num_exonlength_thirty += 1
                        elif item.getLength() < 41:
                            report.num_exonlength_fourty += 1
                        elif item.getLength() < 51:
                            report.num_exonlength_fifty += 1
                        elif item.getLength() < 61:
                            report.num_exonlength_sixty += 1
                        elif item.getLength() < 71:
                            report.num_exonlength_seventy += 1
                        elif item.getLength() < 81:
                            report.num_exonlength_eighty += 1
                        elif item.getLength() < 91:
                            report.num_exonlength_ninety += 1
                        elif item.getLength() < 100:
                            report.num_exonlength_hundred += 1
                        else:
                            report.num_exonlength_other += 1
                        
                        if item.exactEqualsItem(start, end):
                            exonexactmap[item_idx] = 1
                    elif item.startsItem(start, end):
                        exonstartmap[item_idx] = 1
                    elif item.endsItem(start, end):
                        exonendmap[item_idx] = 1
                    
                    expressed_genes[annotation.transcriptname][item_idx] += 1
                    
            if overlap == 0 and (start < annotation.end and end > annotation.start):
                if (start < annotation.end and end > annotation.start):
                    intron_bases_more += qoff
                else:
                    intron_bases_out += qoff

            num_exons = len(annotation.exonitems)
            num_coverd_exons = len([x for x in exonhitmap.values() if x > 0])


            # if num_coverd_exons == num_exons:
            if num_coverd_exons <= num_alignments:
                isCoverAll = True
                report.num_cover_all_exons += 1  #as for transcript

            report.base_level += aligned_bases
            report.intron_bases_overlap += intron_bases_overlap
            report.intron_bases_more += intron_bases_more
            
            num_correct_one_exons = 0
            if num_coverd_exons > 0:
                report.num_cover_some_exons += 1   #as for transcript
            num_correct_one_exons = len([x for x in exonexactmap.values() if x > 0])
            if num_correct_one_exons > 0:
                report.num_correct_one_exons += 1

            if num_coverd_exons == 0:
                report.num_cover_no_exons += 1

            num_approxequal_exons = len([x for x in exoncompletemap.values() if x > 0])
            num_exactequal_exons = len([x for x in exonexactmap.values() if x > 0])
            

            # fw.write("%d--%d--%d: (" %(num_exons, num_alignments, num_approxequal_exons))
            # for idx in exoncompletemap.keys():
            #     #if exonexactmap[idx] > 0:
            #     if exoncompletemap[idx] > 0:
            #         fw.write(str(idx) + ",")
            # fw.write(")\n--------------------------\n")

            report.num_exactequal_exons += num_exactequal_exons
            report.num_approxequal_exons += num_approxequal_exons
            
            isGood, isSpliced = isGoodSplitAlignment(exonhitmap, exoncompletemap, exonstartmap, exonendmap)
        else: #best_match_annotation is none
            report.num_cover_no_exons += 1   #as for alignment

            for op in read.cigartuples:
                op1 = op[0]
                op2 = op[1]
                if (op1 == 0 or op1 == 1):
                    intron_bases_out += op2
            #intron_bases_out += aligned_len

        report.intron_bases_out += intron_bases_out

        ##good alignment: cover all exons and the bases for two end less than 5bp
        if isCoverAll:
            if isGood:
                report.num_coverall_good_alignment += 1
            else:
                # f_all_not_good.write(read.qname+'\n')
                report.num_coverall_bad_alignment += 1
        if isGood:
            report.num_good_alignment += 1 #all the aligned exon are correct
        else:
            report.num_bad_alignment += 1
            

        if badsplit:
            report.num_bad_split_alignments += 1  #the exons number more than the exons in annotation

    for key in expressed_genes.keys():
        genecnt = expressed_genes[key]
        if genecnt[0] > 0:
            report.num_aligned_transcripts += genecnt[0]
        for cnt in genecnt[1:]:
            if cnt > 0:
                report.num_aligned_exons += cnt

    return report
#    report_result(report)



if __name__ == '__main__':

    sam_path = str(sys.argv[1]) # alignment file 
    annotation_path = str(sys.argv[2]) #annotation file
    thread = int(sys.argv[3]) # thread

    eval_mapping(sam_path, annotation_path, thread)

