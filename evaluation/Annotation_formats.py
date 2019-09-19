#! /usr/bin/python

import sys, os

GFF_STRANDFW = '+'
GFF_STRANDRV = '-'
GFF_FRAME = [0, 1, 2]

# While calculating operations on intervals (genes, exons and read mappings)
# some inaccuracy will be allowed. This will have different impact on different operations.
# i.e. For two intervals to overlap, they will have to overlap on at least ALLOWED_INACCURACY bases.

DEFAULT_ALLOWED_INACCURACY = 10 #before 5


class GeneItem:
    def __init__(self):
        self.itemName = ''
        self.start = 0
        self.end = 0
        self.frame = 0

    def getLength(self):
        return self.end - self.start

    def isValidInterval(self):
        if self.start >= self.end:
            return False
        else:
            return True

    # Returns true if a given interval (startpos, endpos) is inside a GeneItem (exon)
    # The interval can extend outside GeneItem at most allowed_inacc bases
    def insideItem(self, startpos, endpos, allowed_inacc = DEFAULT_ALLOWED_INACCURACY):
        if (startpos >= self.start - allowed_inacc) and (endpos <= self.end + allowed_inacc):
            return True
        else:
            return False

    # Returns true if a given interval (startpos, endpos) matches a GeneItem (exon)
    # The interval start and and can differ from GeneItems start and end by
    # at most allowed_inacc bases
    def equalsItem(self, startpos, endpos, allowed_inacc = DEFAULT_ALLOWED_INACCURACY):
        if startpos < self.start - allowed_inacc:
            return False
        if startpos > self.start + allowed_inacc:
            return False
        if endpos < self.end - allowed_inacc:
            return False
        if endpos > self.end + allowed_inacc:
            return False

        return True

    # Returns true if a given interval (startpos, endpos) and a GeneItem (exon)
    # start at the same position (within allowed_inacc)
    # NOTE: Consider if it might be usefull to also look at the end of the interval
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


    # Returns true if a given interval (startpos, endpos) overlaps a GeneItem (exon)
    # The ovelap size must be at least allowed_inacc bases
    # This implementation of overlap, will include inside
    def overlapsItem(self, startpos, endpos, allowed_inacc = DEFAULT_ALLOWED_INACCURACY):
        if (endpos <= self.start + allowed_inacc) or (startpos >= self.end - allowed_inacc):
            return False
        else:
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


class GeneDescription:
    def __init__(self):
        self.seqname = ''
        self.source = ''
        self.genename = ''
        self.transcriptname = ''
        self.strand = GFF_STRANDFW
        self.start = -1
        self.end = -1
        self.score = 0.0
        self.items = []


    def getLength(self):
        return self.end - self.start

    def insideGene(self, startpos, endpos):
        if startpos >= self.start and endpos <= self.end:
            return True
        else:
            return False

    # this implementation of overlap, will include inside
    def overlapsGene(self, startpos, endpos):
        if endpos <= self.start or startpos >= self.end:
            return False
        else:
            return True

    def basesInsideGene(self, startpos, endpos):
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

    def insideItems(self, startpos, endpos, allowed_inacc = DEFAULT_ALLOWED_INACCURACY):
        for item in self.items:
            if item.insideItem(startpos, endpos, allowed_inacc):
                return True

        return False

    def overlapsItems(self, startpos, endpos, allowed_inacc = DEFAULT_ALLOWED_INACCURACY):
        for item in self.items:
            if item.overlapsItem(startpos, endpos, allowed_inacc):
                return True

    def basesInsideItems(self, startpos, endpos):
        count = 0
        for item in self.items:
            bases += item.basesInside(startpos, endpos)

        return count

    # Recallculate gen start and end position from exons
    def calcBoundsFromItems(self):
        if len(self.items) == 0:
            pass
        else:
            self.start = self.items[0].start
            self.end = self.items[0].end
            for item in self.items[1:]:
                if item.start < self.start:
                    self.start = item.start
                if item.end > self.end:
                    self.end = item.end


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

class BEDLine:
    def __init__(self):
        self.chrom = ''
        self.chromStart = -1
        self.chromEnd = -1
        self.name = ''
        self.score = -1
        self.strand = GFF_STRANDFW
        self.thickStart = -1
        self.thickEnd = -1
        self.itemRGB = ''
        self.blockCount = 0
        self.blockSizes = []
        self.blockStarts = []


# Opposed to BED data, annotation data here will have absolute positions
def Annotation_From_BED(bedline):
    genedscp = GeneDescription()
    genedscp.seqname = bedline.chrom
    genedscp.genename = bedline.name
    genedscp.score = bedline.score
    genedscp.strand = bedline.strand
    genedscp.start = bedline.chromStart
    genedscp.end = bedline.chromEnd
    for i in range(bedline.blockCount):
        geneitem = GeneItem()
        geneitem.start = genedscp.start + bedline.blockStarts[i]
        geneitem.end = geneitem.start + bedline.blockSizes[i]
        genedscp.items.append(geneitem)

    return genedscp


# TODO: the implementation is currently Faulty
# It assumes one exon per gene, which is true e.g. for bacteria but not for eucaryota
def Annotation_From_GFF(gffline):
    genedscp = GeneDescription()
    genedscp.seqname = gffline.seqname
    genedscp.source = gffline.source
    genedscp.start = gffline.start
    genedscp.end = gffline.end + 1
    genedscp.score = gffline.score
    genedscp.strand = gffline.strand

    # Extracting from GFF attributes
    # Removing double quotes!

    genedscp.transcriptname = gffline.attribute['transcript_id'][1:-1]
    genedscp.genename = gffline.attribute['gene_id'][1:-1]

    # constructing a single gene item (exon)
    geneitem = GeneItem()
    geneitem.frame = gffline.frame
    geneitem.start = gffline.start
    geneitem.end = gffline.end + 1
    genedscp.items.append(geneitem)

    return genedscp


def Load_Annotation_From_File(filename, check_duplicates = False):

    fname, fext = os.path.splitext(filename)
    if fext == '.gff':
        istype = 'GFF'
    elif fext == '.gtf':
        istype = 'GTF'
    elif fext == '.bed':
        istype = 'BED'
    else:
        raise Exception('Invalid annotation file istype: %s' % fext)

    annotation_dict = {}

    annotations = []

    # Process GFF lines, several lines represent the same transcript
    # Collect the lines with the same transcript name as the same annotation
    # Since there can be more that one annotation with the same name, have to watch out
    # Colect only consequtive enteries.
    if istype == 'GFF' or istype == 'GTF':
        gff_lines = Load_GFF_From_File(filename)
        old_annt_name = ''
        curr_annt = None        # Current collected annotation
        for gffline in gff_lines:
            new_annt = Annotation_From_GFF(gffline)
            new_annt_name = new_annt.transcriptname
            if old_annt_name != new_annt_name:
                if old_annt_name != '':
                    # Store the last collected annotation (calculate start and end position from items first)
                    curr_annt.calcBoundsFromItems()
                    annotations.append(curr_annt)
                # Start a new collected annotation
                curr_annt = new_annt
            else:
                if new_annt.seqname != curr_annt.seqname or \
                   new_annt.source != curr_annt.source or \
                   new_annt.strand != curr_annt.strand or \
                   new_annt.genename != curr_annt.genename or \
                   new_annt.transcriptname != curr_annt.transcriptname:
                    raise Exception('Invalid GFF/GTF line for transcript %s' % new_annt_name)
                # Assuming that new_annt has only one item
                assert len(new_annt.items) == 1
                curr_annt.items.append(new_annt.items[0])

            old_annt_name = new_annt_name

        # Add the last collected annotation
        if old_annt_name != '':
            annotations.append(curr_annt)

    elif istype == 'BED':
        bed_lines = Load_BED_From_File(filename)
        for bedline in bed_lines:
            annt = Annotation_From_BED(bedline)
            annotations.append(annt)

    # Checking annotations for dupicate genenames
    # Raising exception if finding any
    if check_duplicates:
        num_duplicates = 0
        # duplicates = []
        for i in xrange(len(annotations)):
            genename1 = annotations[i].genename
            for j in xrange(i+1, len(annotations)):
                genename2 = annotations[j].genename
                if genename1 == genename2:
                    import pdb
                    pdb.set_trace()
                    num_duplicates += 1
                    # duplicates.append(genename1)

        if num_duplicates > 0:
            raise Exception('Duplicate annotations found (%d)' % num_duplicates)

    return annotations


def Load_GFF_From_File(filename):
    gff_lines = []
    if not (filename.endswith('.gff') or filename.endswith('.gtf')):
        sys.stderr.write('\nWARNING: file %s does not have GFF/GTF extension!\n' % filename)

    fname, fext = os.path.splitext(filename)
    if fext == '.gff':
        type = 'GFF'
    elif fext == '.gtf':
        type = 'GTF'
    else:
        raise Exception('Invalid annotation file type: %s' % fext)

    file = open(filename, 'rU')
    for line in file:
        elements = line.split('\t')
        gffline = GFFLine()

        if elements[0] == '.':
            gffline.seqname = ''
        else:
            gffline.seqname = elements[0]
        if elements[1] == '.':
            gffline.source = ''
        else:
            gffline.source = elements[1]
        if elements[2] == '.':
            gffline.feature = ''
        else:
            gffline.feature = elements[2]
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
            att_line = elements[8]
            att_list = att_line.split(';')          # Separating attribute definitions
            for i in xrange(len(att_list)):
                elements = att_list[i].split()      # Separating key and value for each attribute
                if len(elements) == 2:
                    gffline.attribute[elements[0]] = elements[1]

        # TODO: GFF and GTF contain start and stop codons, CDSs and exons
        # currently using only exons (maybe CDS would be a better choice)
        if gffline.feature == 'exon':
            gff_lines.append(gffline)

    return gff_lines


def Load_BED_From_File(filename):
    bed_lines = []
    if not (filename.endswith('.bed')):
        sys.stderr.write('\nWARNING: file %s does not have BED extension!\n' % filename)

    # Copied from GFF, might be useful in the future
    fname, fext = os.path.splitext(filename)
    if fext == '.bed':
        type = 'BED'
    else:
        raise Exception('Invalid annotation file type: %s' % fext)

    file = open(filename, 'rU')
    for line in file:
        # Ignoring header lines
        if line.startswith('#') or line.startswith('track') or line.startswith('browser'):
            pass
        else:
            elements = line.split()    # splitting with default delimitters
            attcount = len(elements)
            bedline = BEDLine()
            bedline.chrom = elements[0]
            bedline.chromStart = int(elements[1])
            bedline.chromEnd = int(elements[2])
            if attcount >= 4:
                bedline.name = elements[3]
            if attcount >= 5:
                bedline.score = int(elements[4])
            if attcount >= 6:
                bedline.strand = elements[5]
            if attcount >= 7:
                bedline.thickStart = int(elements[6])
            if attcount >= 8:
                bedline.thickEnd = int(elements[7])
            if attcount >= 9:
                bedline.itemRGB = elements[8]
            if attcount >= 10:
                bedline.blockCount = int(elements[9])
            if attcount >= 11:
                if elements[10].endswith(','):
                    elements[10] = elements[10][:-1]
                bedline.blockSizes = [int(el) for el in elements[10].split(',')]
            if attcount >= 12:
                if elements[11].endswith(','):
                    elements[11] = elements[11][:-1]
                bedline.blockStarts = [int(el) for el in elements[11].split(',')]

            bed_lines.append(bedline)

    return bed_lines



if __name__ == "__main__":
	#pass;
        filename = sys.argv[1]
        annotations = Load_Annotation_From_File(filename)

        for tran in annotations:
            print tran.transcriptname, tran.genename, len(tran.items)
