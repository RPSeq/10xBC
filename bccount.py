import sys
import argparse
import pysam
from collections import namedtuple, defaultdict

__author__ = "Ryan Smith (ryanpsmith@wustl.edu)"
__version__ = "$Revision: 0.1.0 $"
__date__ = "$Date: 2017-05-16 $"

def get_args():
    """define command line arguments"""

    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawTextHelpFormatter,
        description = "author: " + __author__ + 
        "\nversion: " + __version__ + 
        "\ndescription: DESCRIPTION"
        )

    parser.add_argument(
        '-i',
        '--input',
        type = str,
        required = True,
        help = "Input .bam (required index)"
        )

    parser.add_argument(
        '-o',
        '--output',
        type = argparse.FileType('w'),
        required = False,
        default = sys.stdout,
        help = "Output file [stdout]"
        )

    parser.add_argument(
        '-s', 
        '--sites', 
        type = argparse.FileType('r'), 
        required = False, 
        default = None, 
        help = "Sites from which to collect BX metrics [stdin]"
        )

    parser.add_argument(
        '-d', 
        '--depth', 
        type = int, 
        required = False, 
        default = 500, 
        help = "Max read depth at site [500]"
        )

    # parse the arguments
    args = parser.parse_args()

    # if no input,
    if args.sites == None:

        # check if part of pipe,
        if sys.stdin.isatty():

            # if not, help
            parser.print_help()
            exit(1)
        
        # if so, read stdin
        else:
            args.sites = sys.stdin

    return args

class BXVar:
    """class representing a collection of BX reads at a given site"""

    def __init__(self, site, al = None):

        # BED record
        self.site = site

        # map of {basecall: {bx: count}}
        self.counts = defaultdict(lambda: defaultdict(int))

        # map of bx: count
        self.barcodes = defaultdict(int)

        # num unique barcodes
        self.n = 0

        # add al if provided
        if al:
            self.add_al(al)

    def add_al(self, pread):
        """add al to collection"""

        # get al from pileupread
        al = pread.alignment

        # get barcode
        bx = get_bx(al)

        # increment unique barcode count
        if bx != "NA" and bx not in self.barcodes:
            self.n += 1

        #get base call from read
        pbase = al.query_sequence[pread.query_position]

        # increment basecall bx count
        self.counts[pbase][bx] += 1

        # increment unique bx count (total)
        self.barcodes[bx] += 1

    def bx_counts(self):
        return ";".join([x+","+str(self.barcodes[x]) for x in self.barcodes])




def get_bx(al):
    """get 10X BX barcode tag from alignment"""

    # return NA if no tag
    try:
        return al.get_tag("BX")
    except:
        return "NA"


def bed(infile):
    """generator for parsing bed files"""

    # define namedtuple for BED records
    BED = namedtuple("BED", "chrom, start, end, ref, alt")

    # yield BED records for each entry in file
    for line in infile:
        line = line.strip().split()

        # len check
        assert len(line) == 5, \
            "Input bed must have chrom, start, end, ref, alt" 

        # set BED fields
        yield BED(
            chrom = line[0],
            start = int(line[1]),
            end = int(line[2]),
            ref = line[3],
            alt = line[4]
            )

def main():
    args = get_args()

    in_bam = pysam.AlignmentFile(args.input, "rb")

    in_bed = bed(args.sites)

    #for each site in coordinates file
    for site in in_bed:

        #generate a truncated pileup
        pile = in_bam.pileup(
            reference = site.chrom,
            start = site.start,
            end = site.end,
            stepper = "all",    #default filter for alignments
            truncate = True,    #only include exact columns of site
            max_depth = args.depth
            )

        #for each column (only 1 for single-base truncated pileups)
        for column in pile:

            # generate BXVar
            var = BXVar(site)

            #iterate over the pileup reads
            for pread in column.pileups:

                # del or ins alignments will not have a called base
                if not pread.is_del and not pread.is_refskip:

                    var.add_al(pread)

            print site.chrom, site.start, var.n, var.bx_counts()



    
if __name__ == '__main__':
    main()