# Debating whether or not to attempt Bowtie2 mapping for speed.
# Local would suck for install.
# Remote would suck for maintenance.

import subprocess
import tempfile
import pysam
from collections import defaultdict
import os

package_directory = os.path.dirname(os.path.abspath(__file__))

indexLookup = {
    'mouse': os.path.join(package_directory,'../indices/mm10/mm10')
}

def genomemask(fasta_string,handleName="tmp",species="mouse"):
    fasta_file = f'{handleName}_reads.fa'
    tmpFasta = open(fasta_file,mode="w")
    #print(tmpFasta.name)
    tmpFasta.write(fasta_string)
    tmpFasta.close()
    sam_file = f'{handleName}.sam'
    print(indexLookup[species])
    res = subprocess.call(["bowtie2", "-x", indexLookup[species], "-f", fasta_file, "-S", sam_file])
    return res

def countHitsFromSam(samFile):
    hitCounts = defaultdict(int)
    sam = pysam.AlignmentFile(samFile,"r")
    for read in sam.fetch():
        #print(read)
        if read.is_unmapped:
            hitCounts[read.query_name] += 0
        else:
            hitCounts[read.query_name] += 1
    return hitCounts

def test():
    fasta_string=">read1\ntacgagcttactggacgagcgtgactctgac\n>read2\nctgagctgatgcgacgtatctgatgctgtacgtgacg\n>read3\ncgagctatcgtactgagcggagcgcgggcgatat"
    handleName = 'test'
    proc_info = genomemask(fasta_string,handleName=handleName,species="mouse")
    res = countHitsFromSam(f'{handleName}.sam')
    print(res)
