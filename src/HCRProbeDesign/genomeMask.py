# Debating whether or not to attempt Bowtie2 mapping for speed.
# Local would suck for install.
# Remote would suck for maintenance.

import subprocess
import tempfile
import pysam
from collections import defaultdict
import os
import urllib.request
import shutil
import pkg_resources
from zipfile import ZipFile
import argparse

package_directory = os.path.dirname(os.path.abspath(__file__))
indices_directory = f'{package_directory}/indices/'
indexLookup = {
    'mouse': os.path.join(indices_directory,'mm10/mm10')
}

#TODO: make genomemask() take transient index argment if not default in species
def genomemask(fasta_string,handleName="tmp",species="mouse",nAlignments = 3, index=None):
    fasta_file = f'{handleName}_reads.fa'
    tmpFasta = open(fasta_file,mode="w")
    #print(tmpFasta.name)
    tmpFasta.write(fasta_string)
    tmpFasta.close()
    sam_file = f'{handleName}.sam'
    print(indexLookup[species])
    res = subprocess.call(["bowtie2", f"-k{nAlignments}", "-x", indexLookup[species], "-f", fasta_file, "-S", sam_file])
    return res

def countHitsFromSam(samFile):
    hitCounts = defaultdict(int)
    sam = pysam.AlignmentFile(samFile,"r")
    for read in sam.fetch():
        #print(read.get_tag(tag="TM:i")) # This doesn't work because pysam doesn't read the XS:i: tag
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

def install_index(url='https://genome-idx.s3.amazonaws.com/bt/mm10.zip',genome="mm10"):
    # Parse arguments if any
    #Argument handling
	parser = argparse.ArgumentParser(description="Bowtie2 index retrieval and installation",formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    args = parser.parse_args()

    print(f'Downloading Bowtie2 index from {url} ...')
    #folder = pkg_resources.resource_filename('HCRProbeDesign','indices/')
    index_folder = indices_directory
    fname = f'{index_folder}{genome}.zip'
    print(fname)
    if not os.path.exists(fname):
        urllib.request.urlretrieve(url, fname)
    with ZipFile(fname, 'r') as zip:
        # printing all the contents of the zip file
        #zip.printdir()
        # extracting all the files
        print(f'Extracting index from {fname} to {index_folder}{genome}...')
        zip.extractall(f'{index_folder}{genome}')
    print('Done!')
