# Debating whether or not to attempt Bowtie2 mapping for speed.
# Local would suck for install.
# Remote would suck for maintenance.

import subprocess

indexLookup = {
    'mouse': '../indices/mm10/mm10'
}

def genomemask(fasta_string,species="mouse"):

    #res = subprocess.call(["bowtie2", "-x", indexLookup[species], "-f", fastqFile1, "S", "Outpufile.sam"])
    pass
