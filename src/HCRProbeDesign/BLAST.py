"""NCBI BLAST utilities for probe sequence validation."""

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from . import utils
import time

def blastProbes(fasta_string,species="mouse",verbose=True):
    """
    Submit a BLASTN job for the given FASTA string.

    :param fasta_string: FASTA-formatted string containing probe sequences.
    :param species: Species key for Entrez restriction (mouse or human).
    :param verbose: Emit progress messages to stderr.
    :return: NCBIWWW result handle for subsequent parsing.
    """
    speciesLookup = {"mouse": "Mus musculus",
                    "human": "Homo sapiens"
                    }
    entrezString = f'{speciesLookup[species]} [ORGN]'
    start = time.time()
    result_handle = NCBIWWW.qblast("blastn", "nt", fasta_string, entrez_query = entrezString)
    end = time.time()
    utils.eprint(f'BLAST took {end - start} seconds...')

    if verbose:
        utils.eprint("BLASTN Success!")
    return result_handle

def getNHits(blast_handle, verbose=True):
    """
    Report the number of hits for each record in a BLAST response.

    :param blast_handle: Handle returned by NCBIWWW.qblast.
    :param verbose: Reserved for future verbosity controls.
    :return: None.
    """
    blast_records = NCBIXML.parse(blast_handle)
    for blast_record in blast_records:
        utils.eprint(blast_record.num_hits)
