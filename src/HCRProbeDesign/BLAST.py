import time

from Bio.Blast import NCBIWWW, NCBIXML

from . import utils


def blastProbes(fasta_string, species="mouse", verbose=True):
    speciesLookup = {"mouse": "Mus musculus", "human": "Homo sapiens"}
    entrezString = f"{speciesLookup[species]} [ORGN]"
    start = time.time()
    result_handle = NCBIWWW.qblast(
        "blastn", "nt", fasta_string, entrez_query=entrezString
    )
    end = time.time()
    utils.eprint(f"BLAST took {end - start} seconds...")

    if verbose:
        utils.eprint("BLASTN Success!")
    return result_handle


def getNHits(blast_handle, verbose=True):
    blast_records = NCBIXML.parse(blast_handle)
    for blast_record in blast_records:
        utils.eprint(blast_record.num_hits)
    for blast_record in blast_records:
        utils.eprint(blast_record.num_hits)
