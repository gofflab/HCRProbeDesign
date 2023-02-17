#!/usr/bin/env python3
import argparse
import re

# from probeDesign import BLAST
import sys
from itertools import product
from string import ascii_uppercase

# from Bio.Seq import Seq
import primer3

from . import HCR, genomeMask, repeatMask, sequencelib, utils
from .tiles import Tile, TileError

# import copy
# import string
from .utils import pp

#######################
# Scan input sequence #
#######################
# TODO: Modify this so that it only gets the window of appropriate size. We will add prefix and suffix afterwards.


def scanSequence(sequence, seqName, tileStep=1, tileSize=52):
    tiles = []
    # Pre-compute number of chunks to emit
    numOfChunks = int(((len(sequence) - tileSize) / tileStep) + 1)

    # Tile across reverse complement of sequence
    for i in range(0, numOfChunks * tileStep, tileStep):
        tile = Tile(
            sequence=sequencelib.reverse_complement(sequence[i : i + tileSize]),
            seqName=seqName,
            startPos=i + 1,
        )
        if not tile.isMasked():
            tiles.append(tile)
    return tiles


###################
# Reporting
###################


def outputTable(tiles, outHandle=sys.stdout):
    """
    Formats tile output and writes to outHandle
    """
    outputKeys = [
        "name",
        "probe",
        "start",
        "length",
        "P1",
        "P2",
        "channel",
        "GC",
        "Tm",
        "dTm",
        "GibbsFE",
    ]
    outHandle.write("\t".join(outputKeys) + "\n")
    for tile in tiles:
        outHandle.write(
            f"{tile.name}\t{tile.sequence}\t{tile.start}\t{len(tile)}\t{tile.P1}\t{tile.P2}\t{tile.channel}\t{tile.GC():.2f}\t{primer3.calcTm(tile.sequence):.2f}\t{tile.dTm:.2f}\t{tile.Gibbs:.2f}\n"
        )


def outputIDT(tiles, outHandle=sys.stdout):
    """
    Formats tile output for direct ordering using IDT template
    """
    # 96-well addressing
    rows = list(ascii_uppercase[0:8])
    columns = [x + 1 for x in range(12)]
    outputKeys = ["name", "start", "length", "P1", "P2", "channel"]
    # Header for IDT plate template
    outHandle.write("\t".join(["Name", "Sequence"]) + "\n")
    # One row per oligo (P1='odd', p2='even')
    odd = []
    even = []
    for tile in tiles:
        odd.append((f"{tile.name}:{tile.channel}:odd", tile.P1))
        even.append((f"{tile.name}:{tile.channel}:even", tile.P2))

    for oligo in odd:
        outHandle.write(f"{oligo[0]}\t{oligo[1]}\n")

    for oligo in even:
        outHandle.write(f"{oligo[0]}\t{oligo[1]}\n")


#
# def alignOutput(inseq,tiles):
#     """Uses tile information to make a nice output w/ probes aligned to inseq
#
#     Returns 3 elements in a list:
#         [0] - the inseq
#         [1] - the oligos
#         [2] - probe # information
#     """
#     nTiles = len(tiles)
#     compoligos  = [seq.complement(inseq[j:(j+len(j))]) for j in tiles]
#     spaceoligos = [' '*(tiles[i]-tiles[i-1]-len(tiles[i])) for i in range(1,nTiles)]
#     spaceoligos = [' '*tiles[0]] + spaceoligos
#     probenum    = [ ('Probe # '+str(i+1)).ljust(len(tiles[i])) for i in range(nTiles)]
#
#     compseq  = ''.join([spaceoligos[i] + compoligos[i] for i in range(nTiles)])
#     probeseq = ''.join([spaceoligos[i] + probenum[i] for i in range(nTiles)])
#
#     compseq = compseq.ljust(len(inseq))
#     probeseq = probeseq.ljust(len(inseq))
#
#     return [inseq, compseq, probeseq]


def calcOligoCost(tiles, pricePerBase=0.19):
    total = 0.0
    for tile in tiles:
        probeSize = len(tile.P1) + len(tile.P2)
        total += probeSize * pricePerBase
    return total


def outputRunParams(args):
    utils.eprint(f"\nParameters:")
    utils.eprint(print(args))


###############
# Test function with hard-coded variables
###############


def test():
    # Set default args
    minGC = 45.0
    maxGC = 55.0
    targetGC = 50.0
    dTmMax = 5.0
    minGibbs = -70.0
    maxGibbs = -50.0
    targetGibbs = -60.0
    tileSize = 52
    verbose = True
    species = "mouse"
    genomemask = True
    num_hits_allowed = 1
    channel = "B1"
    dTmFilter = False

    # Read custom args

    # Parse fasta files and loop over records
    utils.eprint("Reading in Fasta file")
    # fname = "test/eGFP.fa"
    fname = "test/tdTomato.fa"
    targetName = "tdTomato"
    handle = open(fname, "r")
    fastaIter = sequencelib.FastaIterator(handle)

    mySeq = next(fastaIter)  # TODO: make loopable when migrating to main()

    # RepeatMasking
    # TODO: Specify DNA source in input params
    utils.eprint("\nRepeat Masking...")
    mySeq["sequence"] = repeatMask.repeatmask(mySeq["sequence"], dnasource=species)

    # Convert to lowercase
    # mySeq['sequence'] = mySeq['sequence'].lower()

    # Tile over masked sequence record to generate all possible probes of appropriate length that are not already masked
    tiles = scanSequence(
        mySeq["sequence"], mySeq["name"], tileStep=1, tileSize=tileSize
    )
    utils.eprint(f"{len(tiles)} tiles available of length {tileSize}...")

    # Check for invalid characters

    # Crunmask
    utils.eprint("\nChecking for runs of C's")
    tiles = [
        tile
        for tile in tiles
        if not tile.hasRuns(runChar="c", runLength=7, mismatches=2)
    ]
    utils.eprint(f"{len(tiles)} tiles remain")

    # Grunmask
    utils.eprint("\nChecking for runs of G's")
    tiles = [
        tile
        for tile in tiles
        if not tile.hasRuns(runChar="g", runLength=7, mismatches=2)
    ]
    utils.eprint(f"{len(tiles)} tiles remain")

    # Calculate Hairpins
    utils.eprint("\nChecking for hairpins")
    for tile in tiles:
        thermRes = primer3.calcHairpin(tile.sequence)
    tiles = [
        tile for tile in tiles if primer3.calcHairpin(tile.sequence).structure_found
    ]
    utils.eprint(f"{len(tiles)} tiles remain")

    print(len(tiles))

    # PseudogeneMasking.  BLAST?

    # GenomeMasking?  Using bowtie because BLAST over WWW is unpredictable
    if genomemask:
        utils.eprint(
            f"\nChecking unique mapping of remaining tiles against {species} reference genome"
        )
        blast_string = "\n".join([tile.toFasta() for tile in tiles])
        # if verbose:
        # 	utils.eprint("Submitting BLAST request with the following string")
        # 	utils.eprint(blast_string)
        blast_res = genomeMask.genomemask(
            blast_string, handleName=targetName, species=species
        )
        utils.eprint(f"Parsing BLAST output now")
        hitCounts = genomeMask.countHitsFromSam(f"{targetName}.sam")
        # print(hitCounts)
        # Check that keys returned from hitCounts match order of tiles in tiles
        assert all(
            map(
                lambda x, y: x == y,
                [k for k in hitCounts.keys()],
                [tile.name for tile in tiles],
            )
        )
        utils.eprint(
            f"Filtering for <= {num_hits_allowed} alignments to {species} genome..."
        )
        for i in range(len(tiles)):
            k = list(hitCounts.keys())[i]
            tile.hitCount = hitCounts[k]
        tiles = [tile for tile in tiles if tile.hitCount <= num_hits_allowed]
        utils.eprint(f"{len(tiles)} tiles remain")

    # TM filtering

    # GC filtering
    utils.eprint(f"\nChecking for {minGC} < GC < {maxGC}")
    tiles = [tile for tile in tiles if tile.GC() >= minGC]
    tiles = [tile for tile in tiles if tile.GC() <= maxGC]
    utils.eprint(f"{len(tiles)} tiles remain")

    # Gibbs filtering
    utils.eprint(f"\nChecking for {minGibbs} < Gibbs FE < {maxGibbs}")
    [tile.calcGibbs() for tile in tiles]
    tiles = [tile for tile in tiles if tile.Gibbs >= minGibbs]
    tiles = [tile for tile in tiles if tile.Gibbs <= maxGibbs]
    utils.eprint(f"{len(tiles)} tiles remain")

    # Split tile into probeset
    utils.eprint(f"\nSplitting tiles into probesets")
    [tile.splitProbe() for tile in tiles]
    [tile.calcdTm() for tile in tiles]

    # dTm between halves
    if dTmFilter:
        utils.eprint(f"\nChecking for dTm <= {dTmMax} between probes for each tile")
        tiles = [tile for tile in tiles if tile.dTm <= dTmMax]
        utils.eprint(f"{len(tiles)} tiles remain")

    # Break remaining probes into non-overlapping regions
    # TODO: there must be a better way to do this to minimize overlaps.  Perhaps testing overlaps from bestTiles later on?  Would ensure better quality picks make it to the end.
    regions = {}
    regionCount = 0

    regionList = [tiles[0]]
    for i in range(1, len(tiles)):
        if i == len(tiles) - 1:
            regions[ascii_uppercase[regionCount]] = regionList
            break
        if tiles[i].overlaps(tiles[i - 1]):
            regionList.append(tiles[i])
        else:
            regions[ascii_uppercase[regionCount]] = regionList
            regionList = [tiles[i]]
            regionCount = regionCount + 1

    # Select best from each region
    bestTiles = []
    for k, v in regions.items():
        bestTiles.append(
            v[
                min(
                    range(len(v)),
                    key=lambda i: abs([x.Gibbs for x in v][i] - targetGibbs),
                )
            ]
        )

    # Add initator and spacers to split probes
    utils.eprint(
        f"\nAdding spacers and initiator sequences to split probes for channel {channel}"
    )
    [tile.makeProbes(channel) for tile in bestTiles]

    # Print out results
    for tile in bestTiles:
        print(
            f"{tile}\tP1_sequence:{tile.P1}\tP2_sequence:{tile.P2}\tmyTm:{tile.Tm():.2f}\tprimer3-Tm:{primer3.calcTm(tile.sequence):.2f}\tdTm:{tile.dTm:.2f}\tGC%:{tile.GC():.2f}\tGibbs:{tile.Gibbs:.2f}"
        )

    utils.eprint(f"\nThere are a total of {len(bestTiles)} best probes")

    # TODO: dump output to designated file handles


def main():
    """
    Main function for HCR Probe design.  Called when used directly from cmdline
    """
    #######################
    # Variables
    #######################
    # # Set default args
    # tileSize = 52
    # tileStep = 1
    # minGC = 45.0
    # maxGC = 55.0
    # targetGC = 50.0
    # dTmMax = 5.0
    # minGibbs = -70.0
    # maxGibbs = -50.0
    # targetGibbs = -60.0
    # verbose = True
    # species = 'mouse'
    # genomemask = False
    # repeatmask = True
    # blastmask = False
    # pseudogenemask = False
    # channel = "B1"
    # num_hits_allowed = 1
    # dTmFilter = False

    # Argument handling
    parser = argparse.ArgumentParser(
        description="Probe Design Utility for HCR v3.0",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "infile",
        help="Properly formatted fasta file against which to design probes",
        type=argparse.FileType("r"),
    )
    parser.add_argument("-v", "--verbose", help="Verbose output", action="store_true")
    parser.add_argument(
        "-c",
        "--channel",
        help="HCR Channel initiator sequences",
        default="B1",
        choices=HCR.initiators.keys(),
    )
    parser.add_argument(
        "-o",
        "--output",
        help="Output file name",
        nargs="?",
        type=argparse.FileType("w"),
        default=sys.stdout,
    )
    parser.add_argument(
        "--tileSize",
        help="Size of the tiles along the target sequence",
        type=int,
        default=52,
    )
    parser.add_argument(
        "--targetName",
        help="User-friendly name for target sequence (e.g. Gene Name)",
        default="target",
    )
    parser.add_argument(
        "-s", "--species", help="Species for repeatmask and genomemask", default="mouse"
    )
    parser.add_argument("--minGC", help="Min allowable GC", default=45.0, type=float)
    parser.add_argument("--maxGC", help="Max allowable GC", default=55.0, type=float)
    parser.add_argument("--targetGC", help="Target GC", default=50.0, type=float)
    parser.add_argument(
        "--dTmMax",
        help="Max allowable difference in Tm between probes in set",
        default=5.0,
        type=float,
    )
    parser.add_argument(
        "--dTmFilter",
        help="Enable filtering based on dTm between probeset halves.",
        default=False,
        action="store_true",
    )
    parser.add_argument(
        "-g",
        "--no-genomemask",
        help="Disables bowtie2 checking for multiple hits to genome",
        default=True,
        action="store_false",
    )
    parser.add_argument(
        "-i", "--index", help="Location of bowtie2 index file for genomemask analysis"
    )
    parser.add_argument(
        "-r",
        "--no-repeatmask",
        help="Disables repeatmasker masking of target sequence",
        default=True,
        action="store_false",
    )
    parser.add_argument(
        "--minGibbs", help="Min allowable GibbsFE", default=-70.0, type=float
    )
    parser.add_argument(
        "--maxGibbs", help="Max allowable GibbsFE", default=-50.0, type=float
    )
    parser.add_argument(
        "--targetGibbs", help="Target GibbsFE", default=-60.0, type=float
    )
    parser.add_argument(
        "--maxRunLength", help="Max allowable homopolymer run size", default=7, type=int
    )
    parser.add_argument(
        "-n", "--maxProbes", help="Max number of probes to return", default=20, type=int
    )
    parser.add_argument(
        "--maxRunMismatches",
        help="Max allowable homopolymer run mismatches",
        default=2,
        type=int,
    )
    parser.add_argument(
        "--num-hits-allowed",
        help="Number of allowable hits to genome",
        default=1,
        type=int,
    )
    parser.add_argument(
        "--idt",
        help="File name to output tsv format optimized for IDT ordering",
        type=argparse.FileType("w"),
        default=None,
    )
    parser.add_argument(
        "--calcPrice",
        help="Calculate total cost of probe synthesis assuming $0.19 per base",
        default=False,
        action="store_true",
    )
    args = parser.parse_args()

    #########
    # Parse fasta file. Currently not looping over records, only uses first fasta record
    #########
    utils.eprint("Reading in Fasta file")
    handle = args.infile
    fastaIter = sequencelib.FastaIterator(handle)

    mySeq = next(fastaIter)  # TODO: make loopable when migrating to main()

    #############
    # Repeatmask target sequence
    #############
    if args.no_repeatmask:
        # RepeatMasking
        utils.eprint(f"\nRepeat Masking using {args.species} reference...")
        mySeq["sequence"] = repeatMask.repeatmask(
            mySeq["sequence"], dnasource=args.species
        )

    # Check for invalid characters ?

    ###############
    # Tile over masked sequence record to generate all possible probes of appropriate length that are not already masked
    ###############
    utils.eprint(
        f"\nBreaking target sequence into revcomp tiles of size {args.tileSize}..."
    )
    tiles = scanSequence(
        mySeq["sequence"], mySeq["name"], tileStep=1, tileSize=args.tileSize
    )  # Here we remove masked sequences and rev comp for tiles.
    utils.eprint(f"{len(tiles)} tiles available of length {args.tileSize}...")

    ##############
    # Crunmask
    ##############
    utils.eprint("\nChecking for runs of C's")
    tiles = [
        tile
        for tile in tiles
        if not tile.hasRuns(
            runChar="c", runLength=args.maxRunLength, mismatches=args.maxRunMismatches
        )
    ]
    utils.eprint(f"{len(tiles)} tiles remain")

    ##############
    # Grunmask
    ##############
    utils.eprint("\nChecking for runs of G's")
    tiles = [
        tile
        for tile in tiles
        if not tile.hasRuns(
            runChar="g", runLength=args.maxRunLength, mismatches=args.maxRunMismatches
        )
    ]
    utils.eprint(f"{len(tiles)} tiles remain")

    ##############
    # Calculate Hairpins
    ##############
    utils.eprint("\nChecking for hairpins")
    # TODO: add this as a user-selectable parameter. Currently awkward as we don't have a min Tm filter.
    max_Th = 45.0  # This is the maximum tolerated calculated melting temperature of any predicted hairpins
    tiles = [
        tile
        for tile in tiles
        if primer3.calcHairpin(tile.sequence).tm < max_Th
        or not primer3.calcHairpin(tile.sequence).structure_found
    ]
    utils.eprint(f"{len(tiles)} tiles remain")

    ##############
    # GenomeMasking?  Using bowtie because BLAST over WWW is unpredictable
    ##############
    if args.no_genomemask:
        utils.eprint(
            f"\nChecking unique mapping of remaining tiles against {args.species} reference genome"
        )
        blast_string = "\n".join([tile.toFasta() for tile in tiles])
        blast_res = genomeMask.genomemask(
            blast_string, handleName=args.targetName, species=args.species, index=None
        )
        utils.eprint(f"Parsing bowtie2 output now")
        hitCounts = genomeMask.countHitsFromSam(f"{args.targetName}.sam")
        # print(hitCounts)
        # Check that keys returned from hitCounts match order of tiles in tiles
        assert all(
            map(
                lambda x, y: x == y,
                [k for k in hitCounts.keys()],
                [tile.name for tile in tiles],
            )
        )
        utils.eprint(
            f"Filtering for <= {args.num_hits_allowed} alignments to {args.species} genome..."
        )
        for i in range(len(tiles)):
            k = list(hitCounts.keys())[i]
            tiles[i].hitCount = hitCounts[k]
            # utils.eprint(f'{tiles[i].hitCount}')
        tiles = [tile for tile in tiles if tile.hitCount <= args.num_hits_allowed]
        utils.eprint(f"{len(tiles)} tiles remain")

    ###############
    # TM filtering
    ###############

    ###############
    # GC filtering
    ###############
    utils.eprint(f"\nChecking for {args.minGC} < GC < {args.maxGC}")
    tiles = [tile for tile in tiles if tile.GC() >= args.minGC]
    tiles = [tile for tile in tiles if tile.GC() <= args.maxGC]
    utils.eprint(f"{len(tiles)} tiles remain")

    ###############
    # Gibbs filtering
    ###############
    utils.eprint(f"\nChecking for {args.minGibbs} < Gibbs FE < {args.maxGibbs}")
    [tile.calcGibbs() for tile in tiles]
    tiles = [tile for tile in tiles if tile.Gibbs >= args.minGibbs]
    tiles = [tile for tile in tiles if tile.Gibbs <= args.maxGibbs]
    utils.eprint(f"{len(tiles)} tiles remain")

    ###############
    # Split tile into probeset
    ###############
    utils.eprint(f"\nSplitting tiles into probesets")
    [tile.splitProbe() for tile in tiles]
    [tile.calcdTm() for tile in tiles]

    ###############
    # dTm between halves
    ###############
    if args.dTmFilter:
        utils.eprint(
            f"\nChecking for dTm <= {args.dTmMax} between probes for each tile"
        )
        tiles = [tile for tile in tiles if tile.dTm <= args.dTmMax]
        utils.eprint(f"{len(tiles)} tiles remain")

    ###############
    # Break remaining probes into non-overlapping regions
    ###############
    # #TODO: there must be a better way to do this to minimize overlaps.  Perhaps testing overlaps from bestTiles later on?  Would ensure better quality picks make it to the end.
    # regions = {}
    # regionCount = 0
    #
    # regionList = [tiles[0]]
    # for i in range(1,len(tiles)):
    # 	if i == len(tiles)-1:
    # 		regions[ascii_uppercase[regionCount]] = regionList
    # 		break
    # 	if tiles[i].overlaps(tiles[i-1]):
    # 		regionList.append(tiles[i])
    # 	else:
    # 		regions[ascii_uppercase[regionCount]] = regionList
    # 		regionList = [tiles[i]]
    # 		regionCount = regionCount + 1
    #
    # ################
    # # Select best from each region
    # ################
    # bestTiles = []
    # for k,v in regions.items():
    # 	bestTiles.append(v[min(range(len(v)), key=lambda i: abs([x.Gibbs for x in v][i]-args.targetGibbs))])

    ################
    # Select overall best n tiles (regardless of region)
    ################
    # Instead of above 'region-based' approach, Let's start by choosing the best probes (by min distance to targetGibbs and/or min distance to targetGC).
    # As we choose subsequent best probes, test for overlap with any existing tiles (tile.overlaps(tile2)).  If none, then append next best to bestTiles
    # Do this until you are out of tiles or bestTiles reaches a certain number of tiles.

    # TODO: Currently ranking tiles based on min distance to targetGibbs.  Need to make an argument to select targetGC as goal instead.
    utils.eprint(
        f"\nSelecting top {args.maxProbes} tiles based on distance to targetGibbs = {args.targetGibbs}"
    )
    bestTiles = []
    while len(bestTiles) < args.maxProbes and len(tiles) > 0:
        nextBestIdx = min(
            range(len(tiles)),
            key=lambda i: abs([x.Gibbs for x in tiles][i] - args.targetGibbs),
        )
        # print(f'{nextBestIdx}')
        if len(bestTiles) == 0:
            bestTiles.append(tiles.pop(nextBestIdx))
            continue
        if any([tiles[nextBestIdx].overlaps(x) for x in bestTiles]):
            tiles.pop(nextBestIdx)
        else:
            bestTiles.append(tiles.pop(nextBestIdx))

    utils.eprint(f"Selected {len(bestTiles)} tiles for probe design")

    ################
    # Add initator and spacers to split probes
    ################
    utils.eprint(
        f"\nAdding spacers and initiator sequences to split probes for channel {args.channel}"
    )
    [tile.makeProbes(args.channel) for tile in bestTiles]

    ################
    # Print out results
    ################
    # for tile in bestTiles:
    # 	print(f"{tile}\tP1_sequence:{tile.P1}\tP2_sequence:{tile.P2}\tmyTm:{tile.Tm():.2f}\tprimer3-Tm:{primer3.calcTm(tile.sequence):.2f}\tdTm:{tile.dTm:.2f}\tGC%:{tile.GC():.2f}\tGibbs:{tile.Gibbs:.2f}")
    if args.idt is not None:
        outputIDT(bestTiles, outHandle=args.idt)

    outputTable(bestTiles, outHandle=args.output)

    if args.calcPrice:
        utils.eprint(
            f"\nTotal cost to synthesize probe sets ~${calcOligoCost(bestTiles):.2f}"
        )

    outputRunParams(args)


if __name__ == "__main__":
    main()
    main()
