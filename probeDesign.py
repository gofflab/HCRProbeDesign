from probeDesign.tiles import Tile, TileError
from probeDesign import utils
#import copy
#import string
from probeDesign.utils import pp
from probeDesign import sequencelib
from probeDesign import repeatMask
from probeDesign import BLAST
import getopt,sys,re
#from Bio.Seq import Seq
import primer3
from string import ascii_uppercase

#######################
# Scan input sequence #
#######################
#TODO: Modify this so that it only gets the window of appropriate size.  We will add prefix and suffix afterwards.

def scanSequence(sequence,seqName,tileStep=1,tileSize=52):
	tiles = []
	#Pre-compute number of chunks to emit
	numOfChunks = int(((len(sequence)-tileSize)/tileStep) + 1)

	#Tile across reverse complement of sequence
	for i in range(0,numOfChunks*tileStep,tileStep):
		tile = Tile(sequence=sequencelib.reverse_complement(sequence[i:i+tileSize]),seqName=seqName,startPos=i+1)
		if not tile.isMasked():
			tiles.append(tile)
	return tiles

###################
# Reporting
###################

def outputTable(tiles,outHandle=sys.stdout):
	outputKeys=["name","seqName","startPos","oligoSequence","prefix","suffix","tag","GC"]
	print >>outHandle, "\t".join(outputKeys)

	for tile in tiles:
		vals = []
		for k in outputKeys:
			v = getattr(tile,k)
			if callable(v):
				v = v()
			vals.append(str(v))
		print >>outHandle, "\t".join(vals)

def usage():
	utils.eprint('Help Message Goes Here')

#def main():
	# #######################
	# # Variables
	# #######################
    #
	# maxTileSize = 52
	# tileStep = 1
    #
	# #Argument handling
	# try:
	# 	opts,args = getopt.getopt(sys.argv[1:],"ho:vp:s:l:w:te:n:r:",["help","output=","verbose","prefix=","suffix=","max-tile-length=","tile-step=","tag","tag-length=","num-tags-per-tile=","restriction-sites="])
	# except getopt.GetoptError as err:
	# 	print(err)
	# 	usage()
	# 	sys.exit(2)
	# output = None
	# verbose = False
	# for o,a in opts:
	# 	if o == "-v":
	# 		verbose = True
	# 	elif o in ("-h","--help"):
	# 		usage()
	# 		sys.exit()
	# 	elif o in ("-o","--output"):
	# 		output = a
	# 	elif o in ("-p","--prefix"):
	# 		if utils.onlyNucleic(a):
	# 			prefix = a
	# 		else:
	# 			raise(TileError,"Prefix is not a valid nucleic acid sequence.")
	# 	elif o in ("-s","--suffix"):
	# 		if utils.onlyNucleic(a):
	# 			suffix = a
	# 		else:
	# 			raise(TileError,"Suffix is not a valid nucleic acid sequence.")
	# 	elif o in ("-l","--max-tile-length"):
	# 		maxTileSize = int(a)
	# 	elif o in ("-w","--tile-step"):
	# 		tileStep = int(a)
	# 	elif o in ("-t","--tag"):
	# 		addTag = True
	# 	elif o in ("-e","--tag-length"):
	# 		tagLength = int(a)
	# 	elif o in ("-n","--num-tags-per-tile"):
	# 		numTagsPerTile = int(a)
	# 	elif o in ("-r","--restriction-sites"):
	# 		sites = a.strip()
	# 		if sites == "":
	# 			sites = None
	# 	else:
	# 		assert False, "Unhandled option"
	# # Grab fname as remainder argument
	# try:
	# 	fname = str(args[0])
	# 	handle = open(fname,'r')
	# except:
	# 	usage()
	# 	sys.exit(2)
    #
    #
	# #Find window size
	# tileSize = maxTileSize
    #
    # #Fetch all tile sequences
	# fastaIter = sequencelib.FastaIterator(handle)
	# tiles = []
	# for mySeq in fastaIter:
	# 	#Warn about masked regions
	# 	if sites != None:
	# 		warnRestrictionSites(mySeq['sequence'],mySeq['name'],sites)
    #
	# 	#Get tiles from sequence
	# 	tmpTiles = scanSequence(mySeq['sequence'],mySeq['name'],tileStep=tileStep,tileSize=tileSize)
	# 	tiles += tmpTiles
    #
	# # Remove duplicate tile sequences
	# tiles = findUnique(tiles)
    #
	# # Check tiles for restriction sites
	# if sites != None:
	# 	cleanTiles = set()
	# 	for tile in tiles:
	# 		if hasRestrictionSites(tile.sequence, sites):
	# 			continue
	# 		else:
	# 			cleanTiles.add(tile)
	# 	tiles = list(cleanTiles)
    #
	# #determine number of tags needed
	# numTagsReq = len(tiles) * numTagsPerTile
    #
	# tags = buildTags(numTagsReq,tagLength,sites=sites)
    #
	# assert len(tags) == len(tiles) * numTagsPerTile
    #
	# #Create numTagsPerTile tiles for each sequence
	# tmpTiles = set()
	# #Add prefix, suffix, and tag
	# for i in xrange(len(tiles)):
	# 	for j in xrange(numTagsPerTile):
	# 		tmpTile = copy.copy(tiles[i])
	# 		tmpTile.name = "%s:%d" % (tmpTile.name,j)
	# 		tmpTile.prefix = prefix
	# 		tmpTile.suffix = suffix
	# 		tmpTile.tag = tags.pop()
	# 		tmpTiles.add(tmpTile)
    #
	# tiles = list(tmpTiles)
	# tiles.sort()
    #
	# #Just for QC
	# #for i in xrange(10):
	# outputTable(tiles)
    #
	# print >>sys.stderr, "There are a total of %d unique tiles" % len(tiles)

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
	species = 'mouse'
	genomemask = False

	# Read custom args

	# Parse fasta files and loop over records
	utils.eprint("Reading in Fasta file")
	#fname = "test/eGFP.fa"
	fname = "test/mKcnip3.fa"
	handle = open(fname,'r')
	fastaIter = sequencelib.FastaIterator(handle)

	mySeq = next(fastaIter)

	# RepeatMasking
	#TODO: Specify DNA source in input params
	utils.eprint("\nRepeat Masking...")
	mySeq['sequence'] = repeatMask.repeatmask(mySeq['sequence'],dnasource=species)

	#Convert to lowercase
	#mySeq['sequence'] = mySeq['sequence'].lower()

	# Tile over masked sequence record to generate all possible probes of appropriate length that are not already masked
	tiles = scanSequence(mySeq['sequence'],mySeq['name'],tileStep=1,tileSize=tileSize)
	utils.eprint(f'{len(tiles)} tiles available of length {tileSize}...')

	#utils.pp(tiles[0])

	# Check for invalid characters

	# Crunmask
	utils.eprint("\nChecking for runs of C's")
	tiles = [tile for tile in tiles if not tile.hasRuns(runChar='c',runLength=7,mismatches=2)]
	utils.eprint(f'{len(tiles)} tiles remain')

	# Grunmask
	utils.eprint("\nChecking for runs of G's")
	tiles = [tile for tile in tiles if not tile.hasRuns(runChar='g',runLength=7,mismatches=2)]
	utils.eprint(f'{len(tiles)} tiles remain')

	# Calculate Hairpins
	utils.eprint("\nChecking for hairpins")
	for tile in tiles:
		thermRes = primer3.calcHairpin(tile.sequence)
	tiles = [tile for tile in tiles if primer3.calcHairpin(tile.sequence).structure_found]
	utils.eprint(f'{len(tiles)} tiles remain')

	print(len(tiles))

	# PseudogeneMasking.  BLAST?

	# GenomeMasking?  BLAST instead of bowtie
	if genomemask:
		utils.eprint(f"\nBLASTN on remaining tiles against {species} reference database")
		blast_string = "\n".join([tile.toFasta() for tile in tiles[:10]])
		# if verbose:
		# 	utils.eprint("Submitting BLAST request with the following string")
		# 	utils.eprint(blast_string)
		blast_res = BLAST.blastProbes(blast_string, species=species)
		utils.eprint(f'Parsing BLAST output now')
		BLAST.getNHits(blast_res)



	# GC filtering

	# TM filtering
	utils.eprint(f"\nChecking for {minGC} < GC < {maxGC}")
	tiles = [tile for tile in tiles if tile.GC() >= minGC]
	tiles = [tile for tile in tiles if tile.GC() <= maxGC]
	utils.eprint(f'{len(tiles)} tiles remain')

	# Gibbs filtering
	utils.eprint(f"\nChecking for {minGibbs} < Gibbs FE < {maxGibbs}")
	[tile.calcGibbs() for tile in tiles]
	tiles = [tile for tile in tiles if tile.Gibbs >= minGibbs]
	tiles = [tile for tile in tiles if tile.Gibbs <= maxGibbs]
	utils.eprint(f'{len(tiles)} tiles remain')

	# dTm between halves
	utils.eprint(f"\nChecking for dTm <= {dTmMax}")
	[tile.splitProbe() for tile in tiles]
	[tile.calcdTm() for tile in tiles]
	tiles = [tile for tile in tiles if tile.dTm <=dTmMax]
	utils.eprint(f'{len(tiles)} tiles remain')

	# Break remaining probes into non-overlapping regions
	regions = {}
	regionCount = 0

	regionList = [tiles[0]]
	for i in range(1,len(tiles)):
		if i == len(tiles)-1:
			regions[ascii_uppercase[regionCount]] = regionList
			break
		if tiles[i].overlaps(tiles[i-1]):
			regionList.append(tiles[i])
		else:
			regions[ascii_uppercase[regionCount]] = regionList
			regionList = [tiles[i]]
			regionCount = regionCount + 1

	# Select best from region
	hits = {}
	for k,v in regions.items():
		hits[k] = v[min(range(len(v)), key=lambda i: abs([x.Gibbs for x in v][i]-targetGibbs))]

	for k,tile in hits.items():
		print(f'Non-overlapping region {k}')
		#for tile in v:
		print(f"{tile}\tLength:{len(tile)}\tmyTm:{tile.Tm():.2f}\tprimer3-Tm:{primer3.calcTm(tile.sequence):.2f}\tdTm:{tile.dTm:.2f}\tGC%:{tile.GC():.2f}\tGibbs:{tile.Gibbs:.2f}")

	utils.eprint(f'{len(tiles)} tiles total')

# def findProbes(infile, outfile, nProbes, species='mouse',
#                 repeatmask=True,
#                 pseudogenemask=True,
#                 blastmask=True,
#                 miRNAmask=True,
#                 humanfiltermask=True,
#                 genomemask=True,
#                 GCrunmask=True,
#                 GCmask=True,
#                 targetTM = 67.5,
#                 oligoLength = 20,
#                 spacerLength=2,
#                 masksequences=Null,
#                 targetGibbsFE = -23,
#
#  ):
	# # Set default args
	# species = "mouse"
	# dTmMax = 5.0
	# minGibbs = -80.0
	# maxGibbs = -50.0
	# targetGibbs = 60.0
	# tileSize = 52
	# #Argument handling
	# try:
	# 	opts,args = getopt.getopt(sys.argv[1:],"ho:vs:l:w:tfrpbg",["help","output=","verbose","prefix=","species=","oligoLength=","spacerLength=","targetTM=","targetGibbsFE=","repeatmask","pseudogenemask","blastmask","genomemask"])
	# except getopt.GetoptError as err:
	# 	print(err)
	# 	usage()
	# 	sys.exit(2)


if __name__ == "__main__":
    test()
