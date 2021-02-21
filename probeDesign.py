from probeDesign import utils
import copy
#import string
from probeDesign.utils import pp
from probeDesign import sequencelib
from probeDesign import thermo
from probeDesign import repeatMask
import getopt,sys,re
from Bio.Seq import Seq
import primer3

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


class TileError(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class Tile:
	def __init__(self,sequence,seqName,startPos,prefix='',suffix='',tag=''):
		self.sequence = str.lower(sequence)
		self.startPos = startPos
		self.start = startPos
		self.seqName = seqName
		self.name = "%s:%d-%d" % (self.seqName,self.start,self.start+len(self.sequence))
		self.prefix = prefix
		self.suffix = suffix
		self.tag = tag
		self.masked = False
		#self.RajTM = self.calcRajTm()


	def validate(self):
		self.getGC()

	def compiledPrefix(self):
		""" Check prefix for '@' indicating position to add tag'"""
		tagPos = self.prefix.find('@')
		if tagPos == -1:
			return self.prefix
		else:
			return self.prefix[:tagPos]+self.tag+self.prefix[tagPos:]

	def compiledSuffix(self):
		""" Check suffix for '@' indicating position to add tag'"""
		tagPos = self.suffix.find('@')
		if tagPos == -1:
			return self.suffix
		else:
			return self.suffix[:tagPos]+self.tag+self.suffix[tagPos+1:]

	def __repr__(self):
		return "%s:%s" % (self.name,self.oligoSequence())

	def __str__(self):
		#return "%s\t%0.2f\t%d" % (self.__repr__(),self.GC,len(self))
		return "%s" % (self.__repr__())

	def toFasta(self):
		return ">%s\n%s" % (self.name,self.sequence)

	def toBed(self):
		pass

	def GC(self):
		return float(sequencelib.gc_content(self.oligoSequence()))

	def oligoSequence(self):
		return self.compiledPrefix()+self.sequence+self.compiledSuffix()

	def __hash__(self):
		return hash(self.oligoSequence())

	def __eq__(self,other):
		#if self.sequence.upper() == other.sequence.upper():
		if self.oligoSequence() == other.oligoSequence():
			return True
		else:
			return False

	def __len__(self):
		return len(self.oligoSequence())

	def __cmp__(self,other):
		return cmp((self.seqName, self.startPos, self.name),(other.seqName, other.startPos, other.name))

	def toFasta(self):
		return ">%s\n%s" % (self.name,self.oligoSequence())

	def tileFasta(self):
		"""Only write tile sequence to fasta"""
		return ">%s\n%s" % (self.name,self.sequence)

	def calcGibbs(self):
		[dHs,dSs] = thermo.stacks_rna_dna(self.sequence)
		[dHi,dSi] = thermo.init_rna_dna()
		binding_energy = thermo.gibbs(dHs+dHi,dSs+dSi,temp=37)  # cal/mol
		binding_energy = thermo.salt_adjust(binding_energy/1000,len(self.sequence),saltconc=0.33)  # kcal/mol
		self.Gibbs = binding_energy

	def Tm(self):
		return float(sequencelib.getTm(self.sequence))

	def calcRajTm(self):
		[dHs,dSs] = thermo.stacks_rna_dna(self.sequence)
		#rajTm = thermo.melting_temp(dHs,dSs,)
		pass

	def isMasked(self):
		if 'n' in self.sequence:
			self.masked = True
		elif 'N' in self.sequence:
			self.masked = True
		return self.masked

	def hasRuns(self,runChar,runLength,mismatches):
		answer = False
		for i in range(len(self)-runLength+1):
			count = 0
			for j in range(i,i+runLength):
				if self.sequence[j] == runChar:
					count += 1
			if count >= runLength-mismatches:
				self.masked = True
				answer = True
		return answer

	def splitProbe(self):
		self.oddSeq = self.sequence[:int(len(self)/2)]
		self.evenSeq = self.sequence[int(len(self)/2):]
		return

	def calcdTm(self):
		self.dTm = abs(primer3.calcTm(self.oddSeq)-primer3.calcTm(self.evenSeq))


#######################
# Scan input sequence #
#######################
#TODO: Modify this so that it only gets the window of appropriate size.  We will add prefix and suffix afterwards.

def scanSequence(sequence,seqName,tileStep=1,tileSize=52):
	tiles = []
	#Pre-compute number of chunks to emit
	numOfChunks = int(((len(sequence)-tileSize)/tileStep) + 1)

	print(numOfChunks)
	#Tile across sequence
	for i in range(0,numOfChunks*tileStep,tileStep):
		tile = Tile(sequence=sequence[i:i+tileSize],seqName=seqName,startPos=i+1)
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
	dTmMax = 5.0
	minGibbs = -80.0
	maxGibbs = -50.0
	targetGibbs = 60.0
	tileSize = 52

	# Read custom args

	# Parse fasta files and loop over records
	utils.eprint("Reading in Fasta file")
	fname = "test/eGFP.fa"
	handle = open(fname,'r')
	fastaIter = sequencelib.FastaIterator(handle)

	mySeq = next(fastaIter)

	# RepeatMasking
	#TODO: Specify DNA source in input params
	utils.eprint("\nRepeat Masking...")
	mySeq['sequence'] = repeatMask.repeatmask(mySeq['sequence'],dnasource='mouse')

	# Tile over masked sequence record to generate all possible probes of appropriate length that are not already masked
	tiles = scanSequence(mySeq['sequence'],mySeq['name'],tileStep=1,tileSize=tileSize)

	# Check for invalid characters

	# Crunmask
	utils.eprint("\nChecking for runs of C's")
	tiles = [tile for tile in tiles if tile.hasRuns(runChar='c',runLength=7,mismatches=2)]
	utils.eprint(f'{len(tiles)} tiles remain')

	# Grunmask
	utils.eprint("\nChecking for runs of G's")
	tiles = [tile for tile in tiles if tile.hasRuns(runChar='g',runLength=7,mismatches=2)]
	utils.eprint(f'{len(tiles)} tiles remain')

	# Calculate Hairpins
	utils.eprint("\nChecking for hairpins")
	for tile in tiles:
		thermRes = primer3.calcHairpin(tile.sequence)
	tiles = [tile for tile in tiles if primer3.calcHairpin(tile.sequence).structure_found]
	utils.eprint(f'{len(tiles)} tiles remain')

	print(len(tiles))

	# PseudogeneMasking


	# GenomeMasking?  Bowtie hits to >1 regions of the genomemask

	# GC filtering

	# TM filtering

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

	# Select top 20

	for tile in tiles:
		print(f"{tile}\tLength:{len(tile)}\tTm:{tile.Tm():.2f}\tprimer3-Tm:{primer3.calcTm(tile.sequence):.2f}\tdTm:{tile.dTm:.2f}\tGC%:{tile.GC():.2f}\tGibbs:{tile.Gibbs:.2f}")

	print(len(tiles))

if __name__ == "__main__":
    test()
