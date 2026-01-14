"""Tile and probe representation used in the design pipeline."""

from . import utils
from . import thermo
from . import sequencelib
import primer3
from . import HCR


# This class is used to raise exceptions.
class TileError(Exception):
	"""Custom exception type for tile validation and processing."""
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class Tile:
	"""Represents a candidate probe tile extracted from a target sequence."""
	def __init__(self,sequence,seqName,startPos):
		"""
		Initialize a Tile from a sequence and positional metadata.

		:param sequence: Tile sequence (string).
		:param seqName: Source sequence name.
		:param startPos: 1-based start position in the source sequence.
		"""
		self.sequence = str.lower(sequence)
		self.startPos = startPos
		self.start = startPos
		self.end = startPos + len(self.sequence)
		self.seqName = seqName
		self.name = f"{self.seqName}:{self.start}-{self.start+len(self.sequence)}".replace(" ", "_")
		self.masked = False
		self.hitCount = -1 #-1 indicates that genome masking has not yet been performed.
		#self.RajTM = self.calcRajTm()


	def validate(self):
		"""Run lightweight validation checks on the tile."""
		self.GC()

	# def compiledPrefix(self):
	# 	""" Check prefix for '@' indicating position to add tag'"""
	# 	tagPos = self.prefix.find('@')
	# 	if tagPos == -1:
	# 		return self.prefix
	# 	else:
	# 		return self.prefix[:tagPos]+self.tag+self.prefix[tagPos:]

	# def compiledSuffix(self):
	# 	""" Check suffix for '@' indicating position to add tag'"""
	# 	tagPos = self.suffix.find('@')
	# 	if tagPos == -1:
	# 		return self.suffix
	# 	else:
	# 		return self.suffix[:tagPos]+self.tag+self.suffix[tagPos+1:]

	def __repr__(self):
		"""Return a debug-friendly representation of the tile."""
		return f"{self.name}:{self.sequence}"

	def __str__(self):
		"""Return a human-readable representation of the tile."""
		#return "%s\t%0.2f\t%d" % (self.__repr__(),self.GC,len(self))
		return f"{self.__repr__()}"

	def __iter__(self):
		"""Iterate over the tile sequence bases."""
		return iter(self.sequence)

	def __len__(self):
		"""Return the length of the tile interval."""
		return self.end-self.start+1

	def overlaps(self,b):
			"""Return true if b overlaps self"""
			if (self.start <= b.start and b.start <=self.end) or (self.start >= b.start and self.start <= b.end):
				return True
			else:
				return False

	def distance(self,b,enforceStrand=False):
		"""
		Returns absolute distance between self and another interval start positions.
		"""
		return abs(self.start-b.start)

	def toFasta(self):
		"""Return the tile formatted as a FASTA record."""
		return f'>{self.name}\n{self.sequence}'

	def toBed(self):
		"""Placeholder for BED formatting support."""
		pass

	def GC(self):
		"""Return GC percentage for the tile sequence."""
		return float(sequencelib.gc_content(self.sequence))

	#def oligoSequence(self):
	#	return self.compiledPrefix()+self.sequence+self.compiledSuffix()

	def __hash__(self):
		"""Hash tiles by their sequence."""
		return hash(self.sequence)

	def __eq__(self,other):
		"""Compare tiles by their sequences."""
		#if self.sequence.upper() == other.sequence.upper():
		if self.sequence == other.sequence:
			return True
		else:
			return False

	def __len__(self):
		"""Return the length of the tile sequence."""
		return len(self.sequence)

	def __cmp__(self,other):
		"""Legacy comparison for sorting tiles by name and position."""
		return cmp((self.seqName, self.startPos, self.name),(other.seqName, other.startPos, other.name))

	# def tileFasta(self):
	# 	"""Only write tile sequence to fasta"""
	# 	return ">%s\n%s" % (self.name,self.sequence)

	def calcGibbs(self):
		'''
		Calculate the Gibbs free energy of binding for a given sequence
		'''
		[dHs,dSs] = thermo.stacks_rna_dna(self.sequence)
		[dHi,dSi] = thermo.init_rna_dna()
		binding_energy = thermo.gibbs(dHs+dHi,dSs+dSi,temp=37)  # cal/mol
		binding_energy = thermo.salt_adjust(binding_energy/1000,len(self.sequence),saltconc=0.33)  # kcal/mol
		self.Gibbs = binding_energy

	def Tm(self):
		"""Return the basic melting temperature estimate for the tile."""
		return float(sequencelib.getTm(self.sequence))

	def RajTm(self):
		"""Return the SantaLucia-style melting temperature estimate."""
		return thermo.Tm(self.sequence)

	def isMasked(self):
		"""Return True if the tile contains masked bases."""
		if 'n' in self.sequence:
			self.masked = True
		elif 'N' in self.sequence:
			self.masked = True
		return self.masked

	def hasRuns(self,runChar,runLength,mismatches):
		'''
		Given a sequence, a run character, a run length, and a number of mismatches, 
		returns True if the sequence has a run of the specified character of the specified length, 
		with the specified number of mismatches
		
		:param runChar: the character that indicates a run
		:param runLength: the length of the run of the same character
		:param mismatches: the number of mismatches allowed in the run
		:return: A boolean value.
		'''
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
		"""
		Split sequence in half with two bases in the middle removed (flexible gap to help initiator sequence land)
		ie. a 52mer will be split into two 25mers with the middle two bases of the 52mer dropped
		"""
		self.fivePrimeSeq = self.sequence[:int(len(self)/2)-1]
		self.threePrimeSeq = self.sequence[int(len(self)/2)+1:]
		return

	def calcdTm(self):
		'''
		Calculate the difference in melting temperature between the 5' and 3' sequences
		'''
		self.dTm = abs(primer3.calc_tm(self.fivePrimeSeq)-primer3.calc_tm(self.threePrimeSeq))

	#TODO: PLEASE check this to make sure that I'm adding the initiator sequences in the correct position and order
	def makeProbes(self,channel):
		'''
		This function creates the probes for the channel.
		
		:param channel: the channel that the probe is on
		'''
		self.P1 = HCR.initiators[channel]["odd"]+self.threePrimeSeq
		self.P2 = self.fivePrimeSeq + HCR.initiators[channel]["even"]
		self.channel = channel
