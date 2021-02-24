from . import utils
from . import thermo
from . import sequencelib
import primer3


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
		self.end = startPos + len(self.sequence)
		self.seqName = seqName
		self.name = "%s:%d-%d" % (self.seqName,self.start,self.start+len(self.sequence))
		self.prefix = prefix
		self.suffix = suffix
		self.tag = tag
		self.masked = False
		#self.RajTM = self.calcRajTm()


	def validate(self):
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
		return f"{self.name}:{self.sequence}"

	def __str__(self):
		#return "%s\t%0.2f\t%d" % (self.__repr__(),self.GC,len(self))
		return f"{self.__repr__()}"

	def __iter__(self):
		return iter(self.sequence)

	def __len__(self):
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
		return f'>{self.name}\n{self.sequence}'

	def toBed(self):
		pass

	def GC(self):
		return float(sequencelib.gc_content(self.sequence))

	#def oligoSequence(self):
	#	return self.compiledPrefix()+self.sequence+self.compiledSuffix()

	def __hash__(self):
		return hash(self.sequence)

	def __eq__(self,other):
		#if self.sequence.upper() == other.sequence.upper():
		if self.sequence == other.sequence:
			return True
		else:
			return False

	def __len__(self):
		return len(self.sequence)

	def __cmp__(self,other):
		return cmp((self.seqName, self.startPos, self.name),(other.seqName, other.startPos, other.name))

	def toFasta(self):
		return ">%s\n%s" % (self.name,self.sequence)

	# def tileFasta(self):
	# 	"""Only write tile sequence to fasta"""
	# 	return ">%s\n%s" % (self.name,self.sequence)

	def calcGibbs(self):
		[dHs,dSs] = thermo.stacks_rna_dna(self.sequence)
		[dHi,dSi] = thermo.init_rna_dna()
		binding_energy = thermo.gibbs(dHs+dHi,dSs+dSi,temp=37)  # cal/mol
		binding_energy = thermo.salt_adjust(binding_energy/1000,len(self.sequence),saltconc=0.33)  # kcal/mol
		self.Gibbs = binding_energy

	def Tm(self):
		return float(sequencelib.getTm(self.sequence))

	def RajTm(self):
		return thermo.Tm(self.sequence)

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