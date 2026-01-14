"""Sequence parsing and utility functions."""

#/usr/bin/env python
import operator,random,math
from . import prob

######
#Parsers
######
def FastaIterator(handle):
    """
    Generator function to iterate over fasta records in <handle>:
    Use in a loop to apply to each Seq record contained in a .fasta file
    Input: record handle as obtained by handle = open(<file>,'r')
    Returns an iterator across Sequences in file
    """
    #Skip any header text
    while True:
        line = handle.readline()
        if line == "" : return #Premature end of file, or just empty?
        if line [0] == ">":
            break

    while True:
        if line[0] != ">":
            raise ValueError("Records in Fasta files should start with a '>' character")
        name = line[1:].rstrip()
        lines = []
        line = handle.readline()
        while True:
            if not line : break
            if line[0] == ">" : break
            lines.append(line.rstrip().replace(" ",""))
            line = handle.readline()
        #Return record then continue
        newSeq = {'name':name,'sequence':"".join(lines)}
        yield newSeq

        if not line : return #StopIteration
    assert False, "Should not reach this line"

bed_fields = ['chr','start','end','label','score','strand']

###
#Generic Sequence tools
###

def complement(s):
    '''
    Return the complement of a DNA sequence
    
    :param s: sequence
    :return: A list of the complement bases of the sequence.
    '''
    comp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
            'a': 't', 'c': 't', 'g': 'c', 't': 'a', 'n': 'n'
            }
    complseq = [comp[base] for base in s]
    return complseq

def reverse_complement(s):
    '''
    Return the reverse complement of a DNA sequence
    
    :param s: The sequence to be reverse complemented
    :return: The reverse complement of the input sequence.
    '''
    seq = list(s)
    seq.reverse()
    return ''.join(complement(seq))

def rcomp(s):
    """Does same thing as reverse_complement only cooler"""
    return s.translate(str.maketrans("ATCG","TAGC"))[::-1]

def getTm(seq):
    '''
    The function getTm(seq) takes a sequence as an argument and returns the melting temperature of the
    sequence
    
    :param seq: the sequence of interest
    :return: The melting temperature of the sequence.
    '''
    Tm = 79.8 + 18.5*math.log10(0.05) + (58.4 * getGC(seq)) + (11.8 * getGC(seq)**2) - (820/len(seq))
    return Tm

def getGC(seq):
    '''
    The function getGC(seq) takes a string of DNA sequence as input and returns the GC content of the
    sequence
    
    :param seq: the sequence to be analyzed
    :return: A list of tuples. Each tuple contains the name of the gene, the GC content, and the length
    of the gene.
    '''
    return (seq.count('C')+seq.count('G')+seq.count('c')+seq.count('g'))/float(len(seq))

def gc_content(seq):
    '''
    Given a DNA sequence, return the percentage of G's and C's in the sequence
    
    :param seq: the sequence to be analyzed
    :return: The GC content of the sequence.
    '''
    gc = mcount(seq, 'GCgc')
    at = mcount(seq, 'ATUatu')
    return 100*gc/float((gc+at))

def mcount(s, chars):
    '''
    Sums the counts of appearances of each char in chars
    
    :param s: the string to search
    :param chars: a string of characters to count
    :return: The number of times the characters in chars appear in s.
    '''
    # sums the counts of appearances of each char in chars
    count = 0
    for char in chars:
        count = count+str.count(s,char)
    return count

def prob_seq(seq, pGC=.5):
    '''
    Given a sequence and a background GC probability, 
    what is the probability of getting that sequence
    
    :param seq: the sequence of interest
    :param pGC: the probability of a GC base pair
    :return: The probability of the sequence given the background GC probability.
    '''
    # given a GC content, what is the probability
    # of getting the particular sequence

    assert(0<=pGC<=1)
    # the probability of obtaining sequence seq
    # given a background gc probability of .5
    ps = []
    for char in seq:
        if char in 'CG': ps.append(pGC/2)
        elif char in 'AT': ps.append((1-pGC)/2)
        else: raise f"Unexpected char: {char}"
    return reduce(operator.mul, ps, 1)

def transcribe(seq):
    '''
    The function transcribe() takes a DNA sequence and replaces each instance of the nucleotide T with a
    uracil (U) in the transcribed RNA sequence
    
    :param seq: the sequence to be transcribed
    :return: The transcribed RNA sequence.
    '''
    RNA = seq.replace('T', 'U')
    return RNA

def GenRandomSeq(length, type='DNA'):
    '''
    Generate a random sequence of DNA or RNA of a given length
    
    :param length: the length of the random sequence
    :param type: DNA or RNA, defaults to DNA (optional)
    :return: A string of length length consisting of random characters from chars.
    '''
    if type == 'DNA':
        chars = ['A','T','G','C']
    if type == 'RNA':
        chars = ['A','U','G','C']
    return ''.join([random.choice(chars) for i in range(length)])

def seed():
    """Seed the random number generator with system entropy."""
    random.seed()

def draw(distribution):
    '''
    Draw a random value from the distribution, where values with a higher probability are drawn more
    often
    
    :param distribution: a list of positive numbers that sum to 1
    :return: The index of the array that the random number falls into.
    '''
    sum=0
    r = random.random()
    for i in range(0,len(distribution)):
        sum += distribution[i]
        if r< sum:
            return i

def makeDistFromFreqs(freqs):
    '''
    Given a dictionary of character frequencies, return a list of cumulative frequencies
    
    :param freqs: a dictionary of the frequencies of each nucleotide at each position
    :return: a list of cumulative frequencies.
    '''
    res = []
    chars = ['A','T','C','G']
    cum = 0
    res.append(cum)
    for i in chars:
        cum += freqs[i]
        res.append(cum)
    return res

def genRandomFromDist(length,freqs):
    """Generates a random sequence of length 'length' drawing from a distribution of
    base frequencies in a dictionary"""
    myDist = makeDistFromFreqs(freqs)
    chars = ['A','T','C','G']
    return ''.join([chars[draw(myDist)] for i in range(length)])

###########
#Motif Tools
###########
def allindices(string, sub, listindex=[], offset=0):
    '''
    Return a list of all indices of a substring in a string
    
    :param string: the string to be searched
    :param sub: The string you're looking for
    :param listindex: an empty list
    :param offset: the index in the string where you want to start searching, defaults to 0 (optional)
    :return: A list of all the indices where the substring sub is found in string.
    '''
    i = str.find(sub, offset)
    while i >= 0:
        listindex.append(i)
        i = str.find(sub, i + 1)
    return listindex

def find_all(seq, sub):
    '''
    Find all occurences of a substring in a string
    
    :param seq: the string to be searched
    :param sub: The substring to search for
    :return: A list of all the positions where the substring was found.
    '''
    #print "Looking for %s in %s"%(sub,seq)
    found = []
    next = str.find(seq,sub)
    while next != -1:
        found.append(next)
        next = str.find(seq,sub,next+1)
    return found

def kmer_dictionary_counts(seq,k,dic={}):
    """Returns a dictionary of k,v = kmer:'count of kmer in seq'"""
    for i in range(0, len(seq)-k):
        subseq = seq[i:][:k]
        #if not dic.has_key(subseq): dic[subseq] = 1
        #else: dic[subseq] = dic[subseq] + 1
        #OR
        dic[subseq] = 1 + dic.get(subseq,0)
    return dic

def kmer_dictionary(seq,k,dic={},offset=0):
    """Returns dictionary of k,v = kmer:'list of kmer start positions in seq' """
    for i in range(0,len(seq)-k):
        subseq = seq[i:][:k]
        dic.setdefault(subseq,[]).append(i+1)
    return dic

def kmer_stats(kmer,dic,genfreqs):
    """Takes as argument a kmer string, a dictionary with kmers as keys from kmer_dictionary_counts, and a dictionary
        of genomic frequencies with kmers as keys. Returns a dictionary of stats for kmer ("Signal2Noise Ratio, Z-score")
    """
    if not dic: return
    if kmer in dic.keys() and kmer in genfreqs.keys():
        observed = dic[kmer]
        expected = sum(dic.values())*genfreqs[kmer]
        snr = prob.snr(observed,expected)
        zscore = prob.zscore(observed, expected)
        return {'snr':snr,'zscore':zscore}
    else: return

def get_seeds(iter,seeds={}):
    '''
    Given a list of sequences, return a dictionary of the counts of each seed
    
    :param iter: the iterator of sequences
    :param seeds: a dictionary of seeds and their counts
    :return: A dictionary with the seeds as keys and the number of times they occur as values.
    '''
    counter = 0
    for i in iter:
        counter+=1
        if counter%10000==0:
            print(f"{counter}")
        i.CSToDNA()
        seed = i.sequence[1:8]
        seeds[seed] = 1 + seeds.get(seed,0)
    return seeds
