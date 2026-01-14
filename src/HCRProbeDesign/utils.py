"""Miscellaneous utilities for sequence processing and formatting."""

import string
import operator
import random
import math
from Bio import Restriction
import sys
import types

##############
#FastaIterator
##############
def FastaIterator(handle):
    """
    Generator function to iterate over fasta records in <handle>:
    Use in a loop to apply to each sequence record contained in a .fasta file
    Input: record handle as obtained by handle = open(<file>,'r')
    Returns an iterator across sequences in file
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


def eprint(*args, **kwargs):
    """
    Print to stderr with the same signature as print().

    :return: None.
    """
    print(*args, file=sys.stderr, **kwargs)

########
#
#Pretty Printing
#
########
def pretty_print(f, d, level=-1, maxw=0, maxh=0, gap="", first_gap='', last_gap=''):
    """
    Pretty-print nested structures to a file handle.

    :param f: Output file-like object.
    :param d: Data structure to render.
    :param level: Recursion depth (-1 for unlimited).
    :param maxw: Maximum width per line.
    :param maxh: Maximum items per list/dict/tuple.
    :param gap: Base indentation string.
    :param first_gap: Indentation for opening delimiter line.
    :param last_gap: Indentation for closing delimiter line.
    :return: None.
    """
    # depending on the type of expression, it recurses through its elements
    # and prints with appropriate indentation

    # f   is the output file stream
    # d   is the data structure
    #
    # level is the number of allowed recursive calls, the depth at which
    #       the data structure is explored
    #       default: -1 means never stop recursing early
    # maxw  is the maximum width that will be printed from the last element
    #       of the recursion (when no further recursion is possible, or
    #       the maximal depth has been reached)
    #       default: 0 means every line will be printed in its entirety, regardless
    #                of how long it may be
    # maxh  (max height) is the maximum number of elements that will be
    #       printed from a list or a dictionary, at any level or recursion
    #       default: 0 means every list or dictionary will have all its elements
    #                printed, even if it contains thousands of elements
    #
    # gap is the gap to include before every element of a list/dic/tuple
    # first_gap is the opening gap before the opening bracket, parens or curly braces
    # first_gap is the closing gap before the closing bracket, parens or curly braces

    if level == 0:
        if type(d) != str: d = repr(d)

        if maxw and len(d) > maxw:
            final = ifab(maxw > 20, 10, maxw/2)
            f.write(first_gap+d[:maxw-final]+'...'+d[-final:]+' (%s chars)\n' % len(d))
        else: f.write(first_gap+d+'\n')
    elif type(d) == list:
        if not d:
            f.write(first_gap+"[]\n")
            return
        # recurse on lists
        f.write(first_gap+"[\n")
        h = 0
        for el in d:
            pretty_print(f, el, level-1, maxw, maxh, gap+'   ', gap+' ->', gap+'   ')
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(d):
                    f.write(gap+' -> ... (%s in list)\n'%len(d))
                    break
        f.write(last_gap+"]\n")
    elif type(d) == tuple:
        if not d:
            f.write(first_gap+"()\n")
            return
        # recurse on tuples
        f.write(first_gap+"(\n")
        h = 0
        for el in d:
            pretty_print(f, el,
                         level     = level-1,
                         maxw      = maxw,
                         maxh      = maxh,
                         gap       = gap+'   ',
                         first_gap = gap+' =>',
                         last_gap  = gap+'   ')
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(d):
                    f.write(gap+' => ... (%s in tuple)\n'%len(d))
                    break
        f.write(last_gap+")\n")
    elif type(d) == dict:
        if not d:
            f.write(first_gap+"{}\n")
            return
        # recurse on dictionaries
        f.write(first_gap+"{\n")
        keys = d.keys()
        keys.sort()
        key_strings = map(lambda k: ifab(type(k)==str, k, repr(k)), keys)
        maxlen = max(map(len, key_strings))
        h = 0
        for k,key_string in map(None, keys, key_strings):
            key_string = sfill(key_string,maxlen,'.')
            blank_string = ' '*len(key_string)
            pretty_print(f, d[k],
                         level     = level-1,
                         maxw      = maxw,
                         maxh      = maxh,
                         gap       = gap+'    %s'%blank_string,
                         first_gap = gap+'  %s: '%key_string,
                         last_gap  = gap+'    %s'%blank_string)
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(keys):
                    remaining_keys = []
                    for k in keys[h:]:
                        if type(k) == tuple:
                            remaining_keys.append(repr(k))
                        else:
                            remaining_keys.append('%s'%k)
                    remaining_keys = string.join(remaining_keys,',')
                    #f.write(gap+'  %s (%s keys)\n'%(remaining_keys, len(keys)))
                    pretty_print(f, '  %s (%s keys)'%(remaining_keys, len(keys)),0,maxw,0,
                                 gap,gap,'')
                    break

            #gap+' '*(len(key_string)+3), '', gap+' '*(len(key_string)+5))
        f.write(last_gap+"}\n")
    elif type(d) == object:
        fields = dir(d)

        if not fields:
            f.write(first_gap+"*EmptyClass*\n")
            return
        # recurse on classes
        f.write(first_gap+"*ClassInstance %s\n"%d)
        fields.sort()
        key_strings = map(lambda k: ifab(type(k)== str, k, repr(k)), fields)
        maxlen = max(map(len, key_strings))
        h = 0
        for k,key_string in map(None, fields, key_strings):
            key_string = sfill(key_string,maxlen,'.')
            blank_string = ' '*len(key_string)
            pretty_print(f, eval('d.'+k),
                         level     = level-1,
                         maxw      = maxw,
                         maxh      = maxh,
                         gap       = gap+'    %s'%blank_string,
                         first_gap = gap+'  %s: '%key_string,
                         last_gap  = gap+'    %s'%blank_string)
            if maxh:
                h = h+1
                if h >= maxh and maxh<len(keys):
                    remaining_keys = []
                    for k in keys[h:]:
                        if type(k) == type(()):
                            remaining_keys.append(repr(k))
                        else:
                            remaining_keys.append('%s'%k)
                    remaining_keys = string.join(remaining_keys,',')
                    #f.write(gap+'  %s (%s keys)\n'%(remaining_keys, len(keys)))
                    pretty_print(f,
                                 '  %s (%s keys)'%(remaining_keys, len(keys)),
                                 0,
                                 maxw,
                                 0,
                                 gap,
                                 gap,
                                 '')
                    break

            #gap+' '*(len(key_string)+3), '', gap+' '*(len(key_string)+5))
        f.write(last_gap+"*\n")
    elif type(d) == type(""):
        # simply print strings (no quotes)
        if maxw and len(d)>maxw:
            final = ifab(maxw > 20, 10, maxw/2)
            f.write(first_gap+d[:maxw-final]+'..'+d[-final:]+' (%s)\n' % len(d))
        else:
            f.write(first_gap+d+'\n')
    else:
        # string conversion of all other types
        if maxw and len(repr(d))>maxw:
            final = ifab(maxw > 20, 10, maxw/2)
            f.write(first_gap+repr(d)[:maxw-final]+'..'+repr(d)[-final:]+' (%s)\n' % len(repr(d)))
        else:
            f.write(first_gap+repr(d)+'\n')

def pp(d,level=-1,maxw=0,maxh=0,parsable=0):
    """ wrapper around pretty_print that prints to stdout"""
    if not parsable:
        pretty_print(sys.stderr, d, level, maxw, maxh, '', '', '')
    else:
        import pprint
        if maxw: pp2 = pprint.PrettyPrinter(width=maxw, indent=1)#, depth=level
        else: pp2 = pprint.PrettyPrinter(indent=1)#, depth=level
        pp2.pprint(d)

def onlyNucleic(seq,set=['a','c','g','t','u','A','C','G','T','U','n','N','@']):
	"""
	Check whether a sequence contains only nucleic characters.

	:param seq: Input sequence string.
	:param set: Allowed characters.
	:return: True if all characters are allowed.
	"""
	for c in seq:
		if c not in set:
			return False
	return True

def findUnique(tiles):
	"""Return a list of unique Tile objects from the input list."""
	return list(set(tiles))

# def findUnique(tiles):
# 	seen = set()
# 	res = []
# 	for tile in tiles:
# 		if tile not in seen:
# 			seen.add(tile)
# 			res.append(tile)
# 	return res

def estimateAffixLength(sequence,tagLength):
	"""
	Estimate sequence length after tag insertion.

	:param sequence: Sequence containing an optional '@' tag marker.
	:param tagLength: Length of the tag to be inserted.
	:return: Adjusted sequence length.
	:raises TileError: If multiple tag markers are present.
	"""
	tagHits = sequencelib.mcount(sequence, '@')
	if tagHits == 0:
		return len(sequence)
	elif tagHits > 1:
		raise TileError("""You can only have one instance of 'tag' per tile""")
	elif tagHits == 1:
		return len(sequence) + tagLength - 1 #-1 is required upone removal of '@' tag


def buildTags(numTags,tagLength,sites=None):
	"""
	Generate random DNA tags with optional restriction site filtering.

	:param numTags: Number of tags to generate.
	:param tagLength: Length of each tag.
	:param sites: Comma-delimited restriction sites to avoid.
	:return: List of DNA tag strings.
	"""
	tmpTags = set()
	while len(tmpTags)<numTags:
		tmpTag = sequencelib.GenRandomSeq(tagLength,type="DNA")
		if sites != None:
			if hasRestrictionSites(tmpTag,sites):
				continue
		tmpTags.add(tmpTag)
	return list(tmpTags)

def hasRestrictionSites(sequence,sites):
	"""
	Check if a sequence contains restriction sites.

	:param sequence: Sequence to scan.
	:param sites: Comma-delimited restriction site names.
	:return: True if any sites are present.
	"""
	#Parse sites
	sites = sites.split(",")
	rb = Restriction.RestrictionBatch(sites)

	#Get Bio.Seq object
	#amb = IUPACAmbiguousDNA()
	tmpSeq = Seq(sequence)

	#Search for sites
	res = rb.search(tmpSeq)

	#Sum hits
	totalSites = 0
	for v in res.values():
		totalSites += len(v)

	if totalSites > 0:
		return True
	else:
		return False

def warnRestrictionSites(sequence,name,sites):
	"""
	Print a warning if restriction sites are found in a sequence.

	:param sequence: Sequence to scan.
	:param name: Sequence label for logging.
	:param sites: Comma-delimited restriction site names.
	:return: None.
	"""
	sites = sites.split(",")
	rb = Restriction.RestrictionBatch(sites)

	#Get Bio.Seq object
	#amb = IUPACAmbiguousDNA()
	tmpSeq = Seq(sequence)

	#Search for sites
	res = rb.search(tmpSeq)

	#Sum hits
	totalSites = 0
	for v in res.values():
		totalSites += len(v)

	if totalSites > 0:
		print >>sys.stderr, "Warning: The following positions in '%s' will be masked from tiles due to incompatible restictions sites:" % (name)
		pp(res)
	else:
		pass
