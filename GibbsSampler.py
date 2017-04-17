"""
This module implements the gibbs sampler algorithm, which is used to find common motifs in DNA sequences.

GibbsSampler.py
  
Created by Scott Czopek on 1/22/17.
Copyright (c) 2017. All rights reserved.


python 2.7.12
"""
# Need some random number functionality and the operator.add ability.
import random
import operator

def multipleSeedsGibbsSampling(dna, numSeeds, k, N):

	"""A motif finding algorithm that finds one common motif and returns a list bestMotifs containing the closest motif match from each string in "dna".  The library was written to support this function.

	Keyword arguments:
	"dna" -- A list of DNA reads that are the same length. All letters need to be upper case. (no default)
	"numSeeds" -- Is an integer indicating how many times to seed the genetic algorithm.  Each seed represents one monte carlo run. (no default)
	"k" --  An integer indicating the motif length being searched for. (no default)
	"N" -- The number of iterations before returning the best motif. (no default)

	Return value:	
	"bestMotif" -- For each "Dna" reads return the best length "k" motif match in a list BestMotifs.  
					These list entries are the closest scoring matches after running the gibbsSampler 
					algorithm "numSeeds" times.  Each run gets its own unique starting seed, and that seed
					goes through N cycles before returning the local best match.  The final result BestMotifs 
					represent variants of a single DNA motif.  (This assumes a common 
					DNA motif be present in each DNA read.)
	"""
	results = gibbsSampler(dna, k, N)
	bestScore = score(results)
	bestMotifs = list(results)
	for i in range(1, numSeeds):
		results = gibbsSampler(dna, k, N)
		if score(results) < bestScore:
			bestScore = score(results)
			bestMotifs = list(results)
	return bestMotifs
	
def gibbsSampler(dna, k, N):
	"""A motif finding algorithm that finds one common motif in a list of "dna" reads, returns a list bestMotifs that contains the closes motif match from each string in "dna".
	Keyword arguments:
	"dna" -- A list of DNA reads that are the same length.  
				All letters need to be upper case.  (no default)
	"k" -- An integer indicating the motif length being searched for.  (no default)
	"N" -- The number of iterations before returning the best motif.  (no default)

	Return Value:
	"bestMotifs" -- For each "Dna" reads return the best length "k" motif match in a list BestMotifs.  
					These list entries are the closest scoring matches after running the gibbsSampler 
					algorithm "N" times.  They represent variants of a single DNA motif.  (This assumes a common 
					DNA motif be present in each DNA read.)
	"""
	t = len(dna)
	# randomly select k-mers Motifs = (Motif_1, , Motif_t) in each string from Dna
	motifs = []
	for strand in dna:
		i = random.randrange(len(strand)-k+1)
		substr = strand[i:i+k]
		motifs.append(substr)
		
	bestMotifs = list(motifs)
	bestMotifsScore = score(bestMotifs)
	for j in range(1,N):
		i = random.randrange(t)
		subsetMotifs = motifs[0:i]+motifs[i+1:t]
		replacementMotif = singleReplacementMotif(subsetMotifs, dna[i])
		motifs[i] = replacementMotif
		
		if score(motifs) < bestMotifsScore:
			bestMotifs = list(motifs)
			bestMotifsScore = score(bestMotifs)
	return bestMotifs

#
# helper functions
# functions that come after this line are helper functions
#
	
def singleReplacementMotif(motifs, dna_i):
	"""Select a replacement motif for the "motifs" variable in gibbsSampler, then return the string "replacementKmer".


	Detailed function summary:
	Select a replacement motif for the "motifs" in gibbsSampler, then return the string "replacementKmer".  This replacement is a substring
	from the i_th DNA read ("dna_i").  It will later replace the i_th "motifs"'s
	entry.  Once this motif is replaced, the motifs' list converges a little closer
	to the best match.  This best match is the common motif shared between the 
	DNA reads.

	Keyword arguments:	
	"motifs" -- Is the variable from the motifs' function.
					!!! Except that the i_th row has been removed. !!!
					"dna_i" is the i_th DNA read.  (no default)
	Return value:	Return a replacement k-mer "replacementKmer" to replace the i_th "motif" entry.
	"""

	k = len(motifs[0])
	profile = BuildProfile(motifs)
	
	# Calculate probilities for each k-mer in Dna_i
	kmerDensities = [0 for x in range(len(dna_i)-k+1)]
	for i in range(len(dna_i)-k+1):
		prob = 1
		for j in range(k):
			if dna_i[i+j] == 'A':
				prob *= profile[j][0]
			elif dna_i[i+j] == 'C':
				prob *= profile[j][1]
			elif dna_i[i+j] == 'G':
				prob *= profile[j][2]
			elif dna_i[i+j] == 'T':
				prob *= profile[j][3]
		kmerDensities[i] = prob

	# normalize probabilities
	normalizationTot = sum(kmerDensities)
	for i in range(len(dna_i)-k+1):
		kmerDensities[i] = kmerDensities[i]/normalizationTot

	# construct prefix sum for lookup		
	kmerDensities = list(accumulate(kmerDensities))
	
	# randomly select a k-mer
	randVal = random.random()
	for i in range(len(dna_i)-k+1):
		if randVal < kmerDensities[i]:
			break
	replacementKmer = dna_i[i:i+k]

	return replacementKmer

def BuildProfile(motif):
	""" Build the "profile variable for "gibbsSampler", return "profile". Please see special note.

	Special note:
	profile is a nested list
	profile[k][A,C,G,C]			

	Keyword arguments:
	"motif" -- A list of strings.  "motif", each entry on this list corresponds to a 
						DNA read.  It is a substring of that read, k letters long, which represents a 
						common motif between different DNA sequences. (no default)
	Return value:	
	"profile" -- A 4xk probability matrix.  Each of the four rows corresponds to
					a probability that the k-mer is 'A','C','G', or 'T' at the nth position.  The product 
					of one entry from each column is the probability that a k-mer is a given sequence.
					This matrix is a nested list.  The columns are the list, and the four rows are sub-list.
	"""
	k = len(motif[0])
	profile = [[0 for y in range(4)] for x in range(k)]
	for count in range(k):
		A=0
		C=0
		G=0
		T=0
		
		# Add in Laplace counts to avoid
		# prob densities that are zero or one
		# accelerates alg runtime
		A += 1
		C += 1
		G += 1
		T += 1
		for string in motif:
			if string[count]=='A':
				A+=1
			elif string[count]=='C':
				C+=1
			elif string[count]=='G':
				G+=1
			elif string[count]=='T':
				T+=1
		# Insert frequencies if base A
		profile[count][0] = float(A)/(A+C+G+T)
		# Insert frequencies if base C
		profile[count][1] = float(C)/(A+C+G+T)
		# Insert frequencies if base G
		profile[count][2] = float(G)/(A+C+G+T)
		# Insert frequencies if base T
		profile[count][3] = float(T)/(A+C+G+T)
	return profile

def BuildMotifs(profile, dna, k):
	"""  Build the "motifs" variable for the function gibbsSampler, return "motifs".
	Keyword arguments:
	"profile" -- A 4xk matrix containing the probabilities that 
					the nth base in length k motif will be 'A','C','G', or 'T'.
					The prob that the motif is a particular seq is the product of
					one entry from each of the k rows.  (no default)
	"dna" -- Is a list, length t, of Dna reads being searched for a 
					commone motif.  (no default)
	"k" -- An integer indicating the length of the motif being searched for.  (no default)
	Return value:	
	"motif" -- A list of length k Dna sub-strings.  Each list entry 
					corresponds to a substring in Dna.  The first "motif" entry corresponds 
					to the first "dna" entry, the second to the second, and so on.  Each "motif" 
					entry is chosen to be the most probable substring based on the 
					probabilities given in "profile", representing the most probable motif 
					shared between all the "dna" reads.
	"""
	motif = []
	for string in dna:
		bestSubStr = ''
		for i in range(len(string)+1-k):
			substr = string[i:i+k]
			prob = 1
			bestProb = -1
			for j in range(k):
				if substr[j] == 'A':
					prob *= profile[j][0]
				elif substr[j] == 'C':
					prob *= profile[j][1]
				elif substr[j] == 'G':
					prob *= profile[j][2]
				elif substr[j] == 'T':
					prob *= profile[j][3]
			if prob > bestProb:
				bestProb = prob
				bestSubStr = substr
		motif.append(bestSubStr)
	return motif

def score(motifs):
	""" Counts the number of mismatches between strings in the list "motifs", returns the mismatch count as the "score".
	Keyword arguments:
	"motifs" -- A variable representing a collection of DNA sub-reads, length k.  (no defaults)
	Return Value:	
	"score" -- Returns the number of single base mismatches in "motifs".  The higher the value the worse the score.
	"""
	k = len(motifs[0])
	pattern = []
	for i in range(k):
		A=0
		C=0
		G=0
		T=0
		for string in motifs:
			if string[i]=='A':
				A+=1
			elif string[i]=='C':
				C+=1
			elif string[i]=='G':
				G+=1
			elif string[i]=='T':
				T+=1
		
		#"""
		# For tie counts this chooses A over C
		# C over G, and G over T.
		# However, this does not introduce bias.
		# It has the same results as randomly breaking a tie, 
		# but is much simpler to implement. 		
		if A >= C and A >= G and A >= T:
			pattern.append('A')
		elif C >= G and C >= T:
			pattern.append('C')
		elif G >= T:
			pattern.append('G')
		else:
			pattern.append('T')
		#"""
	pattern = "".join(pattern)
			
	score = 0
	for string in motifs:
		score += hammingDistance(string, pattern)
	return score


def d(kmer, dna):
	""" The functiond d calculates the number of mismatches between the string "kmer" and the string list "dna", this number is returned as "totDist".
	Keyword arguments:	
	"dna" -- A list of strings.  Each entry in "dna" must be longer than "kmer".  (no default)
	"kmer" -- A string representing a DNA kmer.

	Return value:
	"totDist" -- An integer which totals number of mismatches between "kmer" and each list entry in "dna".
	"""
	k = len(kmer)
	motif = []
	totDist = 0
	for phrase in dna:
		localDist = len(phrase)+len(kmer)
		word = ""
		for i in range(len(phrase)-k+1):
			subPattern = phrase[i:i+k]
			if localDist > hammingDistance(kmer, subPattern):
				localDist = hammingDistance(kmer, subPattern)
				word = subPattern
		motif.append(word)
		totDist += localDist
	return totDist


def hammingDistance(str1, str2):
	"""
	Keyword arguments:
	"str1" -- A string.  (no default)
	"str2" -- A string.  (no default)

	Return value:	
	"diffs" -- The number of mismatches between "str1" and "str2".
						Strings can be different lengths, but the mismatch count is
						only includes counts on the shortest length.
	"""
        diffs = 0
        for ch1, ch2 in zip(str1, str2):
                if ch1 != ch2:
                        diffs += 1
        return diffs
		
		
def accumulate(iterable, func=operator.add):
	"""Count the # of differences between equal length strings str1 and str2.

	Keyword arguments:	
	"iterable" -- An iterable "iterable" that can be summed through the 
	function operator.add .  (no default)

	Return value:
	yield the total "total".
	"""
	'Return running totals'
	# accumulate([1,2,3,4,5]) --> 1 3 6 10 15
	# accumulate([1,2,3,4,5], operator.mul) --> 1 2 6 24 120
	it = iter(iterable)
	try:
		total = next(it)
	except StopIteration:
		return
	yield total
	for element in it:
		total = func(total, element)
		yield total
