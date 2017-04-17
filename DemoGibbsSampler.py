"""
DemoGibbsSampler.py
Using Python 2.7
  

Created by Scott Czopek on 3/29/17.
Copyright (c) 2017 __MyCompanyName__. All rights reserved.


!!!!
!!!!
Read:
This script demos the Gibbs Sampler Motif finding algorithm.
The test data is a collection of mouse DNA reads.  Each read 
contains a binding site for the Zinc fingered GATA4 promoter.
The algorithm locates and prints each binding site.  You can
compare the algorithm's answer to the real answer, by opening 
the solution file.  In that file the real motif binding site 
appears in capital letters.

This algorithm is an improvement over the brute force algorithm 
reducing the brute force exponential runtime to a probablistic 
genetic algorithm polynomial run time.
!!!!
!!!!
"""

import sys
import os
import inspect
import re

# determine path to script, then
# import GibbsSampler module
filename = inspect.getframeinfo(inspect.currentframe()).filename
path = os.path.dirname(os.path.abspath(filename))
sys.path.insert(0, path)
import GibbsSampler


def DemoMotifFinder(numSeeds, k, N):
	""" This demos the gibbsSampler algorithm, but finding the common motif for the Zinc Fingered GATA4 promoter.""" 
	filename = inspect.getframeinfo(inspect.currentframe()).filename
	path = os.path.dirname(os.path.abspath(filename))
	# Load mm9 DNA from file.
	# Each read contains a motif, specifically a Zinc Fingered GATA4 promoter site.
	# This software will find that promoter site.
	mm9Loc = path + "/mm9Gata4MotifCollection.txt"
	
	#mm9File = open("/Users/scottczopek/Documents/bioPython/pythonMotifFindingDemo/mm9Gata4MotifCollection.txt",'r')
	mm9File = open(mm9Loc, 'r')
	
	Dna = []
	for line in mm9File:                                                    
		if line[0:3] != '>mm':
			line = line.strip()
			Dna.append(line)

	BestMotif = GibbsSampler.multipleSeedsGibbsSampling(Dna, numSeeds, k, N)
	bestScore = GibbsSampler.score(BestMotif)
	
	mm9File.close()
	return BestMotif


numSeeds = 20
k = 11
N = 1000
BestMotif = DemoMotifFinder(numSeeds, k, N)
bestScore = GibbsSampler.score(BestMotif)
print("Gibbs Sampler Motifs")
print("BestScore:   ", bestScore)
#print(BestMotif)

mm9SolLoc = path + "/mm9Gata4Solutions.txt"
mm9Solutions = open(mm9SolLoc, 'r')
SolutionsMotif = []
for line in mm9Solutions:
	if line[0:3] != '>mm':
		line = re.sub('[^A-Z]', '', line)
		SolutionsMotif.append(line)
print("Real Motifs")
print("Real Score:   ", GibbsSampler.score(SolutionsMotif))
print()
print("Algorithms Best Pick", "     Match/Wrong      ", "Real Motif")
for i in range(len(SolutionsMotif)):
	if SolutionsMotif[i] == BestMotif[i]:
		print(BestMotif[i], "     Match     ", SolutionsMotif[i])
	else:
		print(BestMotif[i], "     Wrong     ", SolutionsMotif[i])

	
