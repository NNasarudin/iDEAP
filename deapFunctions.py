
import logging
import re
import random
import stats
import time
import os
import sys
sets = set()
from pathwayStructures import *
import rpy2.robjects as robjects
import threading
import Queue

classKey='class'
def loadEdgesFromFile(filename):
	"""Load edges from a formatted '.edg' file into a list of Edge objects
	Proper '.edg' formatting is a series of lines, each representing a single edge, with the format:
	[Tab delimited list of reactants]::[Tab delimited list of products]::[True if edge is catalytic, False otherwise
	Args:
		filename: absolute path to the .edg file
	Returns:
		list of Edge objects
	"""
	edges=[]
	op=open(filename,'r')
	for line in op:
		line=line.strip()
		if line !='':
			if(len(line.split('::'))==3):
				r, p, t=line.split('::')
				boo=False
				if t=='True' or t=='true' or t=='T' or t=='t':
					boo=True
				reactants=Set()
				for pos in r.split('\t'):
					if pos !='':
						reactants.add(pos)
				products=Set()
				for pos in p.split('\t'):
					if pos!='':
						products.add(pos)
				edges.append(Edge(reactants,products,boo))
			else:
				logging.warning('Warning: error in line in file '+filename+'. Line has '+str(len(line.split('::'))-1)+' sets of \'::\'. Expected 2.')
	return edges

def scoreCalculator(expression_dict=dict(), leafs=[],reactants=dict(), index=0):
	"""scoreCalculator is a thread-safe method for calculating the DEAP score.
	Args:
		expression_dict: should be a dictionary mapping ids used in edges (usually uniprot) to a list of expressionvalues
		leafs: is a list of all proteins which are in an edge's reactant, but no edge's products
		reactants: is a dictionary mapping protein strings to a list of edges containing that protein as a reactant
		index: is the index in the list of expression values corresponding to the set of expression values we are currently looking at
	Returns:
		a list containing the minimum and maximum score for the pathway
	"""
	#Namespace class exists to handle scope issues arising during competing threads
	class Namespace: pass
	ns=Namespace()
	ns.maxscore=0
	ns.minscore=0
	ns.visited=Set()
	ns.edgeScore=dict()
	ns.maxSub=''
	ns.minSub=''
	def genericScore(proteins):
		"""given an iterable of protein name strings, returns the combined score for those proteins 
		intended for use when calculating the score of a single node
		Args:
			proteins: a list of strings of protein names
		Returns:
			a float that is the DEAP score
		"""
		scores=Set()
		for p in proteins:
			for protein in p.split(','):
				protein=protein.strip()
				if protein in expression_dict and len(expression_dict[protein])>index:
					scores.add(expression_dict[protein][index])
		sum=0
		count=0.0
		for score in scores:
			sum=sum+score
			count=count+1.0
		if count==0:
			return 0
		return sum
	def reactantScore(edge):
		return genericScore(edge.reactants)
	def productScore(edge):
		return genericScore(edge.products)
	def edgeTypeMultiplier(edge):
		if edge.typ:
			return 1
		return -1
	def recursiveScore(edge,recdepth=0,alreadyvisited=Set()):
		"""Recursive method for calculating the DEAP score starting at informative edge.
		Implement DMSP algorithm where searching starting from informative node.
		Args:
			edge: object of class Edge (the most informative) where the score being calculated begins
			recdepth: integer representing the current depth of recursion to prevent infinite loops
			alreadyvisited: a set containing a list of all edges which have already been visited
		Returns:
			a list with two float values and two strings. [minvalue, maxvalue, minString, maxString]
				Where minString and maxString 
		"""
		scoreArray=[0,0,'','']
		if recdepth>=100:
			logging.warning('Entered cycle. Recursion depth too great (>100)')
			return scoreArray
		alreadyvisited.add(edge)
		for reactant in edge.reactants:
			ns.visited.add(str(reactant))
		
		#For every one of this edge's products, see if those products are in any other edge's reactants
		#If there are such edges and we have not already visited this edge, then determine the recursive score for the edge
		for product in edge.products:
			if str(product) in reactants:
				for nextedge in reactants[str(product)]:
					if nextedge not in alreadyvisited:
						alreadyvisited.add(nextedge)
						#see if we have already calculated the recursive score for this edge
						#if we have, return the score from out dictionary
						#if not, calculate the score recursively
						if nextedge in ns.edgeScore:
							tempscore=ns.edgeScore[nextedge]
						else:
							tempscore=recursiveScore(nextedge,recdepth+1,alreadyvisited=alreadyvisited.copy())
						if tempscore[0]<scoreArray[0]:
							scoreArray[2]=tempscore[2]
							scoreArray[0]=tempscore[0]
						if tempscore[1]>scoreArray[1]:
							scoreArray[3]=tempscore[3]
							scoreArray[1]=tempscore[1]
		#Calculate the product score only if there is no value returned from the recursive function
		#If there has been a value returned, then product score calculation would be duplicative
		if scoreArray[0]==0:
			scoreArray[2]=str(edge.products)
			scoreArray[0]=productScore(edge)
		if scoreArray[1]==0:
			scoreArray[3]=str(edge.products)
			scoreArray[1]=productScore(edge)
		reactScore=reactantScore(edge)
		edgeMult=edgeTypeMultiplier(edge)
		#calculate the new maximum and minimum scores
		calcWithPrevMin=reactScore+edgeMult*scoreArray[0]
		calcWithPrevMax=reactScore+edgeMult*scoreArray[1]
		scoreArray[1]=max(calcWithPrevMin,calcWithPrevMax)
		scoreArray[0]=min(calcWithPrevMin,calcWithPrevMax)

		#Format the pathwaySubsetText
		arrow='-->' if edgeMult==1 else '--|'
		if scoreArray[0]==calcWithPrevMin:
			scoreArray[2]=str(edge.reactants)+arrow+scoreArray[2]
			scoreArray[3]=str(edge.reactants)+arrow+scoreArray[3]
		else:
			tempStr=scoreArray[2]
			scoreArray[2]=str(edge.reactants)+arrow+scoreArray[3]
			scoreArray[3]=str(edge.reactants)+arrow+tempStr

		#change the namespace max or min if appropriate
		if scoreArray[1]>ns.maxscore:
			ns.maxSub=scoreArray[3]
			ns.maxscore=scoreArray[1]
		if scoreArray[0]<ns.minscore:
			ns.minSub=scoreArray[2]
			ns.minscore=scoreArray[0]
		#reset the max or min score if appropriate
		if scoreArray[1]<0:
			scoreArray[1]=0
			scoreArray[3]=''
		if scoreArray[0]>0:
			scoreArray[0]=0
			scoreArray[2]=''
		#add the edge score to the edge score dictionary
		ns.edgeScore[edge]=scoreArray
		return scoreArray
	startTime=time.clock()
	#For every leaf node, calculate the recursive score
	for leaf in leafs:
		for reactant in reactants[leaf]:
			recursiveScore(reactant)
	#See if any edges have not been visited
	#if there are such edges, then calculate those recursive scores
	allPos=Set(reactants.keys())
	remaining=allPos.difference(ns.visited)
	while len(remaining) != 0:
		logging.info('remaining: '+str(len(remaining)))
		logging.info(str(remaining))
		for reactant in reactants[remaining.pop()]:
			recursiveScore(reactant)
		remaining=allPos.difference(ns.visited)
	return [ns.minscore,ns.maxscore,ns.minSub,ns.maxSub]

def edgeTransforms(edges):
	"""Helper method transforming a list of edges into input for calculateScores
	Args:
		edges: a list of edges
	Returns:
		a list containing two objects [reactants, leafs]
		reactants: a dictionary mapping protein strings to edges in which that protein is a reactant
		leafs: a list of edges which have at least one reactant which is not in any edge's products
	"""
	reactants=dict()
	products=Set()
	for edge in edges:
		for product in edge.products:
			products.add(product)

		for reactant in edge.reactants:
			if reactant not in reactants:
				reactants[reactant]=[edge]
			else:
				reactants[reactant].append(edge)
	leafs=Set()
	for reactant in reactants.keys():
		if reactant not in products:
			leafs.add(reactant)
	if len(leafs)==0 and len(reactants)>0:
		leafs.add(reactants.keys()[0])
	return [reactants,leafs]


def calculateScores(expression_dict, edges):
	"""Calculates scores for an entire expression dictionary and a list of edges
	Args:
		expression_dict: a dictionary whose key is an id and value is a list of expression values
		edges: a list of Edge objects with ids corresponding to expression_dict ids
	Returns:
		a list of lists, each containing the min and max score for every set of expression values
	"""
	reactants,leafs=edgeTransforms(edges)
	return calculateScoresTransformed(expression_dict=expression_dict,reactants=reactants,leafs=leafs)

def calculateScoresTransformed(expression_dict, reactants,leafs):
	"""Calculates scores for an entire expression dictionary
	Args:
		expression_dict: a dictionary whose key is an id and value is a list of expression values
		reactants: a dictionary mapping protein strings to edges in which that protein is a reactant
		leafs: a list of edges which have at least one reactant which is not in any edge's products
	Returns:
		a list of lists, each containing the min and max score for every set of expression values
	"""
	numVals=len(expression_dict[expression_dict.keys()[0]])
	finished = []
	for index in range(0,numVals): 
		finished.append(scoreCalculator(expression_dict=expression_dict,index=index,reactants=reactants,leafs=leafs))
	return finished

def buildRotatedData(expDict, numRotations, isPaired):
	"""Rotate an entire dictionary of data 1000 times
	Args:
		expDict: a dictionary mapping protein ids to a list of float expression values
	globals().update(args.__dict__)
		numRotations: the number of rotations which should be performed
		isPaired: true if data is paired, false otherwise
	Returns:
		a list of dictionaries mapping protein ids to a list of float expression values0
	"""
	rotationList=[]
	for i in range(0,numRotations):
		if i%100==0:
			logging.info(str(i)+" rotations performed")
		rotationList.append(singleRotateData(expDict, isPaired))
	return rotationList

def singleRotateData(expDict,isPaired):
	"""Using R scripts (rotate.r) in combination with rpy2, perform rotations of data
	Args:
		expDict: a dictionary mapping protein ids to a list of float expression values
		isPaired: true if data is paired, false otherwise
	Returns:
		a dictionary with rotated expression values
	"""
	scriptPath= os.path.realpath(__file__)
	rfil=open(os.path.join(os.path.split(scriptPath)[0],'rotate.r'),'r')
	robjects.r(re.sub('\r','',rfil.read()))
	fitResidual=robjects.globalenv['fit.resid']
	rotate=robjects.globalenv['rotate.adjusted.y']
	expressionMatrixValues=[]
	keys=expDict.keys()
	numSamples=len(expDict[keys[0]])
	for key in keys:
		for index in range(0,numSamples):
			expressionMatrixValues.append(expDict[key][index])
	expressionMatrixVector=robjects.FloatVector(expressionMatrixValues)
	expressionMatrix=robjects.r['matrix'](expressionMatrixVector,nrow=numSamples)
	if not isPaired:
		sampleOneVector=robjects.FloatVector([1]*numSamples)
		fit=fitResidual(y=expressionMatrix,x=sampleOneVector)
	else:
		fit=fitResidual(y=expressionMatrix)
	rotatedMatrix=rotate(fit)
	rotatedDic=dict()
	for keyIndex in range(0,len(keys)):
		rotatedList=[]
		for sampleNum in range(0,numSamples):
			rotatedList.append(rotatedMatrix[keyIndex*numSamples+sampleNum])
		rotatedDic[keys[keyIndex]]=rotatedList
	
	#If we have class labels for unpaired data, we don't want to rotate them
	if classKey in expDict:
		rotatedDic[classKey]=expDict[classKey]

	return rotatedDic

def unpairedMean(expVals,classVec):
	"""For unpaired data, transforms data for calculation in DEAP
	Args:
		expVals: list of unpaired expression values
		classVec: labels for each sample. 0 for case and 1 for control.
	Returns:
		list containing a single expression value
	"""
	if len(expVals)==0:
		return expVals
	a=[]
	b=[]
	for i in range(0,len(expVals)):
		if classVec[i]==0:
			a.append(expVals[i])
		else:
			b.append(expVals[i])

	meanA=stats.mean(a)
	meanB=stats.mean(b)

	return [meanA-meanB]

def unpairedTrans(expDict):
	"""Performs transformations on entire dictionary of unparied expression data
	Args:
		expDict: a dictionary containing protein as the key and a list of expression values as the value
	Returns:
		dictionary with protein as the key and a list with a single expression value
	"""
	edt=expDict.copy()
	
	if classKey not in edt:
		raise Exception('Error: no class labels, but data is unpaired. Please include a row with name "'+classKey+'" and values of either 0 for case or 1 for control.')
		
	classVec=edt[classKey]
	
	#Check the class labels are all 0 or 1
	zeroOne=[a==0 or a==1 for a in classVec]
	if len([x for x in zeroOne if x]) != len(zeroOne):
		raise Exception('Error: inappropriate class label. Please ensure that all values in the "'+classKey+'" row are either 0 for case or 1 for control.')

	for key in edt:
		if key != classKey:
			edt[key]=unpairedMean(edt[key],classVec)
	return edt
def pairedTrans(expDict):
	"""Performs transformations on entire dictionary of paried expression data
	Args:
		expDict: a dictionary containing protein as the key and a list of expression values as the value
	Returns:
		dictionary with protein as the key and a list with a single expression value
	"""
	for key in expDict.keys():
		expDict[key]=[stats.mean(expDict[key])]
	return expDict

def meanExpDict(expDict, isPaired):
	"""Calculates the appropriate transformations for an entire dictionary of expression values. Sends off data as either paired on unpaired.
	Args:
		expDict: a dictionary containing protein as the key and a list of expression values as the value
		isPaired: true if data is paired, false otherwise
	Returns:
		dictionary with protein as the key and a list with a single expression value
	"""
	if isPaired:
		return pairedTrans(expDict)
	else:
		return unpairedTrans(expDict)

def meanDictList(permutationList, isPaired):
	"""Performs transformations on a list of dictionaries, which each map proteins to expression values. 
	Args:
		permutationList: a list of dictionaries whose keys are proteins and values are a list of expression values
		isPaired: true if data is paired, false otherwise
	Returns:
		a list of dictionaries whose keys are proteins and values are a list with a single expression value.
	"""
	newList=[]
	for expDict in permutationList:
		newList.append(meanExpDict(expDict,isPaired))
	return newList

def buildNull(permutationList, edges,xmlfile):
	"""Build a null distribution of values for an individual graph
	Args:
		permutationList: a list of dictionaries mapping protein ids to a list of float expression values which will be used to calculate null values for the graph
		edges: a list of Edge objects representing representing the graph
		xmlfile: the name of the xmlfile from which the graph was generated
	Returns:
		a list of lists of pairs of DEAP scores (min, max)
	"""
	scoreset=[]
	#See if a file containing the null values already exists by trying to open it
	#If the file exists, save time by using those scores
	#Otherwise, calculate and write out the null scores
	try:
		openxmlfile=open(xmlfile+'.null','r')		
		for line in openxmlfile:
			line=line.strip()
			pieces=line.split('::')
			curset=[]
			for piece in pieces:
				tabs=piece.split('\t')
				if len(tabs)==2:
					curset.append([float(tabs[0]),float(tabs[1])])
			scoreset.append(curset)
	except IOError:
		outxmlfile=open(xmlfile+'.null','w')
		reactants,leafs=edgeTransforms(edges)
		for i in range(0,len(permutationList)):
			tabs=calculateScoresTransformed(expression_dict=permutationList[i],reactants=reactants,leafs=leafs)
			scoreset.append(tabs)
			for scorePair in tabs:
				outxmlfile.write(str(scorePair[0])+'\t'+str(scorePair[1])+'::')
			outxmlfile.write('\n')
	return scoreset

def createExpressionDict(mappingfil):	
	"""Create expression dictionaries from input file
	Args:
		mappingFil: the full path to a file containing a tab delimited list of expression values. The first tab should represent the protein id, all other tabs should be expression values
	Returns:
		a dictionary mapping protein ids to a list of float expression values
	"""
	openfil=open(mappingfil,'r')
	expressionDict=dict()
	for line in openfil:
		line=re.sub('\n','',line)
		tabs=line.split('\t')
		id=tabs[0]
		val=[]
		for i in range(1, len(tabs)):
			tabs[i]=tabs[i].strip()
			if tabs[i]!='':
				val.append(float(tabs[i]))
		expressionDict[id]=val
	return expressionDict
	
	
def calculateStats(scores, nullvalues):
	"""Calculate statistics based on actual DEAP scores and null value scores
	Args:
		scores: a list of pairs of floats representing expression scores
		nullvalues: a list of lists of pairs of floats representing null expression scores
	Returns:
		A triplet of StatOutput objects (minimum, maximum, absolute)
	"""
	minvalues=[]
	maxvalues=[]
	absvalues=[]
	for pairs in nullvalues:
		minList=[]
		maxList=[]
		for pair in pairs:
			minList.append(pair[0])
			maxList.append(pair[1])
		minVal=stats.mean(minList)
		maxVal=stats.mean(maxList)
		minvalues.append(minVal)
		maxvalues.append(maxVal)
		absvalues.append(max(maxVal,abs(minVal)))

	maxmean=stats.mean(maxvalues)
	maxstdev=stats.stdev(maxvalues)
	minmean=stats.mean(minvalues)
	minstdev=stats.stdev(minvalues)
	absmean=stats.mean(absvalues)
	absstdev=stats.stdev(absvalues)

	realMaxVals=[]
	realMinVals=[]
	realMaxPathSubset=''
	realMinPathSubset=''
	for score in scores:
		realMaxVals.append(score[1])
		realMinVals.append(score[0])
		realMaxPathSubset=score[3]
		realMinPathSubset=score[2]
	maxVal=stats.mean(realMaxVals)
	minVal=stats.mean(realMinVals)
	maxCount=0
	for val in maxvalues:
		if val>=maxVal:
			maxCount+=1
	maxPval=float(maxCount)/float(len(maxvalues))
	maxStats=StatOutput(mean=maxmean,stdev=maxstdev,curval=maxVal,qval=maxPval)
	minCount=0
	for val in minvalues:
		if val<=minVal:
			minCount+=1
	minPval=float(minCount)/float(len(minvalues))
	minStats=StatOutput(mean=minmean,stdev=minstdev,curval=minVal,qval=minPval)
	absVal=max(maxVal,abs(minVal))
	absPathSubset=realMinPathSubset
	if maxVal==absVal:
		absPathSubset=realMaxPathSubset
	absCount=0
	for i in range(0,len(maxvalues)):
		if max(maxvalues[i],abs(minvalues[i]))>=absVal:
			absCount+=1
	absPval=float(absCount)/float(len(maxvalues))
	absStats=StatOutput(mean=absmean,stdev=absstdev,curval=absVal,qval=absPval,pathSubset=absPathSubset)
	return [minStats,maxStats,absStats]

#def setIsPaired(paired):
	"""Sets the global value for paired data for future calculations
	Args:
		isPaired: true if data is paired, false otherwise
	"""
#	global isPaired
#	isPaired=paired
