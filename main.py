
import gc
import re
import random
import stats
import time
import os
import sys
sets = set()
#from sets import Set
from pathwayStructures import *
import rpy2.robjects as robjects
from deapFunctions import *
import threading
import Queue
import logging
import copy
import string
import argparse

#gc.enable()
def DEAP(pwd, subdir, datadir, edgedir, isPaired, numRotations=100):
###Get all expression files
	expressionFiles=[datadir.split('/')[-1]]
	datadir=re.sub(expressionFiles[0],'',datadir)

###Check for existence of all output directories. If they don't exist, make them.
	outputdir=os.path.join(pwd,'output/')
	if not os.path.exists(outputdir):
		os.makedirs(outputdir)
	logdir=os.path.join(outputdir,subdir,'logs/')
	if not os.path.exists(logdir):
		os.makedirs(logdir)
	outdir=os.path.join(outputdir,subdir,'results/')
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	nulldir=os.path.join(outputdir,subdir,'null/')
	if not os.path.exists(nulldir):
		os.makedirs(nulldir)

###Set up logging
	beginTime=re.sub(r'\.','',str(time.time()))
	logging.basicConfig(filename=os.path.join(logdir,subdir+'.log'),level=logging.INFO, format='[%(levelname)s] (%(threadName)-10s) %(message)s',)

###Iterate over all expression files
	for expressionData in expressionFiles:	
		###Open output file, write header line, define function for future writes
		filename=expressionData.split('.')[0]
		outfile=open(os.path.join(outdir,filename)+'.txt','w')
		outfile.write('PathwayName\t')
		outfile.write('||AbsValue\tpVal\tnullMean\tnullStDev\t||PathwaySubset\n')
		def writeResult(result):
			xmlfile=result[0]
			statvals=result[1]
			outfile.write(xmlfile+'\t')
			outfile.write(str(statvals[2].curval)+'\t'+str(statvals[2].qval)+'\t'+str(statvals[2].mean)+'\t'+str(statvals[2].stdev)+'\t')
			outfile.write(str(statvals[2].pathSubset)+'\n')
			outfile.flush()

		###Load expression data
		expdict = createExpressionDict(os.path.join(pwd,datadir,expressionData))

		###Perform data manipulations. In this case, data rotation
		timeb=time.time()
		logging.info('rotating data')
		rotatedData=buildRotatedData(expdict, numRotations, isPaired)
		logging.info('data rotated in '+str(time.time()-timeb)+' secs.')
		
		### Perform transforms on data based on user input

		expdict=meanExpDict(expdict,isPaired)
		rotatedData=meanDictList(rotatedData,isPaired)

		###Define threadable class for performing all operations on each xmlfile
		class xmlFileBuilder(threading.Thread):
			###Initialize with no result
			def __init__(self,xmlfile):
				self.xmlfile=xmlfile
				self.result=None
				threading.Thread.__init__(self,name=xmlfile)

			###Called after thread completion to get the result
			def get_result(self):
				return [self.xmlfile,self.result]

			###The thread processing
			def run(self):
				timea=time.time()
				logging.info('Building for file: '+self.xmlfile)
				protedges=loadEdgesFromFile(os.path.join(pwd,edgedir,self.xmlfile))
				logging.info('Graph built for '+self.xmlfile+' in '+str(time.time()-timea)+' secs')
				timea=time.time()
				scores=calculateScores(expdict,protedges)
				logging.info('Scores for '+self.xmlfile+' calculated in '+str(time.time()-timea)+' secs:')
				for score in scores:
					logging.info('\tmin= '+str(score[0])+'\tmax= '+str(score[1]))
				nullfile=os.path.join(nulldir,filename+'_'+self.xmlfile)
				timea=time.time()
				nullvalues=buildNull(rotatedData,protedges,nullfile)
				logging.info('null distribution built in '+str(time.time()-timea)+' secs')
				self.result=calculateStats(scores,nullvalues)
				logging.info('Stats calculated')
				logging.info('Build completed')

		###Begin a thread for each file and add it to the queue
		def producer(q,files):
			for xmlfile in files:
				thread=xmlFileBuilder(xmlfile)
				thread.start()
				q.put(thread,True)

		###As files finish in the queue, pull them into the list
		finished=[]
		def consumer(q,numFiles):
			while numFiles > len(finished):
				thread=q.get(True)
				thread.join()
				finished.append(thread.get_result())
		
		###Grab list of all graph files, assuming all .edg files in directory are supposed to be parsed
		xmlfiles=[]
		for i in os.listdir(os.path.join(pwd,edgedir)):
			if re.search('\.edg',i):
				xmlfiles.append(i)
		logging.info(xmlfiles)	

		###Begin multithreading of xml file parsing
		timec=time.time()
		try:
			q=Queue.Queue()
			prod_thread=threading.Thread(target=producer, args=(q,xmlfiles))
			cons_thread=threading.Thread(target=consumer, args=(q,len(xmlfiles)))
			prod_thread.start()
			cons_thread.start()
			prod_thread.join()
			cons_thread.join()
		except Exception as inst:
			logging.error(type(inst))
			logging.error(inst.args)
			logging.error(inst)
			sys.exit()
		logging.info('all graphs built in '+str(time.time()-timec)+' secs.')

		###After all threads are completed, sort and write out the results
		for piece in sorted(finished,key=lambda result:result[1][2].qval):
			writeResult(piece)
		outfile.close()
	#	del rotatedData
	#	del expdict
	#	del finished
		gc.collect()

#Set up the command line interactions
if __name__ == "__main__":
	parser=argparse.ArgumentParser(prog="DEAP",description="Perform Differential Expression Analysis for Pathways")
	parser.add_argument("-w","--working_directory", action="store", metavar="workingDirectory", dest="pwd", required=True,
		help="Directory where output files will be written to and relative file paths for input will be calculated from.")
	parser.add_argument("-a","--analysis_name", action="store", metavar="analysisName", dest="subdir", required=True,
		help="Unique name for this analysis. Will be name of directory where results are stored.")
	parser.add_argument("-t","--tsv_file", action="store", metavar="expressionFile", dest="datadir", required=True,
		help="Location of the TSV file storing expression information. Filepath should either be absolute or relative to working directory.")
	parser.add_argument("-e","--edge_file_directory", action="store", metavar="edgeDir", dest="edgedir", required=True,
		help="Directory containing all the *.edg files for biological pathways. Filepath should either be absolute or relative to working directory.")
	parser.add_argument("-p","--is_paired", choices=["Y","y","N","n"], action="store", dest="isPaired", required=True,
		help="Is data paired? [y/N]. Paired data in the TSV file should be in the form of ratio data. Unpaired should have one condition in the first half of columns and the other condition in the second half of the columns.")
	parser.add_argument("-n","--number_rotations", action="store", dest="numRotations", type=int, default=100, required=False,
		help="The number of random rotations to be performed. Default=100. ")
	args=parser.parse_args()
	print args
	if args.isPaired=='Y' or args.isPaired=='y':
		args.isPaired=True
	else:
		args.isPaired=False
	
	DEAP(args.pwd, args.subdir, args.datadir, args.edgedir, args.isPaired, args.numRotations)
