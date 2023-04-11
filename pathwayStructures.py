
class StatOutput:
	"""A class for holding the statistical output created by the DEAP algorithm
	"""
	def __init__(self, mean, stdev, curval, qval,pathSubset=''):
		self.mean=mean
		self.stdev=stdev
		self.curval=curval
		self.qval=qval
		self.pathSubset=pathSubset

class Edge:
	"""A class for representing edges between sets of reactants and products with a relationship between them
	"""
	def __init__(self, react, prod, typed):
		self.reactants=react
		self.products=prod
		self.typ=typed
	def __str__(self):
		return str(self.reactants)+str(self.products)+str(self.typ)
	def __eq__(self,other):
		return str(self)==str(other)
	def __ne__(self,other):
		return str(self)!=str(other)
	def __hash__(self):
		return hash(str(self))

