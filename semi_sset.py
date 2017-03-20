import numpy as np

#-----------------------------------------------------------------------------------------------------------------------------------
class semi_sset():
	def __init__(self, simplices = [], structuretree = []):
		self.simplices = simplices
		self.structuretree = structuretree
		self.topdim = len(self.simplices)-1

	def simplexdim(self,simplex):
		return max(i for i in range(len(self.simplices)) if simplex in self.simplices[i])
		#return depth(self.simplextree(simplex))	

	def Sigma(self,k):
		return list(self.simplices[k])

	def codim1sub(self,i):
		X = semi_sset(self.simplices, self.structuretree[i])
		return X

	def simplexlocate(self,simplex):
		out= ""
		struct = self.structuretree
		if simplex in struct: return str(struct.index(simplex))+"-"
		else:
			for l in struct:
				if isinstance(l,list):
					if isnode(l,simplex):
						out += str(struct.index(l)) + str(self.codim1sub(struct.index(l)).simplexlocate(simplex))
		return out

	def addsimplex(self,name,simplexdim):
			while (simplexdim > len(self.simplices)-1):
				self.simplices.append(set())
		
			self.simplices[simplexdim].add(str(name))

	def simplextree(self,simplex):
		L = self.structuretree
		addresses = self.simplexlocate(simplex)
		firstaddress = addresses.split("0-")[0]
		for direction in firstaddress:
			L = L[int(direction)]
		return L

	
	def d(self,simplex,i):
		if i in range(self.simplexdim(simplex)+1):
			return self.simplextree(simplex)[i+1][0]
	
	def delta(self,simplex):
		n = len(self.Sigma(self.simplexdim(simplex)))
		m = len(self.Sigma(self.simplexdim(simplex)-1))
		out = np.zeros(m)
		
		for i in range(self.simplexdim(simplex)+1):
			out[self.Sigma(self.simplexdim(simplex)-1).index(self.d(simplex,i))] += (-1)**i
		return out

			
	def simplexfaces(self,simplex):
		return set(flatten(self.simplextree(simplex)))

	def simplexboundaries(self,simplex):
		L = []
		for k in range(1,len(self.simplextree(simplex))+1):
			L.append(self.simplextree(simplex)[k][0])
		return L

	def boundary_matrix(self,k):
		self.clean()
		if k in range(self.topdim +1):
			if k == 0: return np.zeros(len(self.Sigma(0)))
			else:
				n = len(self.Sigma(k))
				m = len(self.Sigma(k-1))
				A = np.zeros((m,n))
				for j in range(n):
					A[:,j] = self.delta(self.Sigma(k)[j])
				return A
		elif k == self.topdim+1: return np.zeros(1)	
		else: pass

	def clean(self):
		L = set(flatten(self.structuretree))
		for i in range(len(self.simplices)):
			for x in self.simplices[i]:
				if x not in L:
					self.simplices[i] = self.simplices[i] - {str(x)}
		
		return self	
    
	def betti(self):
		self.clean()

		return [len(self.Sigma(k)) - np.linalg.matrix_rank(self.boundary_matrix(k)) - np.linalg.matrix_rank(self.boundary_matrix(k+1)) for k in range(self.topdim +1)]

	
#--------------------------------------------------------------
def flatten(A):
	if A == []: return A
	if type(A[0]) == list:
		return flatten(A[0]) + flatten(A[1:])
	else: return [A[0]] + flatten(A[1:])

def depth(A):
	if A == []: return -1
	if type(A) == list:
		return 1 + max(depth(a) for a in A)
	else: return -1
	
def isnode(A,x):
	if x in flatten(A): return True
	else: return False

def locate(L,x): 
	out =""
	if x in L: return str(L.index(x))
	else:
		for l in L:
			if isinstance(l,list):
				if isnode(l,x):
					out += str(L.index(l)) + str(locate(L,x))
	return out				

def addressnode(L,address):
	node = L
	for turn in address:
		node = node[int(turn)]
	return node	

def depthcopy(L):
	d = depth(L)
	out = []
	if L == []: return []
	for x in L:
		if type(x) == list: out.append(depthcopy(x))
		else: out.append(str(d) + "-" + str(x))
	return out

def depthcopyreplace(L,A,B):
	d = depth(L)
	out = []
	if L == []: return []
	for x in L: 
		if type(x) == list: out.append(depthcopyreplace(x,A,B))
		else:
			if x in A: out.append(str(B[d]))
			else: out.append(x)
	return out
			

#-----------------------------------------------------------------------------------------------------------

    
X = semi_sset([{"x","y","z"},{"xy","yz","xz"}], [["xy",["x",[]],["y",[]]],["yz",["y",[]],["z",[]]],["xz",["x",[]],["z",[]]]])
print(X.betti())

