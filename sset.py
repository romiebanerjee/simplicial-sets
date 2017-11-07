
import numpy as np

#-----------------------------------------------------------------------------------------------------------------------------------
# CLASS SIMPLICIAL SET

class sset():
	def __init__(self, simplices = [], structuretree = []):
		self.simplices = simplices
		self.structuretree = structuretree
		self.topdim = len(self.simplices)-1

	def simplexdim(self,simplex):
		return max(i for i in range(len(self.simplices)) if simplex in self.simplices[i])
		#return depth(self.simplextree(simplex))	

	def Sigma(self,k):
		return list(self.simplices[k])

	def ndg_Sigma(self,k):
		return [x for x in self.simplices[k] if x[0] not in {"0","1","2"}]

	def ndg_simplices(self,k):
		return [x for x in list(self.simplices[k]) if self.is_degenerate(x) == False]

	def codim1sub(self,i):
		X = sset(self.simplices, self.structuretree[i])
		return X

	def dgtree(self,simplex):
		tree = []
		dim = self.simplexdim(simplex)
		if dim >= self.topdim: return [simplex,[]]
		else: 
			tree = [simplex]
			if simplex[0].isdigit():
				for k in range(int(simplex[0])+1):
					self.simplices[dim+1].add(str(k)+"-"+str(simplex))
					tree.append(self.dgtree(str(k)+"-"+str(simplex)))
			else:
				for k in range(dim+1):
					self.simplices[dim+1].add(str(k)+"-"+str(simplex))
					tree.append(self.dgtree(str(k)+"-"+str(simplex)))
		return tree			
	
	def is_degenerate(self,simplex):
		ans = 1
		for i in range(len(self.simplices)):
			for x in self.simplices[i]:
				if simplex in list(self.simplexdegeneracies(x))[1:]: ans = ans*0
		if ans == 0: return True
		else: return False
		
	def simplexdegeneracies(self,simplex):
		return set(flatten(self.dgtree(simplex)))

	
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


	def add_alldegeneracies(self):
		for i in range(self.topdim + 1):
			for simplex in self.simplices[i]:
			 	self.dgtree(simplex)

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
		#n = len(self.Sigma(self.simplexdim(simplex)))
		m = len(self.Sigma(self.simplexdim(simplex)-1))
		out = np.zeros(m)
		
		for i in range(self.simplexdim(simplex)+1):
			if self.d(simplex,i)[0] not in {"0","1","2"}:
			#if self.is_degenerate(self.d(simplex,i)) == False:
				out[self.Sigma(self.simplexdim(simplex)-1).index(self.d(simplex,i))] += (-1)**i
		return out

			
	def simplexfaces(self,simplex):
		return set(flatten(self.simplextree(simplex)))

	def simplexboundaries(self,simplex):
		L = []
		for k in range(1,len(self.simplextree(simplex))+1):
			L.append(self.simplextree(simplex)[k][0])
		return L

	
	def simplexcollapse(self,simplex):
		print("\n" + "collapsing simplex " + str(simplex) + "...." + "\n")	
		Q = sset(self.simplices, self.structuretree)
		F = Q.simplexfaces(simplex)
		D = set()
		for x in F:
			D = D.union(Q.simplexdegeneracies(x))
		
		collapse = []
		for i in range(len(self.simplices)):
			collapse.append(set())
		for i in range(len(collapse)):
			for x in self.simplices[i]:
				if x not in F.union(D): collapse[i].add(x)
				else: pass
				
		Q1 = sset(collapse,[])
		new_vertex = str(simplex)+"_0"
		Q1.addsimplex(new_vertex,0)
		new_dgs = sorted(list(Q1.simplexdegeneracies(new_vertex)), key = len)	
		

		collapsetree = depthcopyreplace(self.structuretree, F.union(D) , new_dgs)
		return sset(collapse,collapsetree).clean()
		
	
	def clean(self):
		L = set(flatten(self.structuretree))
		for i in range(len(self.simplices)):
			for x in self.simplices[i]:
				if x not in L:
					self.simplices[i] = self.simplices[i] - {str(x)}
		
		return self	
			
																				
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

	def betti(self):
		self.clean()

		return [len(self.ndg_Sigma(k)) - np.linalg.matrix_rank(self.boundary_matrix(k)) - np.linalg.matrix_rank(self.boundary_matrix(k+1)) for k in range(self.topdim +1)]


	def collapse(self):
		for i in range(len(self.simplices)):
			for x in self.simplices[i]:
				sub_sset = sset([],self.simplextree(x))
				if sum(sub_sset.betti()) == 1: #sub_sset has trivial homology over a field
					self = self.simplexcollapse(x)
		return self			
				

#---------------------------------------------------------------------------
#CLASS SEMI SIMPLICIAL SET (or DELTA COMPLEX)

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

	

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

#CLASS SIMPLCIAL COMPLEX 

class scomplex():
	def __init__(self, simplices = []):
		self.simplices = list(set(flatten([self.simplexfaces(simplex) for simplex in simplices])))
		self.vertices = [simplex for simplex in self.simplices if len(simplex)-1 == 0]
		self.dim = max(len(x)-1 for x in self.simplices)
		self.topdimsimplices = [simplex for simplex in self.simplices if len(simplex)-1 == self.dim]
		self.topsimplices = self.topdimsimplices  
		self.facetree = self.facetree()

	def Sigma(self,k):
		return [simplex for simplex in self.simplices if len(simplex)-1 == k]
	
	def d(self,simplex,i):
		if i in range(len(simplex)):
			return simplex[0:i] + simplex[i+1:len(simplex)]
		else:pass

	def simplextree(self, simplex):
		simplextree = []
		if len(simplex) == 1: return [simplex,[]]
		else:
			simplextree = [simplex]
			for i in range(len(simplex)):
				simplextree.append(self.simplextree(self.d(simplex,i)))
		return simplextree
	

	def simplexfaces(self,simplex):
		return flatten(self.simplextree(simplex))


	def facetree(self):
		return [self.simplextree(simplex) for simplex in self.topsimplices]	

	def D(self,k):
		if k in range(self.dim +1):
			if k == 0: return np.zeros(len(self.Sigma(0)))
			else:
				n = len(self.Sigma(k))
				m = len(self.Sigma(k-1))
				A = np.zeros((m,n))
				for i in range(m):
					for j in range(n):
						for l in range(k+1):
							if self.Sigma(k-1)[i] == self.d(self.Sigma(k)[j],l):
								A[i,j] = (-1)**l
							else: pass
				return A
		elif k == self.dim+1: return np.zeros((1,1))	
		else: pass


	def betti(self):
		return [len(self.Sigma(k)) - np.linalg.matrix_rank(self.D(k)) - np.linalg.matrix_rank(self.D(k+1)) for k in range(self.dim +1)]


	
#----------------------------------------------------------------------------------------------------------------------------------------------------------------	

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------



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
			

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------

