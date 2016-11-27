class simplex():
	def __init__(self, dim = None, name = None, faces = []):
		self.name = name
		self.dim = dim
		self.faces = faces
	
	def is_simplex(self):
		c=1
		for face in self.faces:
			if face.dim == self.dim -1:
				c = c*1
			else: c = c*0
		
	
		if c==1 and len(self.faces) == self.dim + 1: 
			return True
		else: 
			return False

	def display(self):
		if self.is_simplex() == True:
			print(self.name)
			print(self.dim)
			for face in self.faces:
				print(face.name)
		else:
			print("not a simplex")

class SimplicialSet(simplex):
	def __init__(self, name = None, NDsimplices = set(), Dsimplices = set()):
		self.name = name
		self.NDsimplcices = NDsimplices
		self.Dsimplices = Dsimplices
	
	


A = simplex(2,"A")
a = simplex(1,"a")

A.faces = [a,a,a]

A.is_simplex()
A.display()
X = SimplicialSet("X", {A,a},set())



