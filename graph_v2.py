class Graph:
	def __init__(self, name, vertices = set(), edges = {}):
		self.name = name
		self.vertices = vertices
		self.edges = edges
		self.num_vertices = len(vertices)
	
	def is_graph(self):
		x = 1
		if type(self.vertices) == set and type(self.edges) == dict: pass
		else: x = 0

		for e,faces in self.edges.items():
			if type(faces) == list and len(faces) == 2 and set(faces) <= self.vertices: 
				x = x*1
			else: x = x*0
		if x == 1: return True
		else: return False
	
	def display(self):
		if self.is_graph() == True: 
			print("The graph " + str(self.name) + " : ")
			print("vertices : " , self.vertices - {"null"})
			#print("edges" , self.edges)
			print("edges : ")
			for e, faces in self.edges.items():
				if faces != ["null","null"]:
					print(str(e) + " : " + str(faces[0]) +  " --> " + str(faces[1]))
			print("\n")
		else: print(str(self.name) + " not a graph :( ")

	def add_vertex(self,node):
		self.num_vertices += 1
		self.vertices.add(node)
		
	def add_edge(self, label, faces):
		if type(faces) == list and len(faces) ==2 and set(faces) <= self.vertices:
			self.edges[label] = faces
		else:pass	

	def remove_vertex_set(self, vertex_set):
		if vertex_set <= self.vertices:

			self.vertices = self.vertices - vertex_set
			new_vertex = ""
			for v in vertex_set:
				new_vertex += str(v)
			self.add_vertex(new_vertex)

			for e,faces in self.edges.items():
				if faces[0] in vertex_set:
					faces[0] = new_vertex
				else: pass
				if faces[1] in vertex_set:
					faces[1] = new_vertex
				else: pass
		else: print(" not a valid vertex set :( ")

	def collapse_edge(self,e):
		self.add_vertex("null")
		if e in self.edges.keys():
			vertex_set = set(self.edges[e])
			self.remove_vertex_set(vertex_set)
			self.edges[e] = ["null","null"]
			
		else: print(" not a valid edge ")
	
	def collapse(self):
	#	print("collapsing " + str(self.name) + " ... " + "\n")
		for e in self.edges.keys():
			if self.edges[e][0] != self.edges[e][1]:
				self.collapse_edge(e)
	
		return self
	def homology(self):
		self = self.collapse()
		betti_number = 0
		for e in self.edges.keys():
			if self.edges[e] == ["null","null"]:
				betti_number += 0
			else: betti_number += 1
		
		components = 0
		for v in self.vertices:
			if v == "null": components += 0
			else: components += 1
		
		print("Computing " + str(self.name) + " homology ..." + "\n")

		print("rank(H_0(" + str(self.name) + ")) = " , components, "\n")
		print("rank(H_1(" + str(self.name) + ")) = " , betti_number, "\n")
	
		return components, betti_number

#--------------------------------------------------------------------------------------------
import numpy as np

def data_to_graph(A, scale):
	edges = []
	#if type(A) == numpy.ndarray:
	m,n = A.shape
	for i in range(m):
		for j in range(i,m):
			if j!=i and np.sum((A[i]-A[j])*(A[i]-A[j]))<= scale**2:
				edges.append([i,j])
	
	vertices = set(range(m))

	edges_dict = {}
	for e in edges:
		edge_name = ""
		for x in e:
			edge_name += str(x)
		edges_dict[edge_name] = e
	
	G = Graph("Out", vertices, edges_dict)
	return G


def persistent_homology(A, scales):
	PH = []
	print("SCALES = ", scales)
	for scale in scales:
		K = data_to_graph(A,scale)
		#print("SCALE = ", scale)
		print("The graph from SCALE " + str(scale) + " has " + str(len(K.vertices)) + " vertices and " + str(len(K.edges)) + " edges " + "\n")
		#K.collapse()
		#K.homology()
		PH.append(K.homology())
		
	print("PH = ", list(PH))

A = 10*np.random.random((100,10))
			

print(A)

scales = np.linspace(0,20,num=10)

persistent_homology(A,scales)

