# Python program to print connected
# components in an undirected graph
class Graph:
	# init function to declare class variables
	def __init__(self,V):
		self.V = V
		self.adj = [[] for i in range(V)]

	def DFSUtil_rec(self, temp, v, visited):
		# Mark the current vertex as visited
		visited[v] = True
		# Store the vertex to list
		temp.append(v)
		# Repeat for all vertices adjacent
		# to this vertex v
		for i in self.adj[v]:
			if visited[i] == False:
				# Update the list
				temp = self.DFSUtil_rec(temp, i, visited)
		return temp
	
	def DFSUtil(self, temp, v, visited):
		stack = []
		stack.append(v)
		while(len(stack) > 0):
			# pop vertex from stack
			v = stack.pop()
			if not visited[v]:
				# Mark the current vertex as visited
				visited[v] = True
				# Store the vertex to list
				temp.append(v)
			# Repeat for all vertices adjacent
			# to this vertex v
			for i in self.adj[v]:
				if not visited[i]:
					# Update the list
					stack.append(i)
		return temp
	
	# method to add an undirected edge
	def addSub(self, v, w):
		self.adj[v].append(w)
		self.adj[w].append(v)

	# Method to retrieve connected components
	# in an undirected graph
	def connectedComponents(self):
		visited = [False]*self.V
		cc = []
		for v in range(self.V):
			if not visited[v]:
				temp = []
				cc.append(self.DFSUtil(temp, v, visited))
		return cc
