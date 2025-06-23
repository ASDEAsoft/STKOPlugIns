from typing import List

# Python program to obtain connected
# components in an undirected graph
class graph:
	
	# init function to declare class variables
	def __init__(self, V : int):
		# number of vertices
		self.V = V
		# adjacency list
		self.adj : List[List[int]] = [[] for i in range(V)]

    # recursive function to perform DFS (Depth First Search)
	# and store the vertices in a list
    # visited[] is a boolean array to keep track of visited vertices
	def _DFSUtil_recursive(self, temp, v, visited):
		# Mark the current vertex as visited
		visited[v] = True
		# Store the vertex to list
		temp.append(v)
		# Repeat for all vertices adjacent
		# to this vertex v
		for i in self.adj[v]:
			if visited[i] == False:
				# Update the list
				temp = self._DFSUtil_recursive(temp, i, visited)
		return temp
	
    # method to perform DFS (Depth First Search)
    # and store the vertices in a list
	def _DFSUtil(self, v : int, visited : List[bool]) -> List[int]:
		output : List[int] = []
		stack : List[int] = []
		stack.append(v)
		while(len(stack) > 0):
			# pop vertex from stack
			v = stack.pop()
			if not visited[v]:
				# Mark the current vertex as visited
				visited[v] = True
				# Store the vertex to list
				output.append(v)
			# Repeat for all vertices adjacent
			# to this vertex v
			for i in self.adj[v]:
				if not visited[i]:
					# Update the list
					stack.append(i)
		return output
	
	# method to add an undirected edge
	# w and v are the indices of the vertices
    # connected by the edge
	def add_edge(self, v : int, w : int) -> None:
		self.adj[v].append(w)
		self.adj[w].append(v)

	# computes the connected components of the graph
    # returns a list of lists (of nodes IDs), where each list is a connected component
	def connected_components(self) -> List[List[int]]:
		visited = [False]*self.V
		cc : List[List[int]] = []
		for v in range(self.V):
			if not visited[v]:
				cc.append(self._DFSUtil(v, visited))
		return cc
