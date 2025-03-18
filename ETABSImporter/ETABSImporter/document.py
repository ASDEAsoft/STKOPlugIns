class vertex:
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z
	def __str__(self):
		return '({:8.3g}, {:8.3g},{:8.3g})'.format(self.x, self.y, self.z)
	def __repr__(self):
		return self.__str__()

class frame:
	def __init__(self, nodes):
		self.nodes = nodes
	def __str__(self):
		return '({:8},{:8})'.format(*self.nodes)
	def __repr__(self):
		return self.__str__()

class area:
	def __init__(self, nodes):
		self.nodes = nodes
	def __str__(self):
		return '({:8},{:8},{:8},{:8})'.format(*self.nodes)
	def __repr__(self):
		return self.__str__()

class document:
	def __init__(self):
		self.vertices = {}
		self.frames = {}
		self.areas = {}
	def __str__(self):
		from io import StringIO
		f = StringIO()
		f.write('Vertices\n')
		for i,v in self.vertices.items():
			f.write(f"{i:8} : {v}\n")
		f.write('Frames\n')
		for i,v in self.frames.items():
			f.write(f"{i:8} : {v}\n")
		f.write('Areas\n')
		for i,v in self.areas.items():
			f.write(f"{i:8} : {v}\n")
		return f.getvalue()
	def __repr__(self):
		return self.__str__()
	
	def plot(self):
		import matplotlib.pyplot as plt
		from shapely.geometry import Polygon, LineString
		from shapely.ops import split
		from descartes import PolygonPatch  # For plotting polygons
		print('done')