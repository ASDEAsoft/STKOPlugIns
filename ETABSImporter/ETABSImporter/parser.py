from ETABSImporter.document import document, vertex, frame, area
from collections import defaultdict

'''
section modifiers parse + STKO
'''

class parser:
	
	def __init__(self, fname):
		# initialize document
		self.doc = document()
		# initialize command map
		self.commands = defaultdict(list)
		# build command map
		with open(fname, 'r') as ifile:
			lines = [pline for pline in (line.strip() for line in ifile.readlines()) if (pline and not pline.startswith('#'))]
			value = None
			for line in lines:
				if line.startswith('*'):
					value = self.commands[line]
				else:
					if value is not None:
						value.append(line)
		# parse
		self._parse()
		print(self.doc)
	
	def _parse(self):
		self._parse_nodes()
		self._parse_frames()
		self._parse_areas()
	
	def _parse_nodes(self):
		for item in self.commands['* JOINT_COORDINATES']:
			words = item.split(',')
			id = int(words[0])
			x = float(words[1])
			y = float(words[2])
			z = float(words[3])
			self.doc.vertices[id] = vertex(x,y,z)
	
	def _parse_frames(self):
		for item in self.commands['* FRAMES_CONNECTIVITY']:
			words = item.split(',')
			id = int(words[0])
			self.doc.frames[id] = frame([int(words[i]) for i in range(1, 3)])
	
	def _parse_areas(self):
		for item in self.commands['* AREAS_CONNECTIVITY']:
			words = item.split(',')
			id = int(words[0])
			self.doc.areas[id] = area([int(words[i]) for i in range(1, 5)])