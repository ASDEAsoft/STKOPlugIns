import importlib
import tclimport.parser
import tclimport.builder
importlib.reload(tclimport.parser)
importlib.reload(tclimport.builder)
from tclimport.parser import parse_tcl
from tclimport.builder import build_cae

# the tcl parsed document
class tcldoc:
	
	# constructor
	def __init__(self):
		# nodes
		self.nodes = {}
		# limit curves
		self.limit_curves = {}
		# geom trans
		self.geom_trans = {}
		# elements
		self.elements = {}
		# links
		self.links = {}
		# materials uniaxial
		self.materials_1d = {}
		# fiber sections
		self.fibers = {}
		# sections
		self.sections = {}
		# aggregators
		self.aggregators = {}
		# beam properties
		self.beam_props = {}
		# plate_sections
		self.plate_sections = {}
		# link_materials
		self.link_materials = {}
		# regions
		self.regions = {}
		# recorders
		self.recorders = []
		# stages
		self.stages = []
	
	# main parsing function
	def parse(self, filecontents):
		parse_tcl(self, filecontents)
	
	# build CAE model
	def build(self, callbacks):
		build_cae(self, callbacks)
	
	# return a string representation of self
	def __str__(self):
		return (
			'TCL DOC {{\n'
			'\tNODES {{\n{}\n\t}}\n'
			'\tLIMIT_CURVES {{\n{}\n\t}}\n'
			'\tGEOM_TRANS {{\n{}\n\t}}\n'
			'\tELEMENTS {{\n{}\n\t}}\n'
			'\tLINKS {{\n{}\n\t}}\n'
			'\tUNIAXIAL_MATERIALS {{\n{}\n\t}}\n'
			'\tFIBERS {{\n{}\n\t}}\n'
			'\tSECTIONS {{\n{}\n\t}}\n'
			'\tAGGREGATORS {{\n{}\n\t}}\n'
			'\tBEAM_PROPERTIES {{\n{}\n\t}}\n'
			'\tPLATE_SECTIONS {{\n{}\n\t}}\n'
			'\tLINK_MATERIALS {{\n{}\n\t}}\n'
			'\tREGIONS {{\n{}\n\t}}\n'
			'\tRECORDERS {{\n{}\n\t}}\n'
			'\tSTAGES {{\n{}\n\t}}\n'
			'}}\n'
			).format(
				'\n'.join(['\t\t{}'.format(i) for i in self.nodes.values()]),
				'\n'.join(['\t\t{}'.format(i) for i in self.limit_curves.values()]),
				'\n'.join(['\t\t{}'.format(i) for i in self.geom_trans.values()]),
				'\n'.join(['\t\t{}'.format(i) for i in self.elements.values()]),
				'\n'.join(['\t\t{}'.format(i) for i in self.links.values()]),
				'\n'.join(['\t\t{}'.format(i) for i in self.materials_1d.values()]),
				'\n'.join(['\t\t{}'.format(i) for i in self.fibers.values()]),
				'\n'.join(['\t\t{}'.format(i) for i in self.sections.values()]),
				'\n'.join(['\t\t{}'.format(i) for i in self.aggregators.values()]),
				'\n'.join(['\t\t{}'.format(i) for i in self.beam_props.values()]),
				'\n'.join(['\t\t{}'.format(i) for i in self.plate_sections.values()]),
				'\n'.join(['\t\t{}'.format(i) for i in self.link_materials.values()]),
				'\n'.join(['\t\t{}'.format(i) for i in self.regions.values()]),
				'\n'.join(['\t\t{}'.format(i) for i in self.recorders]),
				'\n'.join(['\t\t{}'.format(i) for i in self.stages])
				)