import sys
sys.path.append('C:/Develop/Utility/STKO-install/data/win64')

from ETABSImporter.parser import parser
#p = parser('Examples/Ex00/Model_Info.txt')
p = parser('Examples/Ex01/STKO_IMPORT_40_Story_Rc_Core_Wall_Bldng_e4.3.txt')


#p.doc.split_area_members()
print(p.doc)

p.doc.plot(vertices=False, frames=True, areas=True)

'''
Assumptions
areas are generic flat polygons
- they can be split in multiple lines
- holes use another command *AREA_OPENINGS_CONNECTIVITY
'''