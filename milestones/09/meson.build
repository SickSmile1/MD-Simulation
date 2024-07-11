executable(
	'milestone09',
	'main.cpp',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, mpi]
)

fs = import('fs')
fs.copyfile('whisker_large.xyz')
fs.copyfile('whisker_small.xyz')