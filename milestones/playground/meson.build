executable(
	'playground',
	'main.cpp',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, mpi]
)

executable(
	'mackay',
	'ih.C',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, mpi]
)

fs = import('fs')
fs.copyfile('whisker_small.xyz')
