executable(
	'plot',
	'main.cpp',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, mpi]
)

fs = import('fs')
fs.copyfile('lj54.xyz')
