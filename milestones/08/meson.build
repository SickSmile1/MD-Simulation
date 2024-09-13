executable(
	'milestone08',
	'main.cpp',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, mpi]
)

fs = import('fs')
fs.copyfile('cluster_923.xyz')
fs.copyfile('cluster_3871.xyz')
