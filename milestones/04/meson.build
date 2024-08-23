subdir('plot')
executable(
	'milestone04',
	'main.cpp',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, mpi, fmt_dep]
)

fs = import('fs')
fs.copyfile('lj54.xyz')