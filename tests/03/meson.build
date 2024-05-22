tests = executable(
	'tests_3',
	'test_verlet.cpp',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, google_test, mpi]
)

test('Tests3', tests)
