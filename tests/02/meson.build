tests = executable(
	'tests_2',
	'test_pi.cpp',
	include_directories : [lib_incdirs],
	link_with : [lib],
	dependencies : [eigen, google_test, mpi]
)

test('Tests2', tests)
