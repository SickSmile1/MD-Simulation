incdir = include_directories('01','02','03','04')

test_sources = [
    '01/test_hello_world.cpp',
    '01/gtest_mpi_main.cpp',
    '02/test_pi.cpp',
    '03/test_verlet.cpp',
    '04/test_lj_direct_summation.cpp'
]

tests = executable(
    'tests2',
    test_sources,
    include_directories : [incdir,lib_incdirs],
    link_with : [lib],
    dependencies : [eigen,mpi,google_test],
)

test('tests', tests)
