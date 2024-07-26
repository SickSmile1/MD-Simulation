incdir = include_directories('01','02','03','04','05','06','07')

test_sources = [
    '01/test_hello_world.cpp',
    '01/gtest_mpi_main.cpp',
    '02/test_pi.cpp',
    '03/test_verlet.cpp',
    '04/test_lj_direct_summation.cpp',
    '05/test_cubic_latice.cpp',
    '05/test_berendsten_thermostat.cpp',
    '06/test_neighbors.cpp',
    '06/test_updated_lj.cpp',
    '07/test_ducastelle.cpp',
    '07/test_energy_conservation.cpp'
]


tests = executable(
    'tests2',
    test_sources,
    include_directories : [incdir,lib_incdirs],
    link_with : [lib],
    dependencies : [eigen,mpi,google_test],
)

incdir2 = include_directories('09')
test_sources2 = [
    '09/test_whisker_energies_temps.cpp'
]
tests = executable(
    'tests3',
    test_sources2,
    include_directories : [incdir2,lib_incdirs],
    link_with : [lib],
    dependencies : [eigen,mpi,google_test],
)
# test('tests', tests)

fs = import('fs')
fs.copyfile('07/cluster_923.xyz')
fs.copyfile('09/whisker_small.xyz')