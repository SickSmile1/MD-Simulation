eigen = dependency('eigen3', version : '>=3.4.0')
benchmark = dependency('benchmark')

executable(
  'bench',
  'gemv_benchmarks.cpp',
  dependencies : [eigen, benchmark]
)
fs = import('fs')
fs.copyfile('plot_benchmark.py')

executable(
  'bench2',
  'equilibrium_benchmarks.cpp',
  include_directories : [lib_incdirs],
  link_with : [lib],
  dependencies : [eigen, benchmark]
)