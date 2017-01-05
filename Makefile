MAKEFLAGS += --no-builtin-rules                                                                 
.DEFAULT_GOAL := all
.DELETE_ON_ERROR:
.SUFFIXES:
.PHONY: all clean

SOURCES := test-transport.py transport.py flux_divergence.py timestepping.py

all:: build/convergence.svg

clean::
	rm -rf build

build/convergence.svg: \
  convergence.plt \
  build/cf-uniform.dat \
  build/cf-nonuniform.dat \
  build/cf-corr-uniform.dat \
  build/cf-corr-nonuniform.dat
	gnuplot convergence.plt

build/cf-uniform.dat: | build
	pytest test-transport.py -k test_second_order_convergence_on_uniform_mesh_with_cubic_fit

build/cf-nonuniform.dat: | build
	pytest test-transport.py -k test_second_order_convergence_on_nonuniform_mesh_with_cubic_fit
	
build/cf-corr-uniform.dat: | build
	pytest test-transport.py -k test_fourth_order_convergence_on_uniform_mesh_with_corrected_cubic_fit

build/cf-corr-nonuniform.dat: | build
	pytest test-transport.py -k test_second_order_convergence_on_nonuniform_mesh_with_corrected_cubic_fit

build:
	mkdir -p $@
