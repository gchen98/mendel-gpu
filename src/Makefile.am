# This file is part of mendel-gpu.

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

AM_CPPFLAGS = @USE_GPU@ @OPENCL_INCLUDE@

bin_PROGRAMS = impute

impute_SOURCES = \
	init/read_pen.cpp \
	init/read_pen_opencl.cpp \
	init/common.cpp \
	init/common_opencl.cpp \
	init/denovo.cpp \
	init/denovo_glf.cpp \
	init/denovo_glf_opencl.cpp \
	init/denovo_opencl.cpp \
	init/denovo_reads.cpp \
	init/denovo_reads_opencl.cpp \
	init/guide.cpp \
	init/guide_glf.cpp \
	init/guide_glf_opencl.cpp \
	init/guide_opencl.cpp \
	init/guide_reads.cpp \
	init/guide_reads_opencl.cpp \
	penetrance/common.cpp \
	penetrance/common_opencl.cpp \
	penetrance/denovo_glf.cpp \
	penetrance/denovo_glf_opencl.cpp \
	penetrance/denovo_reads.cpp \
	penetrance/guide_glf.cpp \
	penetrance/guide_glf_opencl.cpp \
	penetrance/guide_reads.cpp \
	penetrance/read_parser.cpp \
	penetrance/read_pen.cpp \
	penetrance/read_pen_opencl.cpp \
	imputation/denovo.cpp \
	imputation/denovo_opencl.cpp \
	imputation/guide.cpp \
	imputation/guide_opencl.cpp \
	main.cpp \
	io_manager.cpp \
	cl_templates.cpp \
	$(NULL)

impute_LDADD = \
	-lbam @OPENCL_LIB@ \
	-lstdc++ -lz -lpthread \
	$(NULL)
