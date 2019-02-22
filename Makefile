CCONST=
use_openmp=
MKDIR_P=mkdir -p
dir_guard=@$(MKDIR_P) $(@D)
DOXY=doxygen
srcdir=src
builddir=build
incdir=include
productdir=for_MATLAB
testproductdir=test

c_compiler=clang
cpp_compiler=clang++
dist_dir=distribution
source_dist_dir=source_distribution
extraCFLAGS=-O3
extraCXXFLAGS=-O3 -std=c++11 -stdlib=libc++ -I/usr/local/include
extraLDFLAGS=-L/usr/local/lib -larpack -O3

-include localdefs.mk
ifndef MATLAB_ROOT
$(error MATLAB_ROOT is not defined in localdefs.mk. Please point it to your MATLAB installation, e.g. /Applications/MATLAB_R2016a.app)
endif
ifdef MATLABPATH
$(error MATLABPATH should not be defined in localdefs.mk)
endif


MATLABPATH=$(MATLAB_ROOT)/bin

MEXSUFFIX=$(shell $(MATLABPATH)/mexext)

_DEPS = spikeTrainFun.hpp tSNEExceptions.hpp tSNERunner.hpp tSNERunnerD.hpp tSNERunnerV.hpp ssimsHelper.hpp
DEPS = $(patsubst %,$(incdir)/%,$(_DEPS))


_OBJ_SPK = spikeTrainFun.o
OBJ_SPK = $(patsubst %,$(builddir)/%,$(_OBJ_SPK))

_OBJ_tSNE = tSNEExceptions.o tSNERunner.o tSNERunnerD.o tSNERunnerV.o
OBJ_tSNE = $(patsubst %,$(builddir)/%,$(_OBJ_tSNE))

_OBJ_SSIMS = ssimsHelper.o
OBJ_SSIMS = $(patsubst %,$(builddir)/%,$(_OBJ_SSIMS))


_OBJ_EV = randperm.o
OBJ_EV = $(patsubst %,$(builddir)/%,$(_OBJ_EV))

_OBJ_KNN = knnHelpers.o
OBJ_KNN = $(patsubst %,$(builddir)/%,$(_OBJ_KNN))


_DEPS_TEST = example_data_helper.hpp
DEPS_TEST = $(patsubst %,$(incdir)/%,$(_DEPS_TEST))

_OBJ_TEST = example_data_helper.o
OBJ_TEST = $(patsubst %,$(builddir)/%,$(_OBJ_TEST))

_SSIMS_MEX_PRODUCTS = getSSIMDMatBetweenTimePoints.$(MEXSUFFIX) getSSIMS.$(MEXSUFFIX) \
	runtSNE.$(MEXSUFFIX) getSSIMSprojection.$(MEXSUFFIX) SSIMSToolboxTestParallel.$(MEXSUFFIX) \
	spikeTrainEventWindowCounts.$(MEXSUFFIX) spikeTrainEventWindowSpikes.$(MEXSUFFIX) getSSIMSSweep.$(MEXSUFFIX) \
	vpSpikeTimeDist.$(MEXSUFFIX)
SSIMS_MEX_PRODUCTS = $(patsubst %,$(productdir)/%,$(_SSIMS_MEX_PRODUCTS))

_MAT_TEST_PRODUCTS = SSIMS_test spikeTrain_test
MAT_TEST_PRODUCTS = $(patsubst %,$(testproductdir)/%,$(_MAT_TEST_PRODUCTS))

_LEGACY_MEX_PRODUCTS = KNearestNeighbor_classify.$(MEXSUFFIX)
LEGACY_MEX_PRODUCTS = $(patsubst %,$(productdir)/%,$(_LEGACY_MEX_PRODUCTS))

.PHONY: clean docs test all distribution

CFLAGS=-I$(incdir) -Wall $(extraCFLAGS)

ifneq ($(strip $(use_openmp)),)
CXXFLAGS=-fopenmp=libiomp5 -I$(incdir) -Wall $(extraCXXFLAGS)
else
CXXFLAGS=-I$(incdir) -Wall $(extraCXXFLAGS)
endif

LDFLAGS=-larmadillo $(extraLDFLAGS)
GPP=$(cpp_compiler) $(CXXFLAGS)
GCC=$(c_compiler) $(CFLAGS)

MEX=$(MATLABPATH)/mex -largeArrayDims -lmwblas -larmadillo CFLAGS='$(CFLAGS)' CXXFLAGS='$(CXXFLAGS)' LDFLAGS='$(LDFLAGS) \$$LDFLAGS'

DYLD_LIBRARY_PATH_TEST=$(MATLAB_ROOT)/bin/maci64:$(MATLAB_ROOT)/sys/os/maci64:$(DYLD_LIBRARY_PATH)

all: $(SSIMS_MEX_PRODUCTS) $(LEGACY_MEX_PRODUCTS) $(MAT_TEST_PRODUCTS)

$(builddir)/%.o: $(srcdir)/%.cpp $(DEPS)
	$(dir_guard)
	$(GPP) -c $< -o $@
$(builddir)/%.o: $(srcdir)/%.c
	$(dir_guard)
	$(GCC) -c $< -o $@
$(builddir)/%.o: $(srcdir)/test/%.cpp $(DEPS) $(DEPS_TEST)
	$(dir_guard)
	$(MEX) -c $< -outdir $(builddir)
$(productdir)/getSSIMDMatBetweenTimePoints.$(MEXSUFFIX): $(srcdir)/getSSIMDMatBetweenTimePoints.cpp $(OBJ_SPK)
	$(MEX) $(CCONST) $^ -output $@
$(productdir)/getSSIMDMatBetweenNeurons.$(MEXSUFFIX): $(srcdir)/getSSIMDMatBetweenNeurons.cpp $(OBJ_SPK)
	$(MEX) $(CCONST) $^ -output $@
$(productdir)/getSSIMS.$(MEXSUFFIX): $(srcdir)/getSSIMS.cpp $(OBJ_SPK) $(OBJ_tSNE)
	$(MEX) $(CCONST) $^ -output $@
$(productdir)/getSSIMSSweep.$(MEXSUFFIX): $(srcdir)/getSSIMSSweep.cpp $(OBJ_SPK) $(OBJ_tSNE) $(OBJ_SSIMS)
	$(MEX) $(CCONST) $^ -output $@
$(productdir)/getSSIMSprojection.$(MEXSUFFIX): $(srcdir)/getSSIMSprojection.cpp $(OBJ_SPK)
	$(MEX) $(CCONST) $^ -output $@
$(productdir)/runtSNE.$(MEXSUFFIX): $(srcdir)/runtSNE.cpp $(OBJ_tSNE)
	$(MEX) $(CCONST) $^ -output $@
$(productdir)/SSIMSToolboxTestParallel.$(MEXSUFFIX): $(srcdir)/SSIMSToolboxTestParallel.cpp
	$(MEX) $(CCONST) $^ -output $@
$(productdir)/spikeTrainEventWindowCounts.$(MEXSUFFIX): $(srcdir)/spikeTrainEventWindowCounts.cpp $(OBJ_SPK)
	$(MEX) $(CCONST) $^ -output $@
$(productdir)/spikeTrainEventWindowSpikes.$(MEXSUFFIX): $(srcdir)/spikeTrainEventWindowSpikes.cpp $(OBJ_SPK)
	$(MEX) $(CCONST) $^ -output $@
$(productdir)/vpSpikeTimeDist.$(MEXSUFFIX): $(srcdir)/vpSpikeTimeDist.cpp $(OBJ_SPK)
	$(MEX) $(CCONST) $^ -output $@
$(testproductdir)/SSIMS_test: $(srcdir)/test/SSIMS_test.cpp $(OBJ_SPK) $(OBJ_tSNE) $(OBJ_TEST) $(OBJ_EV) $(OBJ_KNN)
	$(dir_guard)
	$(MEX) -client engine $^ -output $@
$(testproductdir)/spikeTrain_test: $(srcdir)/test/spikeTrain_test.cpp $(OBJ_SPK)  $(OBJ_TEST)
	$(dir_guard)
	$(MEX) -client engine $^ -output $@

$(productdir)/KNearestNeighbor_classify.$(MEXSUFFIX): $(srcdir)/KNearestNeighbor_classify.c $(OBJ_EV) $(OBJ_KNN)
	$(MEX) $(CCONST) $^ -output $@


test: $(MAT_TEST_PRODUCTS)
	DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH_TEST) test/spikeTrain_test examples/SSIMS_demo_data_center_out.mat
	DYLD_LIBRARY_PATH=$(DYLD_LIBRARY_PATH_TEST) test/SSIMS_test examples/SSIMS_demo_data_center_out.mat

docs: Doxyfile
	$(DOXY) Doxyfile

distribution: $(SSIMS_MEX_PRODUCTS) $(LEGACY_MEX_PRODUCTS)
	$(MKDIR_P) $(dist_dir)
	cp -LR $(productdir)/* $(dist_dir)

source_distribution: $(srcdir)/*.cpp $(DEPS) $(srcdir)/*.c
	$(MKDIR_P) $(source_dist_dir)
	$(MKDIR_P) $(source_dist_dir)/$(productdir)
	cp -LR $(incdir) $(source_dist_dir)
	cp -LR $(srcdir) $(source_dist_dir)
	rsync -Lr --exclude='html' --exclude='latex' doc $(source_dist_dir)
	rsync -Lr --exclude='*.mex*' $(productdir) $(source_dist_dir)
	cp build_files_win.m Doxyfile license.txt localdefs_sample.mk Makefile README $(source_dist_dir)

clean:
	rm -f $(builddir)/*.o
	rm -f $(builddir)/*.obj
	rm -f $(MAT_TEST_PRODUCTS)
	rm -f $(LEGACY_MEX_PRODUCTS)
	rm -f $(SSIMS_MEX_PRODUCTS)
	rm -fr doc/html
	rm -fr doc/latex
	rm -rf $(dist_dir)
	rm -rf $(source_dist_dir)

