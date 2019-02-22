% Run this file from its location in MATLAB. MEX functions should be
% compiled and if everything goes right, they will be copied to the
% 'for_MATLAB' folder, which should be included on MATLAB's path.

% Add options you might want to pass to the compiler to this cell array:
mexopts = {};
% build files in this directory:
build_dir_name = 'build';

% Where does armadillo's include directory live (absolute, or relative to this file):
armadillo_include_dir = fullfile('lib', 'armadillo', 'include');

% Additional flags to add to g++:
cxxflags = ['CXXFLAGS="-std=c++11 "'];


cflags = 'CFLAGS="-Wall -O3 -std=c11"';

% Compilation with MS Visual C++ might require the explicit locations of
% BLAS and LAPACK; just uncomment the following lines:

% blaslib = fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft',...
% 'libmwblas.lib');
% lapacklib = fullfile(matlabroot,'extern','lib',computer('arch'),'microsoft',...
%   'libmwlapack.lib');
% mexopts = [mexopts, blaslib, lapacklib];


% Move on, nothing to see below this line.

mexopts = [mexopts, '-outdir', build_dir_name];
armadillo_include = ['-I', armadillo_include_dir];
[~, ~, mkdirmesstr] = mkdir(build_dir_name);
if ~(strcmp(mkdirmesstr,'MATLAB:MKDIR:DirectoryExists')|| strcmp(build_dir_name,''))
	error('Directory %s couldn''t be created.\n', build_dir_name);
end
%%
mex('-c', '-largeArrayDims', mexopts{:}, '-Iinclude', armadillo_include, fullfile('src', 'spikeTrainFun.cpp'), 'CXXFLAGS="-std=c++11"')
mex('-largeArrayDims', mexopts{:}, '-Iinclude', armadillo_include, fullfile('src', 'spikeTrainEventWindowSpikes.cpp'), fullfile(build_dir_name, 'spikeTrainFun.obj'), 'CXXFLAGS="-std=c++0x"')
mex('-largeArrayDims', mexopts{:}, '-Iinclude', armadillo_include, fullfile('src', 'spikeTrainEventWindowCounts.cpp'), fullfile(build_dir_name, 'spikeTrainFun.obj'), 'CXXFLAGS="-std=c++0x"')
mex('-largeArrayDims', mexopts{:}, '-Iinclude', armadillo_include, fullfile('src', 'vpSpikeTimeDist.cpp'), fullfile(build_dir_name, 'spikeTrainFun.obj'), 'CXXFLAGS="-std=c++0x"')
mex('-largeArrayDims', mexopts{:}, '-Iinclude', armadillo_include, fullfile('src', 'SSIMSToolboxTestParallel.cpp'), 'CXXFLAGS="-std=c++0x"')
mex('-largeArrayDims', mexopts{:}, '-Iinclude', armadillo_include, fullfile('src', 'getSSIMDMatBetweenTimePoints.cpp'), fullfile(build_dir_name, 'spikeTrainFun.obj'), 'CXXFLAGS="-std=c++0x"')

mex('-largeArrayDims', mexopts{:}, ...
    '-Iinclude', armadillo_include, fullfile('src', 'getSSIMSprojection.cpp'), fullfile(build_dir_name, 'spikeTrainFun.obj'), cxxflags)
mex('-c', mexopts{:}, '-largeArrayDims', '-Iinclude', armadillo_include, fullfile('src', 'ssimsHelper.cpp'), cxxflags)


mex('-c', mexopts{:}, '-largeArrayDims', '-Iinclude', armadillo_include, fullfile('src', 'tSNEExceptions.cpp'), cxxflags)
mex('-c', mexopts{:}, '-largeArrayDims', '-Iinclude', armadillo_include, fullfile('src', 'tSNERunner.cpp'), cxxflags)
mex('-c', mexopts{:}, '-largeArrayDims', '-Iinclude', armadillo_include, fullfile('src', 'tSNERunnerD.cpp'), cxxflags)
mex('-c', mexopts{:}, '-largeArrayDims', '-Iinclude', armadillo_include, fullfile('src', 'tSNERunnerV.cpp'), cxxflags)


mex('-largeArrayDims', mexopts{:}, '-Iinclude', armadillo_include, fullfile('src', 'getSSIMS.cpp'), ...
    fullfile(build_dir_name, 'spikeTrainFun.obj'), fullfile(build_dir_name, 'tSNEExceptions.obj'), fullfile(build_dir_name, 'tSNERunner.obj'), ...
    fullfile(build_dir_name, 'tSNERunnerD.obj'), fullfile(build_dir_name, 'tSNERunnerV.obj'), ...
    cxxflags)

mex('-largeArrayDims', mexopts{:}, '-Iinclude', armadillo_include, fullfile('src', 'runtSNE.cpp'), ...
    fullfile(build_dir_name, 'tSNEExceptions.obj'), fullfile(build_dir_name, 'tSNERunner.obj'), ...
    fullfile(build_dir_name, 'tSNERunnerD.obj'), fullfile(build_dir_name, 'tSNERunnerV.obj'), ...
    cxxflags)

mex('-largeArrayDims', mexopts{:}, '-Iinclude', armadillo_include, fullfile('src', 'getSSIMSSweep.cpp'), ...
    fullfile(build_dir_name, 'spikeTrainFun.obj'), fullfile(build_dir_name, 'tSNEExceptions.obj'), fullfile(build_dir_name, 'tSNERunner.obj'), ...
    fullfile(build_dir_name, 'tSNERunnerD.obj'), fullfile(build_dir_name, 'tSNERunnerV.obj'), fullfile(build_dir_name, 'ssimsHelper.obj'), ...
    cxxflags)
%%
mex('-c', mexopts{:}, '-largeArrayDims', '-Iinclude', fullfile('src', 'knnHelpers.c'), cflags)
mex('-c', mexopts{:}, '-largeArrayDims', '-Iinclude', fullfile('src', 'randPerm.c'), cflags)
mex('-largeArrayDims', mexopts{:}, '-Iinclude', fullfile('src', 'KNearestNeighbor_classify.c'), ...
    fullfile(build_dir_name, 'knnHelpers.obj'), fullfile(build_dir_name, 'randPerm.obj'), ...
    cflags)

%%
copyfile(fullfile(build_dir_name, ['*.' mexext]), 'for_MATLAB')

