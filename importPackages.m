function importPackages()
%IMPORTPACKAGES Import the necessary packages. Change the path variables
%here!

CAIMAN_PATH = 'D:\Codebase\CaImAn-MATLAB\';
CELLSORT_PATH = 'D:\Codebase\CellSort\';
NORMCORRE_PATH = 'D:\Codebase\NoRMCorre\';
MATLAB_2P_PATH = 'D:\Codebase\matlab-2p\';

addpath(genpath(CAIMAN_PATH)); %CaImAn replacing ca_source_extraction
addpath(genpath(CELLSORT_PATH));
addpath(genpath(NORMCORRE_PATH));
addpath(genpath(MATLAB_2P_PATH));
end

