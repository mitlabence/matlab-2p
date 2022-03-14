function importPackages()
%IMPORTPACKAGES Import the necessary packages. Change the path variables
%here!

CAIMAN_PATH = 'D:\Software\CaImAn\CaImAn-MATLAB\';
CELLSORT_PATH = 'D:\Software\CellSort\CellSort\';
NORMCORRE_PATH = 'D:\Software\NoRMCorre\NoRMCorre\';
MATLAB_2P_PATH = 'D:\PhD\matlab-2p\';

addpath(genpath(CAIMAN_PATH)); %CaImAn replacing ca_source_extraction
addpath(genpath(CELLSORT_PATH));
addpath(genpath(NORMCORRE_PATH));
addpath(genpath(MATLAB_2P_PATH));
end

