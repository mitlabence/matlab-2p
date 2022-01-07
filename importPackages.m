function importPackages()
%IMPORTPACKAGES Import the necessary packages. Change the path variables
%here!

CAIMAN_PATH = 'D:\Software\CaImAn\CaImAn-MATLAB\';
CELLSORT_PATH = 'D:\Software\CellSort\CellSort\';
NORMCORRE_PATH = 'D:\Software\NoRMCorre\NoRMCorre\';
OTHER_EXTERNAL_PATH = 'D:\PhD\matlab-2p\matlab-2p\external\other\';
MARTIN_FUNCTIONS_PATH = 'D:\PhD\matlab-2p\matlab-2p\external\Martin Code\'; 
FILEHANDLING_PATH = 'D:\PhD\matlab-2p\matlab-2p\FileHandling\';
PROCESSING_PATH = 'D:\PhD\matlab-2p\matlab-2p\Processing\';


addpath(genpath(CAIMAN_PATH)); %CaImAn replacing ca_source_extraction
addpath(genpath(CELLSORT_PATH));
addpath(genpath(NORMCORRE_PATH));
addpath(genpath(OTHER_EXTERNAL_PATH));
addpath(genpath(MARTIN_FUNCTIONS_PATH));
addpath(genpath(FILEHANDLING_PATH));
addpath(genpath(PROCESSING_PATH));
end

