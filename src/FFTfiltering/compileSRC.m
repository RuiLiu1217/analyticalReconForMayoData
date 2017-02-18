%Compile ParallelRebinningCBCurve
%------------------------------------------------------------------------- 
% Wake Forest University Health Sciences
% Univeristy of Massachusetts Lowell
% Date: 2016-08-04
% Routine: compileCPP.m
% Author
%	Rui Liu
% Email: liurui1217@gmail.com
%-------------------------------------------------------------------------
arch = computer;
if strcmp(arch, 'GLNXA64') % LINUX
    NVCC = '/usr/local/cuda/bin/nvcc';
    MATLAB_INC = [matlabroot, '/extern/include'];
    CUDA_INC = '/usr/local/cuda/include';
    CUDA_SAMP_INC = '/usr/local/cuda/samples/common/inc';
elseif  strcmp(arch,'PCWIN64') % WINDOWS
    CUDA_PATH = getenv('CUDA_PATH');
    NVCC = ['"',CUDA_PATH,'\bin\nvcc.exe"'];
    MATLAB_INC = [getenv('MATLAB_PATH'),'\extern\include'];
    CUDA_INC = [getenv('CUDA_PATH'),'\include'];
    CUDA_LIB = ['"',getenv('CUDA_PATH'),'\lib\x64"'];
    CUDA_SAMP_INC = [getenv('NVCUDASAMPLES_ROOT'),'\common\inc'];
elseif strcmp(arch, 'MACI64') % MAC (IS NOT IMPLEMENTED)

end


filterMEXFile = 'reWeigAdFiltr_GPU.cpp';
filterGPUFile = 'filtering.cu';
filterGPUOFile = 'filtering.o';

if strcmp(arch, 'GLNXA64') % LINUX
    filterSys = [NVCC, ' -Xcompiler -O3 --use_fast_math --compile -o ', filterGPUOFile, '  --compiler-options -fPIC -I"', MATLAB_INC, '" -I"', CUDA_INC, '" -I"', CUDA_SAMP_INC, '" ', filterGPUFile];
    system(filterSys);
    mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -fopenmp -Wall -std=c++11" -L"/usr/local/cuda/lib64" -lcudart -lcufft  reWeigAdFiltr_GPU.cpp filtering.o;
elseif strcmp(arch, 'PCWIN64') % WINDOWS
    % One GPU
    mex -v -largeArrayDims COMPFLAGS="$COMPFLAGS -Wall -std=c++11" -L"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v8.0\lib\x64" -lcudart -lcufft -lFFTfiltering reWeigAdFiltr_GPU.cpp;
    
    
    % Rui Liu's Note: This is the link to the multi-GPU based FFT, However, one GPU is
    % GTX670. The behavior of this GPU is wrong, Tesla K10 and Titan X
    % performs good. Newer version will be added to address problem later.
    % mex -v -largeArrayDims COMPFLAGS="$COMPFLAGS -Wall -std=c++11" -L"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v8.0\lib\x64" -lcudart -lcufft -lmultiGPU_FFTfiltering reWeigAdFiltr_GPU.cpp;
elseif strcmp(arch, 'MACI64')
    
end
