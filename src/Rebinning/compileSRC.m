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

curveCPUFile = 'ParallelRebinningCBCurve_CPU.cpp';
curveCPUOFile = 'ParallelRebinningCBCurve_CPU.o';

curveMEXFile = 'ParallelRebinningCBCurve.cpp';
curveGPUFile = 'ParallelRebinningCBCurve_GPU.cu';
curveGPUOFile = 'ParallelRebinningCBCurve_GPU.o';
curveSys = [NVCC, ' -Xcompiler -O3 --use_fast_math --compile -o ', curveGPUOFile, '  --compiler-options -fPIC -I"', MATLAB_INC, '" -I"', CUDA_INC, '" -I"', CUDA_SAMP_INC, '" ', curveGPUFile];
system(curveSys);
!g++ -std=c++11 -fopenmp -lgomp -c -fPIC ParallelRebinningCBCurve_CPU.cpp';


if strcmp(arch, 'GLNXA64') % LINUX
    mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -fopenmp -Wall -std=c++11" -L"/usr/local/cuda/lib64" -lcudart -lgomp ParallelRebinningCBCurve.cpp ParallelRebinningCBCurve_GPU.o ParallelRebinningCBCurve_CPU.o ;
elseif strcmp(arch, 'PCWIN64') % WINDOWS
    mex -v -largeArrayDims COMPFLAGS="$COMPFLAGS -Wall -std=c++11" -L"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v8.0\lib\x64" -lcudart -lParallelRebinningCBCurve_call ParallelRebinningCBCurve.cpp;
elseif strcmp(arch, 'MACI64')
    
end
