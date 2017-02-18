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
    CUDA_LIB = ['"',getenv('CUDA_PATH'),'\lib64"'];
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
!g++ -std=c++11 -fopenmp -lgomp -c -fPIC ParallelRebinningCBCurve_CPU.cpp;


filterMEXFile = 'reWeigAdFiltr_GPU.cpp';
filterGPUFile = 'filtering.cu';
filterGPUOFile = 'filtering.o';
filterSys = [NVCC, ' -Xcompiler -O3 --use_fast_math --compile -o ', filterGPUOFile, '  --compiler-options -fPIC -I"', MATLAB_INC, '" -I"', CUDA_INC, '" -I"', CUDA_SAMP_INC, '" ', filterGPUFile];
system(filterSys);

% 
% backMEXFile = 'Backprojection3DWeighting_mex.cpp';
% backCPUFile = 'Backprojection3DWeighting_CPU.cpp';
% backGPUFile = 'Backprojection3DWeighting_GPU.cu';
% backGPUOFile = 'Backprojection3DWeighting_GPU.o';
% backSys = [NVCC, ' -Xcompiler -O3 --use_fast_math --compile -o ', backGPUOFile, '  --compiler-options -fPIC -I"', MATLAB_INC, '" -I"', CUDA_INC, '" -I"', CUDA_SAMP_INC, '" ', backGPUFile];
% system(backSys);
% !g++ -std=c++11 -fopenmp -lgomp -c -fPIC Backprojection3DWeighting_CPU.cpp;


% ResultName1 = 'ParallelFDKHelical3DWeightingGPU.o';
% SourceFile1 = 'ParallelFDKHelical3DWeightingGPU.cu';
% CPPFile1 = 'Parallel_FDK_Helical_3DWeighting.cpp';
% sys1 = [NVCC, ' -Xcompiler -O3 --use_fast_math --compile -o ', ResultName1, '  --compiler-options -fPIC -I"', MATLAB_INC, '" -I"', CUDA_INC, '" -I"', CUDA_SAMP_INC, '" ', SourceFile1];
% system(sys1);


if strcmp(arch, 'GLNXA64') % LINUX
    mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -fopenmp -Wall -std=c++11" -L"/usr/local/cuda/lib64" -lcudart -lgomp ParallelRebinningCBCurve.cpp ParallelRebinningCBCurve_GPU.o ParallelRebinningCBCurve_CPU.o ;
    mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -fopenmp -Wall -std=c++11" -L"/usr/local/cuda/lib64" -lcudart -lcufft  reWeigAdFiltr_GPU.cpp filtering.o;
%     mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -fopenmp -Wall -std=c++11" -L"/usr/local/cuda/lib64" -lcudart -lgomp Backprojection3DWeighting_mex.cpp Backprojection3DWeighting_CPU.o Backprojection3DWeighting_GPU.o;
%     mex -v -largeArrayDims  COMPFLAGS="$COMPFLAGS -Wall -std=c++11" -L"/usr/local/cuda/lib64" -lcudart -lgomp Parallel_FDK_Helical_3DWeighting.cpp ParallelFDKHelical3DWeightingGPU.o;
elseif strcmp(arch, 'PCWIN64') % WINDOWS
    
    %mex -v -largeArrayDims COMPFLAGS="$COMPFLAGS -Wall -std=c++11" -L"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v8.0\lib\x64" -lcudart Parallel_FDK_Helical_3DWeighting.cpp ParallelFDKHelical3DWeightingGPU.o;
elseif strcmp(arch, 'MACI64')
    
end
