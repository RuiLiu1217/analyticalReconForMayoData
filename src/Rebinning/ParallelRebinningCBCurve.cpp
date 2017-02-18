/////////////////////////////////////////////////////
// Imaging and Informatics Lab
// Department of Electrical and Computer Engineering
// University of Massachusetts Lowell
// Graduate School of Arts and Sciences
// Wake Forest University
// \brief Rebinning the projection data
// \author Rui Liu, Hengyong Yu
// \date Jan. 15, 2016
// \version 1.0
//////////////////////////////////////////////////////
#include <mex.h>
#include "ParallelRebinningCBCurve.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // The projection is ordered in ZL, YL, ViewN order
    float* Proj = (float*)mxGetPr(prhs[0]);
    int YL = *((int*)(mxGetPr(prhs[1])));
    int ZL = *((int*)(mxGetPr(prhs[2])));
    int ViewN = *((int*)(mxGetPr(prhs[3])));
    float YLC = *((float*)(mxGetPr(prhs[4])));
    float dYA = *((float*)(mxGetPr(prhs[5])));
    float DeltaTheta = *((float*)(mxGetPr(prhs[6])));
    float PLC = *((float*)(mxGetPr(prhs[7])));
    float DeltaT = *((float*)(mxGetPr(prhs[8])));
    float DeltaFai = *((float*)(mxGetPr(prhs[9])));
    float SO = *((float*)(mxGetPr(prhs[10])));
    int useGPU = *((int*)(mxGetPr(prhs[11])));
    
    const mwSize dims[]={static_cast<mwSize>(ZL),static_cast<mwSize>(YL),static_cast<mwSize>(ViewN)};
    plhs[0] = mxCreateNumericArray(3,dims,mxSINGLE_CLASS,mxREAL);
    float* outputProj = (float*)mxGetPr(plhs[0]);
    if(useGPU == 1)
    {
          ParallelRebinningCBCurve_GPU(outputProj, Proj,
              YL, ZL, ViewN, YLC - 1.0f, dYA, DeltaTheta,
              PLC - 1.0f, DeltaT, DeltaFai, SO);
    }
    else
    {
          ParallelRebinningCBCurve_CPU(outputProj, Proj,
             YL, ZL, ViewN, YLC - 1.0f, dYA, DeltaTheta,
             PLC - 1.0f, DeltaT, DeltaFai, SO);
    }
}
