/**
 * Wake Forest University Health Science
 * University of Massachusetts Lowell
 *
 * Organization:
 * 	Wake Forest University
 *
 * 	reWeiAdFiltr.cpp
 * 	Matlab mex routine for the GPU based reweighting
 * 	and filtering for analytical helical CT
 * 	reconstruction
 *
 * 	Author: Rui Liu
 * 	Email: liurui1217@gmail.com
 * 	Date: 2016-08-04
 *
 * 	Version 1.0
 */

#include "mex.h"
#include "matrix.h"

#include "filtering.h"
#include <iostream>
extern "C"
void filtering(float* hfpwd,
		const float* hProj,
		const int YL, const int ZL, const int ViewN,
		const float PLC, const float ZLC,
		const float dYL, const float dZL,
		const float SO);


// The function routine should be called
// FPWD = reWeihAdFiltr(Proj, PLC, ZLC, dYL, dZL, SO);
void mexFunction(
		int nlhs, mxArray* plhs[],
		int nrhs, const mxArray* prhs[])
{
	// YL, ZL, ViewN should be calculated from Proj
	if(nlhs != 1 || nrhs != 6)
	{
		std::cerr<<"Wrong number of parameters\n";
		std::cerr<<"This function requires 6 input parameter and provides one return value\n";
		exit(-1);
	}
	const mwSize* siz = mxGetDimensions(prhs[0]);
	const int YL = siz[0];
	const int ZL = siz[1];
	const int ViewN = siz[2];
	// Generate a new array
	plhs[0] = mxCreateNumericArray(3,siz,mxSINGLE_CLASS,mxREAL);
	float* fpwd = (float*)mxGetPr(plhs[0]);
	const float* Proj = (float*)mxGetPr(prhs[0]);
	const float PLC = *((float*)mxGetData(prhs[1]));
	const float ZLC = *((float*)mxGetData(prhs[2]));
	const float dYL = *((float*)mxGetData(prhs[3]));
	const float dZL = *((float*)mxGetData(prhs[4]));
	const float SO = *((float*)mxGetData(prhs[5]));
    

	filtering(fpwd, Proj, YL, ZL, ViewN,
			PLC - 1.0f, ZLC - 1.0f, dYL, dZL, SO);
}














