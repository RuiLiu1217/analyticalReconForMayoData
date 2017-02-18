#include "ParallelRebinningCBCurve.h"
#include <thrust/device_vector.h>
#include <thrust/copy.h>

__global__ void ParallelRebinningCBCurve_GPU_ker(
		float* outputProj,
		const float* Proj,
		const int YL, const int ZL, const int ViewN,
		const float YLC, 
        const float dYA, // Detector corresponding 
		const float DeltaTheta, // view step
        const float PLC,
		const float DeltaT, // ideal detector size
        const float DeltaFai,
        const float SO)
{
	int k = threadIdx.x + blockIdx.x * blockDim.x;
	int j = threadIdx.y + blockIdx.y * blockDim.y;
	int i = threadIdx.z + blockIdx.z * blockDim.z;
	if(i < ViewN && j < YL && k < ZL)
	{
        float Theta = i * DeltaTheta;
        float t = (j - PLC) * DeltaT;
		float Beta = asinf(t / SO);
		float Fai = Theta + Beta;
		float a = atanf(t / sqrtf(SO*SO-t*t));
		float FaiIndex = (Fai / DeltaFai);
		float UIndex = a / dYA + YLC;
		int FI = ceilf(FaiIndex);
		int UI = ceilf(UIndex);
		float coeXB = FI - FaiIndex;
		float coeXU = 1.0f - coeXB;
		float coeYB = UI - UIndex;
		float coeYU = 1.0f - coeYB;

		int IndexXU(0);
		int IndexXB(0);
		int IndexYU(0);
		int IndexYB(0);

		if(FI <= 0)
		{
			IndexXU = 0;
			IndexXB = 0;
		}
		else if(FI >= ViewN - 1)
		{
			IndexXU = ViewN - 1;
			IndexXB = ViewN - 1;
		}
		else
		{
			IndexXU = FI;
			IndexXB = FI - 1.0;
		}

		if(UI <= 0)
		{
			IndexYU = 0;
			IndexYB = 0;
		}
		else if(UI >= YL - 1)
		{
			IndexYU = YL - 1;
			IndexYB = YL - 1;
		}
		else
		{
			IndexYU = UI;
			IndexYB = UI - 1;
		}
		outputProj[(i * YL + j) * ZL + k] =
			coeXB * coeYB * Proj[(IndexXB * YL + IndexYB) * ZL + k] +
			coeXU * coeYB * Proj[(IndexXU * YL + IndexYB) * ZL + k] +
			coeXB * coeYU * Proj[(IndexXB * YL + IndexYU) * ZL + k] +
			coeXU * coeYU * Proj[(IndexXU * YL + IndexYU) * ZL + k];
	}
}


void ParallelRebinningCBCurve_GPU_fcn(
		float* outputProj,
		const float* Proj,
		const int YL, const int ZL, const int ViewN,
		const float YLC, 
        const float dYA, // Detector corresponding 
		const float DeltaTheta, // view step
        const float PLC,
		const float DeltaT, // ideal detector size
        const float DeltaFai,
        const float SO)
{
	dim3 blk(64,16,1);
	dim3 gid(
		(ZL + blk.x - 1) / blk.x,
		(YL + blk.y - 1) / blk.y,
		(ViewN + blk.z - 1) / blk.z);
	thrust::device_vector<float> dProj(Proj, Proj + YL * ZL * ViewN);
	thrust::device_vector<float> doutputProj(dProj.size(),0);
	ParallelRebinningCBCurve_GPU_ker<<<gid,blk>>>(
			thrust::raw_pointer_cast(&doutputProj[0]),
			thrust::raw_pointer_cast(&dProj[0]), YL, ZL, ViewN,
			YLC, dYA, DeltaTheta, PLC, DeltaT, DeltaFai, SO);

	thrust::copy(doutputProj.begin(),doutputProj.end(),outputProj);
	dProj.clear();
	doutputProj.clear();
}


extern "C"
void ParallelRebinningCBCurve_GPU(
		float* outputProj,
		const float* Proj,
		const int YL, const int ZL, const int ViewN,
		const float YLC, 
        const float dYA, // Detector corresponding 
		const float DeltaTheta, // view step
        const float PLC,
		const float DeltaT, // ideal detector size
        const float DeltaFai,
        const float SO)
{
    cudaDeviceReset();
    cudaSetDevice(0);
    ParallelRebinningCBCurve_GPU_fcn(outputProj, Proj,
		YL, ZL, ViewN, YLC, dYA, DeltaTheta, PLC, DeltaT, 
        DeltaFai, SO);
    cudaDeviceReset();
}





