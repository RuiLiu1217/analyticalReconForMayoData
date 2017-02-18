#include <vector>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <cuda_runtime.h>
#include <iostream>
#include <fstream>
#include "Backprojection3DWeighting.h"

#ifndef nullptr
#define nullptr NULL
#endif

#define TWOPI 6.2831853f

__global__ void ParallekFDKHelical3DWeightingKer_texture_grouped(
		float* ctimage,
		cudaTextureObject_t* texObj, const int prjPerGroup,
		float2* __restrict cossinV,
		float* __restrict w,
		const int XN, const int YN, const int ZN, // cannot change
		float x, float y, float z, // used
		float dx, float dy, float dz,
		//const int YL, const int ZL,
		const int PN,
		float YLC, float ZLC, //cannot change
		float invdYL, float invdZL,
		float h_div_twopi,
		float BetaS, const float DeltaFai,
		float N_pi,
		int NHSCOEF, float SOSO, float SOSD, float k1)
{
	int zi = threadIdx.x + blockIdx.x * blockDim.x;
	int xi = threadIdx.y + blockIdx.y * blockDim.y;
	int yi = threadIdx.z + blockIdx.z * blockDim.z;
	if(zi < ZN && xi < XN && yi < YN)
	{

		x = (xi - x) * dx; // z,y,z cannot change
		y = (yi - y) * dy;
		z = (zi - z) * dz;

		// Calculate Beta0 according to the position
		int s0 = ceilf((z / h_div_twopi - BetaS) / DeltaFai - 0.5); //Center projection index ;
		int s1 = s0 - NHSCOEF;
		int s2 = s0 + NHSCOEF;
		dx = 0;

		s1 = (s1 < 0)?0:s1;
		s2 = (s2 > PN-1)?PN-1:s2;
		BetaS = BetaS + s1 * DeltaFai;
		for(; s1 <= s2; ++s1)
		{

			int d1 = N_pi - s0 + s1; // s0, s2 NO!
			int d2 = (s1 < s0)?(d1 + N_pi) : (d1 - N_pi);
			float weight1 = w[d1];
			float weight2 = w[d2];

			dy = (z - BetaS * h_div_twopi);

			float UU = -x * cossinV[s1].x - y * cossinV[s1].y; // x,y NO!
			float Yr = -x * cossinV[s1].y + y * cossinV[s1].x;
			dz = sqrtf(SOSO - Yr * Yr);

			float Zr = dy * (SOSD) / ((UU + dz) * dz);

			float U1 = Yr * invdYL + YLC;
			float V1 = Zr * invdZL + ZLC;
			d1 = s1 / prjPerGroup;
			d2 = s1 - prjPerGroup * d1;
			float touying = tex3D<float>(texObj[d1],U1+0.5, V1+0.5, d2+0.5);

			float m1 = powf(fabsf(dy / (dz + UU)), k1) * weight2;
			float m2 = powf((s1 < s0) ? fabsf((dy - 0.5) / (dz - UU)) : fabsf((dy + 0.5) / (dz - UU)), k1) * weight1;

			float weight = m2 / (m1 + m2);

			dx += weight * touying;
			BetaS += DeltaFai;
		}
		ctimage[(yi*XN+xi)*ZN+zi] = dx * DeltaFai;
	}

}

void createTextureObject(
	cudaTextureObject_t& texObj,
	cudaArray* d_prjArray,
	int Width, int Height, int Depth,
	float* sourceData,
	cudaMemcpyKind memcpyKind,
	cudaTextureAddressMode addressMode,
	cudaTextureFilterMode textureFilterMode,
	cudaTextureReadMode textureReadMode,
	bool isNormalized)
{
	cudaExtent prjSize;
	prjSize.width = Width;
	prjSize.height = Height;
	prjSize.depth = Depth;
	cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc<float>();

	cudaMalloc3DArray(&d_prjArray, &channelDesc, prjSize);
	cudaMemcpy3DParms copyParams = { 0 };
	copyParams.srcPtr = make_cudaPitchedPtr(
		(void*) sourceData, prjSize.width * sizeof(float),
		prjSize.width, prjSize.height);
	copyParams.dstArray = d_prjArray;
	copyParams.extent = prjSize;
	copyParams.kind = memcpyKind;
	cudaMemcpy3D(&copyParams);
	cudaResourceDesc resDesc;
	memset(&resDesc, 0, sizeof(resDesc));
	resDesc.resType = cudaResourceTypeArray;
	resDesc.res.array.array = d_prjArray;
	cudaTextureDesc texDesc;
	memset(&texDesc, 0, sizeof(texDesc));
	texDesc.addressMode[0] = addressMode;
	texDesc.addressMode[1] = addressMode;
	texDesc.addressMode[2] = addressMode;
	texDesc.filterMode = textureFilterMode;
	texDesc.readMode = textureReadMode;
	texDesc.normalizedCoords = isNormalized;

	cudaCreateTextureObject(&texObj, &resDesc, &texDesc, nullptr);
}



void destroyTextureObject(cudaTextureObject_t& texObj, cudaArray* d_array)
{
	cudaDestroyTextureObject(texObj);
	cudaFreeArray(d_array);
}


// use texture to accelerate the backprojection
void ParallekFDKHelical3DWeighting_GPU_texture_grouped(
		float* hctimage,
		const float* hfp,
		const float* hw, const int N_2pi,
		const int XN, const int YN, const int ZN,
		const float XNC, const float YNC, const float ZNC,
		const float dx, const float dy, const float dz,
		const int YL, const int ZL, const int PN,
		const float YLC, const float ZLC,
		const float dYL, const float dZL,
		const float h, const float BetaS, const float DeltaFai,
		const float N_pi, const float HSCoef, const float SO, const float SD, const float k1)
{
	dim3 blk(64,4,3); // It can be configured as you like!
	dim3 gid(
			(ZN + blk.x - 1) / blk.x,
			(XN + blk.y - 1) / blk.y,
			(YN + blk.z - 1) / blk.z);

	thrust::device_vector<float> ctimage(XN * YN * ZN, 0);
	thrust::device_vector<float> w(hw, hw + N_2pi);

	// Divide the projection into several groups
	int prjNumPerGroup = 8192; // We need to change this value according to the requirement, the maximum dim
	int groupNum = PN / prjNumPerGroup; // first N groups
	int lasGrpNum = PN - prjNumPerGroup * groupNum;
	if (lasGrpNum != 0)
	{
		++groupNum;
	}
	thrust::host_vector<thrust::device_vector<float> > fp(groupNum);
	for(int i = 0; i != groupNum - 1; ++i)
	{
		fp[i].resize(YL * ZL * prjNumPerGroup);
		thrust::fill(fp[i].begin(),fp[i].end(), 0);
		thrust::copy(hfp + i * YL * ZL * prjNumPerGroup,
				hfp + (i+1) * YL * ZL * prjNumPerGroup, fp[i].begin());
	}
	fp[groupNum - 1].resize(YL * ZL * prjNumPerGroup);
	thrust::fill(fp[groupNum - 1].begin(),fp[groupNum - 1].end(), 0);
	if(lasGrpNum != 0)
	{
		thrust::copy(hfp + (groupNum - 1) * YL * ZL * prjNumPerGroup,
				hfp + (groupNum - 1) * YL * ZL * prjNumPerGroup
				+ YL * ZL * lasGrpNum, fp[groupNum - 1].begin());
	}
	else
	{
		thrust::copy(hfp + (groupNum - 1) * YL * ZL * prjNumPerGroup,
				hfp + groupNum * YL * ZL * prjNumPerGroup,
				fp[groupNum - 1].begin());
	}

	std::vector<cudaTextureObject_t> texObj(groupNum);
	std::vector<cudaArray*> d_prjArray(groupNum);
	cudaTextureAddressMode addressMode = cudaAddressModeBorder;
	cudaTextureFilterMode textureFilterMode = cudaFilterModeLinear;
	cudaTextureReadMode textureReadMode = cudaReadModeElementType;
	for(int i = 0; i != groupNum; ++i)
	{
		createTextureObject(texObj[i], d_prjArray[i],
			YL, ZL, prjNumPerGroup,
			thrust::raw_pointer_cast(&fp[i][0]),
			cudaMemcpyDeviceToDevice,
			addressMode,
			textureFilterMode,
			textureReadMode,
			false);
	}
	thrust::device_vector<cudaTextureObject_t> dtexObj = texObj;
	thrust::host_vector<float2> cossinV(PN);
	for(int ProjInd = 0; ProjInd != PN; ++ProjInd)
	{
		const float View = BetaS + ProjInd * DeltaFai;
		cossinV[ProjInd] = make_float2(cosf(View), sinf(View));

	}
	thrust::device_vector<float2> dcossinV = cossinV;

	int NHSCOEF = ceilf(N_pi * HSCoef);
	float h_div_twopi = h / TWOPI;
	ParallekFDKHelical3DWeightingKer_texture_grouped<<<gid,blk>>>(
			thrust::raw_pointer_cast(&ctimage[0]),
			thrust::raw_pointer_cast(&dtexObj[0]), prjNumPerGroup,
			thrust::raw_pointer_cast(&dcossinV[0]),
			thrust::raw_pointer_cast(&w[0]), XN, YN, ZN, XNC, YNC, ZNC,
			dx, dy, dz,
			PN, YLC, ZLC, 1.0 / dYL, 1.0 / dZL, h_div_twopi, BetaS, DeltaFai, N_pi,
			NHSCOEF, SO * SO, SO * SD, k1);


	for(int i = 0; i != groupNum; ++i)
	{
		destroyTextureObject(texObj[i], d_prjArray[i]);
	}
	thrust::copy(ctimage.begin(),ctimage.end(),hctimage);

}

///////////////////////////////////////////////////////////////////////////////



extern "C"
void Backprojection3DWeightingGPU(
		float* hctimage, const float* hfp, const float* hw, const int N_2pi,
		const int XN, const int YN, const int ZN,
		const float XNC, const float YNC, const float ZNC,
		const float dx, const float dy, const float dz,
		const int YL, const int ZL, const int PN,
		const float YLC, const float ZLC, const float dYL, const float dZL,
		const float h, const float BetaS, const float DeltaFai,
		const float N_pi, const float HSCoef, const float SO, const float SD, const float k1)
{
  cudaDeviceReset();
  cudaSetDevice(0);

  ParallekFDKHelical3DWeighting_GPU_texture_grouped(hctimage, hfp, hw, N_2pi,
			XN, YN, ZN, XNC, YNC, ZNC, dx, dy, dz,
			YL, ZL, PN, YLC, ZLC, dYL, dZL, h, BetaS, DeltaFai,
			N_pi, HSCoef, SO, SD, k1);
  cudaDeviceReset();
}


