#include <omp.h>
#include <cmath>
#include <assert.h>
#include "ParallelRebinningCBCurve.h"
// Let's make the order of the projection data ZL, YL, ViewN
extern "C"
void ParallelRebinningCBCurve_CPU(
		float* outputProj,
		const float* Proj,
		const int YL, const int ZL, const int ViewN,
		const float YLC, const float dYA,
		const float DeltaTheta,
		const float PLC, const float DeltaT, const float DeltaFai,
		const float SO)
{
#pragma omp parallel for
	for(int i = 0; i < ViewN; i++)
	{
		float Theta = i * DeltaTheta; // The View for the parallel projection
		for(int j = 0; j != YL; ++j)
		{
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
			else if(FI > ViewN - 1)
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
			else if(UI > YL - 1)
			{
				IndexYU = YL - 1;
				IndexYB = YL - 1;
			}
			else
			{
				IndexYU = UI;
				IndexYB = UI - 1;
			}

			for(int k = 0; k != ZL; ++k)
			{
				outputProj[(i * YL + j) * ZL + k] =
						coeXB * coeYB * Proj[(IndexXB * YL + IndexYB) * ZL + k] +
						coeXU * coeYB * Proj[(IndexXU * YL + IndexYB) * ZL + k] +
						coeXB * coeYU * Proj[(IndexXB * YL + IndexYU) * ZL + k] +
						coeXU * coeYU * Proj[(IndexXU * YL + IndexYU) * ZL + k];
			}
		}
	}
}
