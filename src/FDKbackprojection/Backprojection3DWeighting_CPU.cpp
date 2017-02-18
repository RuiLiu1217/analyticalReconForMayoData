#include "Backprojection3DWeighting.h"
#include <cmath>
#include <complex>
#include <vector>
#include <thread>
#include <omp.h>


const static float pi = 3.141592654f;
const static float TWOPI = pi * 2.0f;

extern "C"
void Backprojection3DWeightingCPU(float* ctimage, const float* fp, const float* w, const int N_2pi,
		const int XN, const int YN, const int ZN,
		const float XNC, const float YNC, const float ZNC,
		const float dx, const float dy, const float dz,
		const int YL, const int ZL, const int PN,
		const float YLC, const float ZLC, const float dYL, const float dZL,
		const float h, const float BetaS, const float DeltaFai,
		const float N_pi, const float HSCoef, const float SO, const float SD, const float k1)
{   
    const float invPI = 1.0f / TWOPI;
    
    int nProcessors=omp_get_max_threads();                  
    //std::cout<<nProcessors<<std::endl;
    omp_set_num_threads(nProcessors);            
#pragma omp parallel for
    for(int yi=0;yi<YN;yi++)
    {
        const float y = (yi-YNC)*dy;
        for(int xi=0;xi<XN;xi++)
        {
            const float x  = (xi-XNC)*dx;

            for(int zi = 0; zi < ZN; ++zi)
            {
                ///compute the projection position for every grid on the image plane
                const float z = (zi-ZNC) * dz;
                const float Beta0 = 2.0f * pi * z / h;
                const int s0 = ceilf((Beta0-BetaS) / DeltaFai-0.5f);
                int s1 = s0 - ceilf(N_pi * HSCoef);
                int s2 = s0 + ceilf(N_pi * HSCoef);
                float res = 0; // used to accumulate the results
                if((s1 < PN) && (s2 > 0))
                {
                    s1 = (s1 < 0)?0:s1;
                    s2 = (s2 > PN - 1)?PN-1:s2;
                    for(int ProjInd = s1; ProjInd <= s2; ++ProjInd)
                    {
                            const float View = BetaS + ProjInd * DeltaFai;
                            const int d1   = N_pi-(s0-ProjInd); //d1 = ProjInd;
                            const int d2 = (ProjInd < s0)?(d1 + N_pi) : (d1 - N_pi);
                            const float UU = -x * cosf(View) - y * sinf(View);
                            const float Yr = -x * sinf(View) + y * cosf(View);
                            const float temp1 = sqrtf(SO * SO - Yr * Yr);
                            const float temp2 = (z-h * View * invPI);
                            const float Zr = temp2*(SO*SD)/(UU * temp1+SO * SO - Yr * Yr);
                            const float U1 = Yr/dYL+YLC;
                            const int U  = ceilf(U1);
                            const float V1 = Zr/dZL+ZLC;
                            const int V  = ceilf(V1);
                            const float Dey = U-U1;
                            const float Dez = V-V1;
                            //Linear interploate
                            if ((U>0)&&(U<YL)&&(V>0)&&(V<ZL))
                            {
                                const float touying =
                                        Dey *          Dez * fp[(ProjInd * ZL + (V-1)) * YL + (U-1)] +
                                        Dey * (1.0f - Dez) * fp[(ProjInd * ZL + (V)) * YL + (U-1)] +
                                        (1.0f - Dey) * Dez * fp[(ProjInd * ZL + (V-1)) * YL + U] +
                                        (1.0f - Dey) * (1.0f - Dez) * fp[(ProjInd * ZL + V) * YL + U];
                                const float weight1 = w[d1];
                                const float weight2 = w[d2];

                                const float Gama   = fabsf( temp2 / ( temp1 + UU));
                                const float Gama_C = (ProjInd < s0) ? fabsf((temp2 - 0.5f) / (temp1 - UU)) : fabsf((temp2 + 0.5f) / (temp1 - UU));
                                const float m1 = powf(Gama,  k1);    //m1     = std::real(std::pow(Gama,k1*h));
                                const float m2 = powf(Gama_C, k1);  //m2     = std::real(std::pow(Gama_C,k1*h));
                                const float tw = (weight2*m1+weight1*m2);
                                float weight = 0;
                                if(tw != 0)
                                {
                                    weight = (weight1*m2) / tw;
                                }                              

                                res += weight*touying*DeltaFai;
                            }// end if linear interpolation
                    }// end for projection
                }// end if range
                ctimage[(yi*XN+xi)*ZN+zi] = res;
            } // end zi
        }	// xi
    }// yi
}
