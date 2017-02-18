//////////////////////////////////
//  Imaging and Informatics Lab
//  Department of Electrical and Computer Engineering
//  University of Massachusetts Lowell
//  Parallel-cone backprojection
//  May 28, 2016
//////////////////////////////////

#include <mex.h>
#include <cmath>
#include <complex>
#include <vector>
#include <cstdlib>
#include "omp.h"
#include <thread>
#include <iostream>
#include <fstream>

#define TYPE double
const static TYPE pi = 3.14159265358979;
const static float TWOPI = pi * 2.0f;
#include "Backprojection3DWeighting_CPU.h"
#include "Backprojection3DWeighting_GPU.h"



extern "C"
void ParallelFDKHelical3DWeighting(
		double* ctimage, // Image to be reconstructed
		const double* fp,
		const double SO,
		const double DO,
		const int YL,
		const int ZL,
		const double DecWidth,
		const double DecHeigh,
		const double YLC,
		const double ZLC,
		const double h,
		const double BetaS,
		const int PN,
		const int N_2pi,
		const double ObjR,
		const int XN,
		const int YN,
		const int ZN,
		const double XNC,
		const double YNC,
		const double ZNC,
		const int delta,
		const double HSCoef,
		const double k1,
        const int useGPU)
{
        const double PI = 3.141592653589793;
        const double N_pi = N_2pi / 2.0;
        const double dx = 2.0 * ObjR / XN;
        const double dy = 2.0 * ObjR / YN;
        const double dz = dx;
        const double dYL = DecWidth / YL;
        const double dZL = DecHeigh / ZL;
        const double DeltaFai = 2.0 * PI / N_2pi;
        const double inv2PI = 1.0 / (2.0 * PI);
        const double SD = SO + DO;
        const double SO_square = SO * SO;
        std::vector<double> w(N_2pi,0);
        const int L = 2.0 * ceil(N_pi * HSCoef) + 1.0;
        const int Shift = N_pi - ceil(N_pi * HSCoef);
    
        ////Producing the weighting function

        for (int k=0;k<L;k++)
        {
            if (0 <= k && k<delta)
                w[k+Shift]= pow(cos((pi/2)*(delta-k-0.5)/delta),2);
            else if(L-delta<=k && k < L)
                w[k+Shift]= pow(cos((pi/2)*(k-(L-delta)+0.5)/delta),2);
            else
                w[k+Shift] = 1;
        }

        
        if(useGPU==1)
        {

            // use single floating data to implement the backprojection
             std::vector<float> fw(w.begin(),w.end());
             float* fctimage = new float[XN * YN * ZN];

             float* ffp = new float[YL * ZL * PN];
             std::copy(fp, fp + YL * ZL * PN, ffp);

             Backprojection3DWeightingGPU(fctimage,ffp, &fw[0],N_2pi,
                     XN, YN, ZN, XNC, YNC, ZNC, dx, dy, dz,
                      YL, ZL, PN, YLC, ZLC, dYL, dZL, h, BetaS, DeltaFai,
                      N_pi, HSCoef, SO, SD, k1);
             std::copy(fctimage, fctimage + XN * YN * ZN, ctimage);
             delete[] fctimage;
             delete[] ffp;
        }
        else
        {
            std::vector<float> fctimage(ctimage, ctimage + XN * YN * ZN);
            std::vector<float> ffp(fp, fp + YL * ZL * PN);
            std::vector<float> fw(w.begin(), w.end());

            
            Backprojection3DWeightingCPU(&fctimage[0], &ffp[0], &fw[0], N_2pi,
                    XN, YN, ZN, XNC, YNC, ZNC, dx, dy, dz,
                     YL, ZL, PN, YLC, ZLC, dYL, dZL, h, BetaS, DeltaFai,
                     N_pi, HSCoef, SO, SD, k1);
            std::copy(fctimage.begin(),fctimage.end(), ctimage);
        }
        
        //////end of the main code
}


static const TYPE *geom[4];
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    TYPE ObjR,SO,DO,YLC,ZLC,XNC,YNC,ZNC,DecWidth,DecHeigh,h,dYL,dZL,dx,dy,dz,DeltaFai,BetaS,k1,HSCoef;
    int  YL,ZL,XN,YN,ZN,N_2pi,delta,N_pi,PN,i,j,k;
    TYPE *ctimage,*w;
    const TYPE *fp;
    /* Check for proper number of arguments */
    if (nrhs != 4) {
        mexErrMsgTxt("Backward projection needs 4 inputs.");
    }
    if (nlhs != 0) {
        mexErrMsgTxt("Backward projection does not need any output.");
    }
    
    geom[0] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 0));
    SO            = TYPE(geom[0][0]);
    DO            = TYPE(geom[0][1]);
    YL            = int(TYPE(geom[0][2])+0.5);
    ZL            = int(TYPE(geom[0][3])+0.5);
    DecWidth      = TYPE(geom[0][4]);
    DecHeigh      = TYPE(geom[0][5]);
    YLC           = TYPE(geom[0][6]);
    ZLC           = TYPE(geom[0][7]);
    h             =     TYPE(geom[0][8])*DecHeigh;
    
    geom[1] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 1));
    BetaS         =     TYPE(geom[1][0]);
    PN            =     int(geom[1][2]);
    N_2pi         = int(geom[1][3]);
    
    geom[2] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 2));
    ObjR          =    TYPE(geom[2][0]);
    XN            = int(TYPE(geom[2][1])+0.5);
    YN            = int(TYPE(geom[2][2])+0.5);
    ZN            = int(TYPE(geom[2][3])+0.5);
    XNC           =     TYPE(geom[2][4]);
    YNC           =     TYPE(geom[2][5]);
    ZNC           =     TYPE(geom[2][6]);
    
    geom[3] = mxGetPr(mxGetFieldByNumber(prhs[0], 0, 3));
    delta         = int(TYPE(geom[3][0])+0.5);
    HSCoef        =     TYPE(geom[3][1]);// It is the value of Half_Scan/2*pi;
    k1            =     TYPE(geom[3][2])*TYPE(geom[0][8]);
    
    fp            =     mxGetPr(prhs[1]);
    ctimage       =     mxGetPr(prhs[2]);
    int useGPU = *((int*)mxGetPr(prhs[3]));
    
    
    N_pi = N_2pi/2.0;
    dx = 2.0*ObjR/XN;
    dy = 2.0*ObjR/YN;
    dz = dx;
    dYL = DecWidth/YL;
    dZL = DecHeigh/ZL;
    DeltaFai = 2.0*pi/N_2pi;
    
    ParallelFDKHelical3DWeighting(ctimage, fp, SO, DO, YL, ZL,
		DecWidth, DecHeigh, YLC - 1.0f, ZLC - 1.0f, h, BetaS, PN, N_2pi,
		ObjR, XN, YN, ZN, XNC - 1.0f, YNC - 1.0f, ZNC - 1.0f, delta, HSCoef, k1,useGPU);

}

