extern "C"
void Backprojection3DWeightingCPU(float* ctimage, const float* fp, const float* w, const int N_2pi,
	const int XN, const int YN, const int ZN,
	const float XNC, const float YNC, const float ZNC,
	const float dx, const float dy, const float dz,
	const int YL, const int ZL, const int PN,
	const float YLC, const float ZLC, const float dYL, const float dZL,
	const float h, const float BetaS, const float DeltaFai,
	const float N_pi, const float HSCoef, const float SO, const float SD, const float k1);

extern "C"
void Backprojection3DWeightingGPU(float* hctimage, const float* hfp, const float* hw, const int N_2pi,
	const int XN, const int YN, const int ZN,
	const float XNC, const float YNC, const float ZNC,
	const float dx, const float dy, const float dz,
	const int YL, const int ZL, const int PN,
	const float YLC, const float ZLC, const float dYL, const float dZL,
	const float h, const float BetaS, const float DeltaFai,
	const float N_pi, const float HSCoef, const float SO, const float SD, const float k1);