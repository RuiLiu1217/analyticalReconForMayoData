extern "C"
void ParallelRebinningCBCurve_CPU(
		float* outputProj,
		const float* Proj,
		const int YL, const int ZL, const int ViewN,
		const float YLC, const float dYA,
		const float DeltaTheta,
		const float PLC, const float DeltaT, const float DeltaFai,
		const float SO);

extern "C"
void ParallelRebinningCBCurve_GPU(float* outputProj, 
    const float* Proj, const int YL, const int ZL, const int ViewN,
	const float YLC, const float dYA, const float DeltaTheta,
	const float PLC, const float DeltaT, const float DeltaFai,	const float SO);
// 
// extern "C"
// void ParallelRebinningCBCurve_CPU(
// 		float* outputProj,
// 		const float* Proj,
// 		const int YL, const int ZL, const int ViewN,
// 		const float YLC, 
//         const float dYA, // Detector corresponding 
// 		const float DeltaTheta, // view step
// 		const float DeltaT, // ideal detector size
//         const float rho0, // ideal S2O
//         const float rho1, // real S2O
//         const float deltaPhi, // offset of the phi
// 		const float S2D);
// 
// extern "C"
// void ParallelRebinningCBCurve_GPU(
// 		float* outputProj,
// 		const float* Proj,
// 		const int YL, const int ZL, const int ViewN,
// 		const float YLC, 
//         const float dYA, // Detector corresponding 
// 		const float DeltaTheta, // view step
// 		const float DeltaT, // ideal detector size
//         const float rho0, // ideal S2O
//         const float rho1, // real S2O
//         const float deltaPhi, // offset of the phi
// 		const float S2D);