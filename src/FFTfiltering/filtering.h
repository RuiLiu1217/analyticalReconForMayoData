/*
 * filtering.h
 *
 *  Created on: Aug 3, 2016
 *      Author: liurui
 */

#ifndef FILTERING_H_
#define FILTERING_H_


extern "C"
void filtering(float* hfpwd,
		const float* hProj,
		const int YL, const int ZL, const int ViewN,
		const float PLC, const float ZLC,
		const float dYL, const float dZL,
		const float SO);

#endif /* FILTERING_H_ */
