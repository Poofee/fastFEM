// Copyright 2020 Poofee (https://github.com/Poofee)
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------
/*****************************************************************************
 *                                                                           *
 *                                                                           *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *  Authors: Poofee                                                          *
 *  Email:   poofee@qq.com                                                   *
 *  Address:                                                                 *
 *  Original Date: 2020-07-15                                                *
 *                                                                           *
 *****************************************************************************/
#include "datatype.h"
#include <stdio.h>
#include <stdlib.h>
#define PI 3.14159265358979323846
CMaterial::CMaterial() {
    miu = 1;//*4*PI*1e-7;
	BHpoints = 0;
    Hdata = nullptr;
    Bdata = nullptr;
	H_c = 0;
	Theta_m = 0;
    Jr = 0;
    I = 0;
    tau = 0;
    sigma = 0;
}

CMaterial::~CMaterial() {
	if (Bdata != NULL)
		free(Bdata);
	if (Hdata != NULL)
		free(Hdata);
}

/** 这个求解出来的应该是相对磁导率 **/
double CMaterial::getMiu(double B) {
	double slope, H;
	if (B < 1e-3){
		B += 1e-3;
	}
	if (BHpoints == 0) {
		return miu;
	} else {
		for (int i = 0; i < BHpoints - 2; i++) {
			if (B > Bdata[i] && B <= Bdata[i + 1]) {
				slope = (Hdata[i + 1] - Hdata[i]) / (Bdata[i + 1] - Bdata[i]);
				H = Hdata[i] + slope*(B - Bdata[i]);
				return B / H;
			}
		}
	}
    /** 如果曲线最后是平的，就会有无穷大的bug **/
    /** 但是，这时候我如果只知道B，该如何得到H？ **/
	int i = BHpoints - 3;
	slope = (Hdata[i + 1] - Hdata[i]) / (Bdata[i + 1] - Bdata[i]);
	H = Hdata[i] + slope*(B - Bdata[i]);
	return B / H ;
}

/*!
 \brief 这里的处理仅仅使用了线性处理

 \param B
 \return double
*/
double CMaterial::getdvdB(double B) {
	double slope, H, b;
    if (B < 1e-9){/** 防止B=0 **/
		return 0;
	} else if (B > 3) {
		//B = 2;
	}
	if (BHpoints == 0) {
		return 0;
	} else {
		for (int i = 0; i < BHpoints - 2; i++) {
			if (B >= Bdata[i] && B <= Bdata[i + 1]) {
				slope = (Hdata[i + 1] - Hdata[i]) / (Bdata[i + 1] - Bdata[i]);
				H = Hdata[i] + slope*(B - Bdata[i]);
				b = Hdata[i] - slope * Bdata[i];
				return -b/(B*B);
			}
		}
	}
	int  i = BHpoints - 3;
	slope = (Hdata[i + 1] - Hdata[i]) / (Bdata[i + 1] - Bdata[i]);
	H = Hdata[i] + slope*(B - Bdata[i]);
	b = Hdata[i] - slope * Bdata[i];

	return -b / (B*B);
}
