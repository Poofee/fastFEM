#include "datatype.h"
#include <stdio.h>
#include <stdlib.h>
#define PI 3.14159265358979323846
CMaterial::CMaterial() {
	miu = 1*4*PI*1e-7;
	BHpoints = 0;
	Hdata = NULL;
	Bdata = NULL;
	H_c = 0;
	Theta_m = 0;
	Jr = 0;
}

CMaterial::~CMaterial() {
	if (Bdata != NULL)
		free(Bdata);
	if (Hdata != NULL)
		free(Hdata);
}

//这个求解出来的应该是相对磁导率
double CMaterial::getMiu(double B) {
	double slope, H;
	
	if (BHpoints == 0) {
		return miu;
	} else {
		for (int i = 0; i < BHpoints - 2; i++) {
			if (B >= Bdata[i] && B <= Bdata[i + 1]) {
				slope = (Hdata[i + 1] - Hdata[i]) / (Bdata[i + 1] - Bdata[i]);
				H = Hdata[i] + slope*(B - Bdata[i]);
				return B / H ;
			}
		}
	}
	//如果曲线最后是平的，就会有无穷大的bug
	//但是，这时候我如果只知道B，该如何得到H？
	int i = BHpoints - 3;
	slope = (Hdata[i + 1] - Hdata[i]) / (Bdata[i + 1] - Bdata[i]);
	H = Hdata[i] + slope*(B - Bdata[i]);
	return B / H ;
}
//这里的处理仅仅使用了线性处理
double CMaterial::getdvdB(double B) {
	double slope, H, b;
	if (B < 1e-9){//防止B=0
		return 0;
	}
	if (BHpoints == 0) {
		return 0;
	} else {
		for (int i = 0; i < BHpoints - 1; i++) {
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
