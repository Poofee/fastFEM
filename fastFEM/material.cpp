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

//�����������Ӧ������Դŵ���
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
	//������������ƽ�ģ��ͻ���������bug
	//���ǣ���ʱ�������ֻ֪��B������εõ�H��
	int i = BHpoints - 3;
	slope = (Hdata[i + 1] - Hdata[i]) / (Bdata[i + 1] - Bdata[i]);
	H = Hdata[i] + slope*(B - Bdata[i]);
	return B / H ;
}
//����Ĵ������ʹ�������Դ���
double CMaterial::getdvdB(double B) {
	double slope, H, b;
	if (B < 1e-9){//��ֹB=0
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
