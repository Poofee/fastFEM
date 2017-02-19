#include "datatype.h"
#include <stdio.h>
#include <stdlib.h>

CMaterial::CMaterial() {
	miu = 1;
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

double CMaterial::GetMiu(double B) {
	double slope, H;
	if (BHpoints == 0) {
		return 1;
	} else {
		for (int i = 0; i < BHpoints - 1; i++) {
			if (B >= Bdata[i] && B <= Bdata[i + 1]) {
				slope = (Hdata[i + 1] - Hdata[i]) / (Bdata[i + 1] - Bdata[i]);
				H = Hdata[i] + slope*(B - Bdata[i]);
				return B / H;
			}
		}
	}
	int i = BHpoints - 2;
	slope = (Hdata[i + 1] - Hdata[i]) / (Bdata[i + 1] - Bdata[i]);
	H = Hdata[i] + slope*(B - Bdata[i]);
	return B / H;
}
