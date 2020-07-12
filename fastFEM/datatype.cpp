#include "datatype.h"
#include <math.h>

//double CNode::length()
//{
//	return sqrt(x*x + y*y);
//}
FEMface &FEMface::operator=(const FEMface &f)
{
    this->n1 = f.n1;
    this->n2 = f.n2;
    this->n3 = f.n3;
    return *this;
}

/*!
 \brief 三个点的编号大小都要比较大小。字典序算法。

 \param f
 \return bool FEMface::operator
*/
bool FEMface::operator <(const FEMface &f)
{
    if(this->n1 > f.n1)
        return false;

    if(this->n1 < f.n1){
        return true;
    }else{
        if(this->n2 > f.n2)
            return false;
        if(this->n2 < f.n2){
            return true;
        }else{
            if(this->n3 > f.n3)
                return false;
            if(this->n3 < f.n3){
                return true;
            }else{
                return false;
            }
        }
    }
}

bool FEMface::operator ==(const FEMface &f)
{
    return this->n1 == f.n1 && this->n2 == f.n2 && this->n3 == f.n3;
}

bool FEMface::operator !=(const FEMface &f)
{
    return !(this->n1 == f.n1 && this->n2 == f.n2 && this->n3 == f.n3);
}
