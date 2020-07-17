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
