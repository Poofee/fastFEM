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
#ifndef RELAY1250_H
#define RELAY1250_H

#include "relay.h"

/*!
 \brief 1250模型的特殊性。在稳定状态，有永磁的保持力和反力。
 一般来说，在稳定位置，都是机械反力使衔铁稳定的，不需要先计算磁场。

*/
class Relay1250 : public Relay
{
public:
    Relay1250();
    ~Relay1250();

    void init();
    void calcMagForce(int index);
    void makeTriangle(int index);
    void AxisNRsolve();
    void AxisTLMsolve();
    void run();
    void outputResults();

protected:


};

#endif // RELAY1250_H
