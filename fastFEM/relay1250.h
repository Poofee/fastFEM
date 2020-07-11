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
