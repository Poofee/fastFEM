#ifndef VALVERELAY_H
#define VALVERELAY_H

#include "relay.h"
/*!
 \brief 仿真flux的电磁阀静态模型

*/
class ValveRelay :public Relay
{
public:
    ValveRelay();
    ~ValveRelay();

    void init();
    void calcMagForce(int index);
    void run();
    void outputResults();

private:

};

#endif // VALVERELAY_H
