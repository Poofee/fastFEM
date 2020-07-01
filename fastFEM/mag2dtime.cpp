#include "mag2dtime.h"

Mag2DTime::Mag2DTime()
{

}

/*!
 \brief 用来初始化模型的各种数据

*/
void Mag2DTime::init()
{
    /** 初始化模型数据 **/
    MAX_NONLINEARSTEPS = 20;

    /** 初始化时间步 **/

}

/*!
 \brief 将分网在某个方向上产生形变

 \param dx
 \param dy
*/
void Mag2DTime::reMeshParallelShift(double dx, double dy)
{

}

/*!
 \brief 求解线性稀疏矩阵

*/
void Mag2DTime::solve()
{

}

/*!
 \brief 调用所有的步骤进行求解

*/
void Mag2DTime::run()
{

    /** 时间步循环 **/
    for(int i = 0; i < data.timesteps.size();i++){

        /** remesh **/

        /** 处理分网 **/

        /** 非线性迭代求解 **/
        for(int non_iter = 0;non_iter < MAX_NONLINEARSTEPS;non_iter++){

            /** 有限元装配 **/

        }


        /** 处理结果 **/

    }
    /** 输出信息 **/

}
