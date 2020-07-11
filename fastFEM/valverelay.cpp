#include "valverelay.h"

ValveRelay::ValveRelay()
{

}

ValveRelay::~ValveRelay()
{

}

/*!
 \brief 初始化变量信息。

*/
void ValveRelay::init()
{
    /** 文件名，需不带扩展名 **/
    sprintf(fileName,"%s","valve/valve3dnew");
    /** 初始化时间步长 **/
    current_step = 0;
    min_time = 0;
    max_time = 0;

    /** 初始化变量 **/
    current_step = 0;
    Ddisplacements.resize(timesteps.size()+1);
    Ddisplacements.at(current_step) = 0;
    displacements.resize(timesteps.size()+1);
    displacements.at(current_step) = 0;
    velocities.resize(timesteps.size()+1);
    velocities.at(current_step) = 0;
    accelerations.resize(timesteps.size()+1);
    accelerations.at(current_step) = 0;
    PhiCoil.resize(timesteps.size()+1);
    PhiCoil.at(current_step) = 0;
    ICoil.resize(timesteps.size()+1);
    ICoil.at(current_step) = 0;
    magforces.resize(timesteps.size()+1);
    magforces.at(current_step) = 0;
    UCoil.resize(timesteps.size()+1);
    UCoil.at(current_step) = 420;

    MAX_NONLINEARSTEPS = 100;
    /** 设置材料参数 **/
    CMaterial* air_material = new CMaterial;
    materialList.emplace_back(air_material);

    CMaterial* pm_material = new CMaterial;
    materialList.emplace_back(pm_material);
    pm_material->H_c = 1e6;
    pm_material->miu = 1.12;

    CMaterial* dt4e_material = new CMaterial;
    materialList.emplace_back(dt4e_material);
    dt4e_material->sigma = 1/(2e-7);
    dt4e_material->BHpoints = 20;
    double* Hdata = (double*)malloc(20*sizeof(double));
    double* Bdata = (double*)malloc(20*sizeof(double));
    /** valve的BH曲线 **/


    CMaterial* coil_material = new CMaterial;
    materialList.emplace_back(coil_material);
    coil_material->tau = 60/(3e-3 * 17e-3);

    /** 设置映射关系 **/
    materialMap[1] = dt4e_material;/** 衔铁 **/
    materialMap[2] = coil_material;/** 线圈 **/
    materialMap[3] = dt4e_material;/** 铁芯 **/
    materialMap[4] = dt4e_material;/** 铁芯 **/
    materialMap[5] = dt4e_material;/** 铁芯 **/
    materialMap[6] = pm_material;/** 永磁 **/
    materialMap[7] = dt4e_material;/** 铁芯（导磁环） **/
    materialMap[8] = dt4e_material;/** 铁芯 **/
    materialMap[9] = air_material;/** 外部空气 **/
    materialMap[10] = air_material;/** 可压缩空气 **/
    /** 设置形变区域 **/
    tag_xiantie = 1;
    tag_air = 10;

    Precision = 1e-6;
}

void ValveRelay::calcMagForce(int index)
{

}

/*!
 \brief 使用牛顿非线性迭代法计算

*/
void ValveRelay::NRSolve()
{

}

void ValveRelay::run()
{

}

void ValveRelay::outputResults()
{

}
