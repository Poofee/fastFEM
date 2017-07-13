#pragma once
#include "datatype.h"
#include "plot.h"
#include "slu_mt_ddefs.h"

class CFastFEMcore
{
public:
	int num_pts;//�ڵ���Ŀ
	int num_ele;//��Ԫ��Ŀ
	int numDomain;//����Ŀ
	CNode * pmeshnode;
    CElement * pmeshele;
    CElement4 * pmeshele4;
	char filename[256];//

	double  Precision;//���㾫��
	double  Relax;
	int		LengthUnits;//���ȵ�λ����������
	CMaterial* materialList;//���������

	Plot * thePlot;
	CFastFEMcore();
	~CFastFEMcore();	
	// load mesh
	int Load2DMeshCOMSOL(const char fn[]);//���������Ϣ����
    int LoadQ4MeshCOMSOL(const char fn[]);//�����Ľڵ㵥Ԫ
	bool StaticAxisymmetricTLM();//ʹ��TLM��̬��Գƴų��ļ��㺯��
	double CalcForce();//��������㺯��
	int openProject(QString proFile);//�򿪹����ļ�����
	int preCalculation();//Ԥ���㺯��
	int solve();//��⺯��
	void readProjectElement(QXmlStreamReader &reader);//��ȡproject�ڵ�
	void readDomainElement(QXmlStreamReader &reader, int i);//��ȡDomain�ڵ�
	void readBHElement(QXmlStreamReader &reader,int i);//��ȡBH�ڵ�
	int StaticAxisymmetricNR();//ʹ��NR��̬��Գƴų��ļ��㺯��
    void myTriSolve(int ncore, SuperMatrix *L, SuperMatrix *U,
                    int_t *perm_r, int_t *perm_c, SuperMatrix *B, int_t *info);
};

