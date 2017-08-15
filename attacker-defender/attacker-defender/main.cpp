// bug detection.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "computation.h"

void init(){
	if (!(ep = engOpen(NULL))) //�����Ƿ�����Matlab����ɹ���
	{
		cout << "Can't start Matlab engine!" << endl;
		exit(1);
	}
	engEvalString(ep, "clear;");
	engSetVisible(ep, true);

	outstuf.open("results.csv", ios::out);//��csv�ļ���׼��д��

	initPathAndPoint();  //��ʼ�������·����Ǳ�ڲ����
	initPureStrategies();  //��ʼ��˫���Ĵ����Լ���

	//�����ʼʱ˫����best
	outstuf << "����" << ",";//����
	//defender
	outstuf << "f1" << "," << "f2" << "," << "f3" << "," << "f4" << "," << "f5" << "," << "f6" << "," << "," << "Ud(a d)" << "," << "Ud(a D)" << ",";

	outstuf << "," << ",";//������

	//attacker
	outstuf << "a1" << "," << "a2" << "," << "a3" << "," << "a4" << "," << "a5" << "," << "a6" << "," << "," << "Ua(a d)" << "," << "Ua(A d)";
	outstuf << endl;
}

//���˫����best response
void outputBestResponse(){

	outstuf << totalRound << ",";//����
	//defender
	for (int i = 0; i < defenderBestStrategy.size(); i++)
		outstuf << defenderBestStrategy[i] << ",";
	outstuf << "," << defenderUtility_da << "," << defenderUtility_Da << ",";
	outstuf << "," << ",";//������
	//attacker
	for (int i = 0; i < attackerBestStrategy.size(); i++)
		outstuf << attackerBestStrategy[i] << ",";
	outstuf << "," << attackerUtility_ad << "," << attackerUtility_Ad;
	outstuf << endl;
}

//������ս��
void outputResult(){
	
	outstuf << endl;
	outstuf << endl;
	//���ʵ�����ü�����
	outstuf << "�����Ĵ����Ը���" << "," << "," << LEFT_NUMBER << endl;
	outstuf << "��������" << "," << "," << totalRound << endl;
	outstuf << "���뷽����" << "," << "," << attackerUtility_ad << endl;
	outstuf << endl;
	//�����ⷽ����
	outstuf << "����" << "," << "," << "����" << "," << "," << "��ⷽ��ϲ���" << endl;
	for (int i = 0; i < defenderPureStrategies.size(); i++){
		outstuf << i << "," << "," << defenderMixedStrategy[i] << "," << ",";
		for (int j = 0; j < TOTAL_PATH_NUMBER; j++){
			outstuf << defenderPureStrategies[i][j] << ",";
		}
		outstuf << endl;
	}
	outstuf << endl;

	//������뷽����
	outstuf << "����" << "," << "," << "����" << "," << "," << "���뷽��ϲ���" << endl;
	for (int i = 0; i < attackerPureStrategies.size(); i++){
		outstuf << i << "," << "," << attackerMixedStrategy[i] << "," << ",";
		for (int j = 0; j < TOTAL_PATH_NUMBER; j++){
			outstuf << attackerPureStrategies[i][j] << ",";
		}
		outstuf << endl;
	}
}


int main()
{
	init();//��ʼ��

	//��ʼѭ������double oracle�㷨
	while (true){
		cout << ++totalRound << endl;

		computeMixedStrategy();  //����˫���ڵ�ǰ�����Լ��µĻ�ϲ���
		
		computeAttackerBestResponse();  //����best response
		computeDefenderBestResponse();

		outputBestResponse();	//��˫�����µ�best responseд��csv�ļ�

		if (payoff_convergence())  //�ж���������˫�������õĴ����Ծ����ڴ����Լ����ڣ�������
			break;
		else{
			//deleteBadStrategy_average();  //ɾ��˫�������Լ����еĽϻ�����
			attackerPureStrategies.push_back(attackerBestStrategy);  //�����Ŵ����Է��봿���Լ���
			defenderPureStrategies.push_back(defenderBestStrategy);		
		}

		if (totalRound > 100)//��ֹ���ִ���ʱ��һֱѭ��
			cout<<"more than 100 rounds";
	}
	 
	outputResult();	//������ս��
	outstuf.close();	  //�ر��ļ�д����
	engClose(ep);	//�ر�matlab����

	return 0;
}


