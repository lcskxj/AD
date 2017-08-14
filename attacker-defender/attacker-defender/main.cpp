// bug detection.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "computation.h"

int main()
{

	if (!(ep = engOpen(NULL))) //�����Ƿ�����Matlab����ɹ���
	{
		cout << "Can't start Matlab engine!" << endl;
		exit(1);
	}
	engEvalString(ep, "clear;");
	engSetVisible(ep, true);

	fstream outstuf;
	outstuf.open("results.csv", ios::out);

	initPathAndPoint();  //��ʼ�������·����Ǳ�ڲ����
	initPureStrategies();  //��ʼ��˫���Ĵ����Լ���

	while (true){
		cout << ++totalRound << endl;

		computeMixedStrategy();  //����˫���ڵ�ǰ�����Լ��µĻ�ϲ���
		
		computeAttackerBestResponse();  //�������Ŵ�����
		computeDefenderBestResponse();

		if (payoff_convergence())  //�ж���������˫�������õĴ����Ծ����ڴ����Լ����ڣ�������
			break;
		else{

			//deleteBadStrategy_average();  //ɾ��˫�������Լ����еĽϻ�����

			attackerPureStrategies.push_back(attackerBestStrategy);  //�����Ŵ����Է��봿���Լ���
			defenderPureStrategies.push_back(defenderBestStrategy);		
		}

		//defender
		for (int i = 0; i < defenderBestStrategy.size();i++)
			outstuf << defenderBestStrategy[i] << ",";
		outstuf << "," <<defenderUtility_da << ",";

		outstuf << "," << ",";

		//attacker
		for (int i = 0; i < attackerBestStrategy.size(); i++)
			outstuf << attackerBestStrategy[i] << ",";
		outstuf << "," << attackerUtility_ad << ",";


		outstuf << endl;
		cout << endl;
		
		if (totalRound > 200)
			break;
	}

	//�������ʱ��best
	//defender
	for (int i = 0; i < defenderBestStrategy.size(); i++)
		outstuf << defenderBestStrategy[i] << ",";
	outstuf << "," << defenderUtility_da << ",";

	outstuf << "," << ",";//������

	//attacker
	for (int i = 0; i < attackerBestStrategy.size(); i++)
		outstuf << attackerBestStrategy[i] << ",";
	outstuf << "," << attackerUtility_ad << ",";
	outstuf << endl;
	outstuf << endl;
	outstuf << endl;

	outstuf << "�����Ĵ����Ը���" << "," << "," << LEFT_NUMBER << endl;
	outstuf << "��������" << "," << "," << totalRound << endl;
	outstuf << "���뷽����" << "," << "," << attackerUtility_ad << endl;
	outstuf << endl;

	outstuf << "��ⷽ��ϲ���" << "," << "," << "," << endl;
	for (int i = 0; i < defenderPureStrategies.size(); i++){
		outstuf << defenderMixedStrategy[i] << "," << ",";
		for (int j = 0; j < TOTAL_PATH_NUMBER; j++){
			outstuf << defenderPureStrategies[i][j] << ",";
		}
		outstuf << endl;
	}

	outstuf << "���뷽��ϲ���" << "," << "," << "," << endl;
	for (int i = 0; i < attackerPureStrategies.size(); i++){
		outstuf << attackerMixedStrategy[i] << "," << ",";
		for (int j = 0; j < TOTAL_PATH_NUMBER; j++){
			outstuf << attackerPureStrategies[i][j] << ",";
		}
		outstuf << endl;
	}

	outstuf.close();//�رս�������
	engClose(ep);//�ر�matlab����

	return 0;
}


