#include <iostream>
#include <vector>
#include <time.h>
#include<algorithm>
#include <math.h>
#include "engine.h"
#include <fstream>
using namespace std;

Engine *ep; //����Matlab����ָ��

#define TOTAL_PATH_NUMBER  6   //�����е�·������

#define CROSS_POINT_NUMBER 3  //�����н����ĸ���

#define ATTACKER_CAPACITY   TOTAL_PATH_NUMBER/3  //attacker����ܲ����bug��

#define REWARD  5   //��һ��bug��������attacker��õ�����

#define LEFT_NUMBER 5  //�����Կռ��б����ĸ���

double percentage[TOTAL_PATH_NUMBER];  //������нṹ�����ĵ����·���İٷֱ�

vector<vector<double>> pathRelation;  //�����и�·��֮��Ĺ�ϵ

vector<vector<double>> attackerPureStrategies, defenderPureStrategies;//˫���Ĵ����Լ���

vector<double> attackerMixedStrategy, defenderMixedStrategy;  //˫���Ļ�ϲ���

vector<double> attackerBestStrategy, defenderBestStrategy;  //˫����õ����Ŵ�����

double attackerAverageMixedLog[LEFT_NUMBER], defenderAverageMixedLog[LEFT_NUMBER];  //˫����ϲ��Ե�ϵ����ʷ��¼

double attackerUtility_ad;
double attackerUtility_Ad;

double defenderUtility_da;
double defenderUtility_Da;

int totalRound = 0;


//���ó����е�·����Ǳ�ڵ�
//������Ϊ5��·����ǰ����ΪX�ͣ������ΪV2�����м�����ΪY�ͣ��غϵ�ΪV9��V10��V11�������һ��Ϊֱ����,��18��Ǳ�ڲ���㡣
void initPathAndPoint(){
	//�ó���·���еĽ�������Ⱥͳ�������ʾ·��֮��Ĺ�ϵ�����ơ�
	//�þ����ʾ��ÿ�д���һ������㣬ÿ�д���һ��·�������Ϊ-1������Ϊ1�����һ�б�ʾ�����С�
	//3������㣬6��·��
	vector<double> point1 = { 1, 1, 0, 0, 0, 0 };
	pathRelation.push_back(point1);
	vector<double> point2 = { -1, 0, 1, 1, 0, 0};
	pathRelation.push_back(point2);
	vector<double> point3 = { 0, -1, 0, -1, 1, 1};
	pathRelation.push_back(point3);


	//������нṹ�����ĵ����·���İٷֱ�
	percentage[0] = 0.5;
	percentage[1] = 0.5;
	percentage[2] = 0.2;
	percentage[3] = 0.3;
	percentage[4] = 0.6;
	percentage[5] = 0.2;
}

//��ʼ��˫���Ĵ����Լ���
void initPureStrategies(){
	attackerPureStrategies.clear();  //��մ����Լ���
	defenderPureStrategies.clear();
	attackerBestStrategy.clear();
	defenderBestStrategy.clear();
	attackerBestStrategy.resize(TOTAL_PATH_NUMBER);

	//��ʼ��defender�Ĵ����ԣ�ƽ������
	defenderBestStrategy = { 0.5, 0.5, 0.25, 0.25, 0.375, 0.375 };
	defenderPureStrategies.push_back(defenderBestStrategy);


	//��ʼ��attacker�Ĵ�����,���ѡȡ�㣬���в���bug
	attackerBestStrategy[2] = 0.3;
	attackerBestStrategy[4] = 0.3;


	//srand((unsigned)time(0));
	//int tmp = 0;
	//for (int i = 0; i < ATTACKER_CAPACITY;i++){
	//	tmp = (rand() % TOTAL_PATH_NUMBER);
	//	attackerBestStrategy[tmp] = 0.3;
	//}

	attackerPureStrategies.push_back(attackerBestStrategy);
}

//���㴿�����µ�����UAD
double attackerUtilityForPureStrategy(vector<double> d, vector<double> a){
	double utility = 0;
	for (int i = 0; i < TOTAL_PATH_NUMBER; i++){
		utility += exp(-5*d[i]*a[i]) * (1- exp(-5*percentage[i]*a[i])) * REWARD;
	}
	return utility;
}

//����˫���ڵ�ǰ�����Լ��µĻ�ϲ���
void computeMixedStrategy(){

	//attackerĿ�꺯��
	vector<double> f_am;
	for (int i = 0; i < attackerPureStrategies.size(); i++){
		f_am.push_back(1);
	}
	mxArray *f_a = mxCreateDoubleMatrix(attackerPureStrategies.size(), 1, mxREAL);
	double* pa = mxGetPr(f_a);

	for (int i = 0; i < attackerPureStrategies.size(); i++){
		pa[i] = f_am[i];
	}
	engPutVariable(ep, "f_am", f_a);
	mxDestroyArray(f_a);


	//attacker���������,��Ϊ�����߲��ԣ���Ϊ�����߲���
	vector<vector<double>> attackerUtilityMatrix;
	//������������Զ�Ӧ��utility
	vector<double> tmp;
	double utility;
	for (int i = 0; i < defenderPureStrategies.size(); i++){
		attackerUtilityMatrix.push_back(tmp);
		for (int j = 0; j < attackerPureStrategies.size(); j++){
			//����utility
			utility = attackerUtilityForPureStrategy(defenderPureStrategies[i], attackerPureStrategies[j]);
			attackerUtilityMatrix[i].push_back(utility);
		}
	}

	mxArray *a_am = mxCreateDoubleMatrix(attackerUtilityMatrix.size(), attackerUtilityMatrix[0].size(), mxREAL);
	pa = mxGetPr(a_am);
	for (int i = 0; i < attackerUtilityMatrix[0].size(); i++)//����
	{
		for (int j = 0; j < attackerUtilityMatrix.size(); j++)//����
		{
			pa[i*attackerUtilityMatrix.size() + j] = attackerUtilityMatrix[j][i];
		}
	}
	engPutVariable(ep, "a_am", a_am);
	mxDestroyArray(a_am);




	//attacker��Լ�������ĳ����1��
	vector<double> b_am;
	for (int i = 0; i < defenderPureStrategies.size(); i++){
		b_am.push_back(-1);
	}
	mxArray *b_a = mxCreateDoubleMatrix(defenderPureStrategies.size(), 1, mxREAL);
	pa = mxGetPr(b_a);

	for (int i = 0; i < defenderPureStrategies.size(); i++){
		pa[i] = b_am[i];
	}
	engPutVariable(ep, "b_am", b_a);
	mxDestroyArray(b_a);

	//attacker��Լ�������0��
	vector<double> z_am;
	for (int i = 0; i < attackerPureStrategies.size(); i++){
		z_am.push_back(0);
	}
	mxArray *z_a = mxCreateDoubleMatrix(attackerPureStrategies.size(), 1, mxREAL);
	pa = mxGetPr(z_a);
	for (int i = 0; i < attackerPureStrategies.size(); i++){
		pa[i] = z_am[i];
	}
	engPutVariable(ep, "z_am", z_a);
	mxDestroyArray(z_a);

	//����attacker�Ļ�ϲ���
	engEvalString(ep, "[x_am,fval_am]=linprog(f_am,-a_am,b_am,[],[],z_am,f_am);x_am=x_am/fval_am; fval_am = 1/fval_am;");
	mxArray *output_matlab = mxCreateDoubleMatrix(attackerPureStrategies.size(), 1, mxREAL);
	output_matlab = engGetVariable(ep, "x_am");
	pa = mxGetPr(output_matlab);

	attackerMixedStrategy.clear();
	for (int i = 0; i < attackerPureStrategies.size(); i++){
		attackerMixedStrategy.push_back(pa[i]);
	}
	mxDestroyArray(output_matlab);


	mxArray *fval_matlab = mxCreateDoubleMatrix(1, 1, mxREAL);
	fval_matlab = engGetVariable(ep, "fval_am");
	pa = mxGetPr(fval_matlab);
	attackerUtility_ad = pa[0];

	mxDestroyArray(fval_matlab);
	
	

	//defender��Լ�������0��
	vector<double> z_dm;
	for (int i = 0; i < defenderPureStrategies.size(); i++){
		z_dm.push_back(0);
	}
	mxArray *z_d = mxCreateDoubleMatrix(defenderPureStrategies.size(), 1, mxREAL);
	pa = mxGetPr(z_d);
	for (int i = 0; i < defenderPureStrategies.size(); i++){
		pa[i] = z_dm[i];
	}
	engPutVariable(ep, "z_dm", z_d);
	mxDestroyArray(z_d);

	//����defender�Ļ�ϲ���
	engEvalString(ep, "[x_dm,fval]=linprog(b_am,a_am',f_am,[],[],z_dm,-b_am);fval = -fval; x_dm=x_dm/fval; fval = 1/fval; fval = -fval;");
	mxArray *output_matlab1 = mxCreateDoubleMatrix(defenderPureStrategies.size(), 1, mxREAL);
	output_matlab1 = engGetVariable(ep, "x_dm");
	pa = mxGetPr(output_matlab1);

	defenderMixedStrategy.clear();
	for (int i = 0; i < defenderPureStrategies.size(); i++){
		defenderMixedStrategy.push_back(pa[i]);
	}
	mxDestroyArray(output_matlab1);

	mxArray *fval_matlab1 = mxCreateDoubleMatrix(1, 1, mxREAL);
	fval_matlab1 = engGetVariable(ep, "fval");
	pa = mxGetPr(fval_matlab1);
	defenderUtility_da = pa[0];

	mxDestroyArray(fval_matlab1);	
}


//attacker��best response
void computeAttackerBestResponse(){
	//��������·������bug֮���õ���������
	double a;
	double best;
	double tmp;
	double bestUtility[TOTAL_PATH_NUMBER] = { 0 };
	for (int i = 0; i < TOTAL_PATH_NUMBER; i++){
		a = 0;
		best = 0;

		while (a < 1){

			tmp = 0;
			for (int j = 0; j < defenderPureStrategies.size(); j++){
				tmp += exp(-5 * a * defenderPureStrategies[j][i]) * (1 - exp(-5 * a * percentage[i])) * REWARD * defenderMixedStrategy[j];
			}

			if (tmp>best){
				best = tmp;
				attackerBestStrategy[i] = a;
			}
			a = a + 0.001;
		}
		bestUtility[i] = best;
	}

	//������·���������������ѡ��ǰn���
	double bestPath[TOTAL_PATH_NUMBER] = {0,1,2,3,4,5};
	for (int i = 0; i < ATTACKER_CAPACITY; i++){
		for (int j = 0; j < TOTAL_PATH_NUMBER - 1; j++){
			if (bestUtility[j]>bestUtility[j + 1]){  //�����������ŵ����
				//����utility
				tmp = bestUtility[j];
				bestUtility[j] = bestUtility[j + 1];
				bestUtility[j + 1] = tmp;

				//����·�����
				tmp = bestPath[j];
				bestPath[j] = bestPath[j + 1];
				bestPath[j + 1] = tmp;
				//����a
				tmp = attackerBestStrategy[j];
				attackerBestStrategy[j] = attackerBestStrategy[j+1];
				attackerBestStrategy[j + 1] = tmp;
			}
		}
	}

	//��ǰn����������·������bug������a=0
	for (int i = 0; i < ATTACKER_CAPACITY; i++){
		double pathNumber = bestPath[TOTAL_PATH_NUMBER - i - 1];
		attackerBestStrategy[pathNumber] = attackerBestStrategy[TOTAL_PATH_NUMBER - i - 1];
	}
	for (int i = ATTACKER_CAPACITY; i < TOTAL_PATH_NUMBER; i++){
		double pathNumber = bestPath[TOTAL_PATH_NUMBER - i - 1];
		attackerBestStrategy[pathNumber] = 0;
	}

	//�ڵ�ǰattacker��ȡbest response��defender��ȡ��ϲ�������£�����Ϊ��
	attackerUtility_Ad = bestUtility[0] + bestUtility[1];
}

//defender��best response
void computeDefenderBestResponse(){


	//0625 1200 A�ǹ����ߴ����Լ��ϣ�6 x n�� testFun
	mxArray *attackerPures = mxCreateDoubleMatrix(TOTAL_PATH_NUMBER, attackerPureStrategies.size(), mxREAL);
	double* pa = mxGetPr(attackerPures);

	for (int i = 0; i < attackerPureStrategies.size(); i++){
		for (int j = 0; j < attackerPureStrategies[0].size(); j++){
			pa[i*attackerPureStrategies[0].size() + j] = attackerPureStrategies[i][j];
		}
	}

	engPutVariable(ep, "A", attackerPures);
	mxDestroyArray(attackerPures);

	//M
	mxArray *M_mat = mxCreateDoubleMatrix(1, attackerPureStrategies.size(), mxREAL);
	pa = mxGetPr(M_mat);
	for (int i = 0; i < attackerPureStrategies.size(); i++){
		pa[i] = attackerMixedStrategy[i];
	}
	engPutVariable(ep, "M", M_mat);
	mxDestroyArray(M_mat);

	//P
	mxArray *P_mat = mxCreateDoubleMatrix(1, TOTAL_PATH_NUMBER, mxREAL);
	pa = mxGetPr(P_mat);
	for (int i = 0; i < TOTAL_PATH_NUMBER; i++){
		pa[i] = percentage[i];
	}
	engPutVariable(ep, "P", P_mat);
	mxDestroyArray(P_mat);

	//R
	double R_c[1] = { REWARD };
	mxArray *R_mat = mxCreateDoubleMatrix(1, 1, mxREAL);
	memcpy(mxGetPr(R_mat), R_c, sizeof(double));
	engPutVariable(ep, "R", R_mat);
	mxDestroyArray(R_mat);


	//Լ������ϵ��

	//Aeq ϵ������
	mxArray *Aeq_mat = mxCreateDoubleMatrix(pathRelation.size(), TOTAL_PATH_NUMBER, mxREAL);
	pa = mxGetPr(Aeq_mat);

	for (int i = 0; i < TOTAL_PATH_NUMBER; i++){
		for (int j = 0; j <pathRelation.size(); j++){
			pa[i*pathRelation.size() + j] = pathRelation[j][i];
		}
	}
	engPutVariable(ep, "Aeq", Aeq_mat);
	mxDestroyArray(Aeq_mat);

	//beq ��ʽ��
	double beq_c[CROSS_POINT_NUMBER] = {1,0,0};
	mxArray *beq_mat = mxCreateDoubleMatrix(CROSS_POINT_NUMBER, 1, mxREAL);
	pa = mxGetPr(beq_mat);
	for (int i = 0; i < CROSS_POINT_NUMBER; i++){
		pa[i] = beq_c[i];
	}
	engPutVariable(ep, "beq", beq_mat);
	mxDestroyArray(beq_mat);


	//Ŀ�꺯��
	engEvalString(ep, "addpath('D:/Program Files (x86)/MATLAB/R2016a/bin');");

	engEvalString(ep, "[row,col]=size(A);");
	engEvalString(ep, "x0 = zeros(row,1);");//x0 ��ʼֵ
	engEvalString(ep, "lb = zeros(row,1);");//lb  ȫ0��
	engEvalString(ep, "ub = ones(row,1);");//ub  ȫ1��

	engEvalString(ep, "[x_db, fval_db] = fmincon(@(x)testFun(x, A, M, P, R), x0, [], [], Aeq, beq, lb, ub);");
	

	double output[TOTAL_PATH_NUMBER];
	mxArray *output_matlab = mxCreateDoubleMatrix(1, TOTAL_PATH_NUMBER, mxREAL);
	output_matlab = engGetVariable(ep, "x_db");
	memcpy((void*)output, (void*)mxGetPr(output_matlab), sizeof(output));
	mxDestroyArray(output_matlab);

	for (int i = 0; i < TOTAL_PATH_NUMBER; i++){
		defenderBestStrategy[i]=output[i];
	}

	//�ڵ�ǰdefender��ȡbest response��attacker��ȡ��ϲ��Ե�����£�defender����Ϊ��
	mxArray *fval_matlab = mxCreateDoubleMatrix(1, 1, mxREAL);
	fval_matlab = engGetVariable(ep, "fval_db");
	pa = mxGetPr(fval_matlab);
	defenderUtility_Da = pa[0];
	mxDestroyArray(fval_matlab);
}

//ɾ��˫�������Լ����еĽϻ�����
//��Ϊ���֣�ÿ��ɾ�����Ĳ��ԣ�n��ƽ��ɾ�����Ĳ��ԣ�n�ּ��ۿ�����ɾ�����
void deleteBadStrategy_average(){

	if (totalRound < LEFT_NUMBER-1){  //20��֮ǰ���ۼ�
		for (int i = 0; i < totalRound; i++){
			attackerAverageMixedLog[i] = attackerAverageMixedLog[i] + attackerMixedStrategy[i];
			defenderAverageMixedLog[i] = defenderAverageMixedLog[i] + defenderMixedStrategy[i];
		}
	}
	else if (totalRound == LEFT_NUMBER-1){
		for (int i = 0; i < LEFT_NUMBER-1; i++){  //20��ʱ������ƽ������
			attackerAverageMixedLog[i] = attackerAverageMixedLog[i] / (totalRound - i);
			defenderAverageMixedLog[i] = defenderAverageMixedLog[i] / (totalRound - i);
		}
	}
	else{
		for (int i = 0; i < LEFT_NUMBER; i++){  //20��֮��ÿ����ƽ������ɾ��һ��������
			attackerAverageMixedLog[i] = (attackerAverageMixedLog[i] + attackerMixedStrategy[i]) / (totalRound - LEFT_NUMBER + 2);
			defenderAverageMixedLog[i] = (defenderAverageMixedLog[i] + defenderMixedStrategy[i]) / (totalRound - LEFT_NUMBER + 2);
		}
		
		//����ð��˼�룬����С���ŵ����ɾ������һ������
		double tmpLog;
		vector<double> tmpStratery;
		for (int i = 0; i < LEFT_NUMBER-1; i++){
			if (attackerAverageMixedLog[i] < attackerAverageMixedLog[i + 1]){  
				//����ƽ����ʷ��¼
				tmpLog = attackerAverageMixedLog[i];  
				attackerAverageMixedLog[i] = attackerAverageMixedLog[i+1];
				attackerAverageMixedLog[i + 1] = tmpLog;

				//���������Լ����е�λ��
				tmpStratery = attackerPureStrategies[i];
				attackerPureStrategies[i] = attackerPureStrategies[i+1];
				attackerPureStrategies[i + 1] = tmpStratery;
			}
		}

		//ɾ����С��
		attackerAverageMixedLog[LEFT_NUMBER-1] = 0;
		attackerPureStrategies.pop_back();



		//����ð��˼�룬����С���ŵ����ɾ������һ������
		double tmpLog2;
		vector<double> tmpStratery2;
		for (int i = 0; i < LEFT_NUMBER-1; i++){
			if (defenderAverageMixedLog[i] < defenderAverageMixedLog[i + 1]){
				//����ƽ����ʷ��¼
				tmpLog2 = defenderAverageMixedLog[i];
				defenderAverageMixedLog[i] = defenderAverageMixedLog[i + 1];
				defenderAverageMixedLog[i + 1] = tmpLog2;

				//���������Լ����е�λ��
				tmpStratery2 = defenderPureStrategies[i];
				defenderPureStrategies[i] = defenderPureStrategies[i + 1];
				defenderPureStrategies[i + 1] = tmpStratery2;
			}
		}

		//ɾ����С��
		defenderAverageMixedLog[LEFT_NUMBER-1] = 0;
		defenderPureStrategies.pop_back();

	}
}

//�ж���������˫�������õĴ����Ծ����ڴ����Լ����ڣ�������
bool strategy_convergence(){

	bool attackerFlag = false;
	bool defenderFlag = false;

	for (int i = 0; i < attackerPureStrategies.size(); i++){
		attackerFlag = true;
		for (int j = 0; j < attackerBestStrategy.size(); j++){
			if (abs(attackerBestStrategy[j] - attackerPureStrategies[i][j]) >0.001){  //����ǰ����������һ��ֵ��best�����е�ֵ��ͬ���򷵻�false,Ȼ��break�������Ա���һ��������
				attackerFlag = false;
				break;
			}
		}
		if (attackerFlag == true)  //����һ���������뵱ǰ���Ų��Էǳ��������break����ʾ���ҵ�
			break;
	}

	for (int i = 0; i < defenderPureStrategies.size(); i++){
		defenderFlag = true;
		for (int j = 0; j < defenderBestStrategy.size(); j++){
			if (abs(defenderBestStrategy[j] - defenderPureStrategies[i][j]) >0.001){
				defenderFlag = false;
				break;
			}
		}
		if (defenderFlag == true)  
		    break;
	}

	return (attackerFlag && defenderFlag);
}

//����payoff�ж��Ƿ�����
bool payoff_convergence(){
	//�ڶԷ����Թ̶�������£��Է���ȡ��ϲ��ԣ����ҷ������Ĵ����Ի�õ����治�Ȼ�ϲ��Ի�õ������ʱ����Ϊ����
	if ((attackerUtility_Ad <= attackerUtility_ad) && (defenderUtility_Da <= defenderUtility_da))
		return true;
	else
		return false;
}
