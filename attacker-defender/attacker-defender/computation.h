#include <iostream>
#include <vector>
#include <time.h>
#include<algorithm>
#include <math.h>
#include "engine.h"
#include <fstream>
using namespace std;

Engine *ep; //定义Matlab引擎指针

#define TOTAL_PATH_NUMBER  6   //程序中的路径数量

#define CROSS_POINT_NUMBER 3  //程序中交叉点的个数

#define ATTACKER_CAPACITY   TOTAL_PATH_NUMBER/3  //attacker最多能插入的bug数

#define REWARD  5   //若一个bug被触发，attacker获得的收益

#define LEFT_NUMBER 5  //纯策略空间中保留的个数

double percentage[TOTAL_PATH_NUMBER];  //程序固有结构决定的到达各路径的百分比

vector<vector<double>> pathRelation;  //程序中各路径之间的关系

vector<vector<double>> attackerPureStrategies, defenderPureStrategies;//双方的纯策略集合

vector<double> attackerMixedStrategy, defenderMixedStrategy;  //双方的混合策略

vector<double> attackerBestStrategy, defenderBestStrategy;  //双方获得的最优纯策略

double attackerAverageMixedLog[LEFT_NUMBER], defenderAverageMixedLog[LEFT_NUMBER];  //双方混合策略的系数历史记录

double attackerUtility_ad;
double attackerUtility_Ad;

double defenderUtility_da;
double defenderUtility_Da;

int totalRound = 0;


//设置程序中的路径和潜在点
//现设置为5条路径，前两条为X型（交叉点为V2），中间两条为Y型（重合点为V9，V10，V11），最后一条为直线型,共18个潜在插入点。
void initPathAndPoint(){
	//用程序路径中的交叉点的入度和出度来表示路径之间的关系和限制。
	//用矩阵表示，每行代表一个交叉点，每列代表一个路径。入度为-1，出度为1。最后一列表示增广列。
	//3个交叉点，6条路径
	vector<double> point1 = { 1, 1, 0, 0, 0, 0 };
	pathRelation.push_back(point1);
	vector<double> point2 = { -1, 0, 1, 1, 0, 0};
	pathRelation.push_back(point2);
	vector<double> point3 = { 0, -1, 0, -1, 1, 1};
	pathRelation.push_back(point3);


	//程序固有结构决定的到达各路径的百分比
	percentage[0] = 0.5;
	percentage[1] = 0.5;
	percentage[2] = 0.2;
	percentage[3] = 0.3;
	percentage[4] = 0.6;
	percentage[5] = 0.2;
}

//初始化双方的纯策略集合
void initPureStrategies(){
	attackerPureStrategies.clear();  //清空纯策略集合
	defenderPureStrategies.clear();
	attackerBestStrategy.clear();
	defenderBestStrategy.clear();
	attackerBestStrategy.resize(TOTAL_PATH_NUMBER);

	//初始化defender的纯策略，平均分配
	defenderBestStrategy = { 0.5, 0.5, 0.25, 0.25, 0.375, 0.375 };
	defenderPureStrategies.push_back(defenderBestStrategy);


	//初始化attacker的纯策略,随机选取点，进行插入bug
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

//计算纯策略下的收益UAD
double attackerUtilityForPureStrategy(vector<double> d, vector<double> a){
	double utility = 0;
	for (int i = 0; i < TOTAL_PATH_NUMBER; i++){
		utility += exp(-5*d[i]*a[i]) * (1- exp(-5*percentage[i]*a[i])) * REWARD;
	}
	return utility;
}

//计算双方在当前纯策略集下的混合策略
void computeMixedStrategy(){

	//attacker目标函数
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


	//attacker的收益矩阵,上为攻击者策略，左为防御者策略
	vector<vector<double>> attackerUtilityMatrix;
	//计算各个纯策略对应的utility
	vector<double> tmp;
	double utility;
	for (int i = 0; i < defenderPureStrategies.size(); i++){
		attackerUtilityMatrix.push_back(tmp);
		for (int j = 0; j < attackerPureStrategies.size(); j++){
			//计算utility
			utility = attackerUtilityForPureStrategy(defenderPureStrategies[i], attackerPureStrategies[j]);
			attackerUtilityMatrix[i].push_back(utility);
		}
	}

	mxArray *a_am = mxCreateDoubleMatrix(attackerUtilityMatrix.size(), attackerUtilityMatrix[0].size(), mxREAL);
	pa = mxGetPr(a_am);
	for (int i = 0; i < attackerUtilityMatrix[0].size(); i++)//列数
	{
		for (int j = 0; j < attackerUtilityMatrix.size(); j++)//行数
		{
			pa[i*attackerUtilityMatrix.size() + j] = attackerUtilityMatrix[j][i];
		}
	}
	engPutVariable(ep, "a_am", a_am);
	mxDestroyArray(a_am);




	//attacker的约束条件的常数项（1）
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

	//attacker的约束项（大于0）
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

	//计算attacker的混合策略
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
	
	

	//defender的约束项（大于0）
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

	//计算defender的混合策略
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


//attacker的best response
void computeAttackerBestResponse(){
	//计算所有路径插入bug之后获得的期望收益
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

	//对所有路径的收益进行排序，选出前n大的
	double bestPath[TOTAL_PATH_NUMBER] = {0,1,2,3,4,5};
	for (int i = 0; i < ATTACKER_CAPACITY; i++){
		for (int j = 0; j < TOTAL_PATH_NUMBER - 1; j++){
			if (bestUtility[j]>bestUtility[j + 1]){  //将收益最大的排到最后
				//交换utility
				tmp = bestUtility[j];
				bestUtility[j] = bestUtility[j + 1];
				bestUtility[j + 1] = tmp;

				//交换路径序号
				tmp = bestPath[j];
				bestPath[j] = bestPath[j + 1];
				bestPath[j + 1] = tmp;
				//交换a
				tmp = attackerBestStrategy[j];
				attackerBestStrategy[j] = attackerBestStrategy[j+1];
				attackerBestStrategy[j + 1] = tmp;
			}
		}
	}

	//将前n个收益最大的路径插入bug，其余a=0
	for (int i = 0; i < ATTACKER_CAPACITY; i++){
		double pathNumber = bestPath[TOTAL_PATH_NUMBER - i - 1];
		attackerBestStrategy[pathNumber] = attackerBestStrategy[TOTAL_PATH_NUMBER - i - 1];
	}
	for (int i = ATTACKER_CAPACITY; i < TOTAL_PATH_NUMBER; i++){
		double pathNumber = bestPath[TOTAL_PATH_NUMBER - i - 1];
		attackerBestStrategy[pathNumber] = 0;
	}

	//在当前attacker采取best response且defender采取混合策略情况下，收益为：
	attackerUtility_Ad = bestUtility[0] + bestUtility[1];
}

//defender的best response
void computeDefenderBestResponse(){


	//0625 1200 A是攻击者纯策略集合（6 x n） testFun
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


	//约束条件系数

	//Aeq 系数矩阵
	mxArray *Aeq_mat = mxCreateDoubleMatrix(pathRelation.size(), TOTAL_PATH_NUMBER, mxREAL);
	pa = mxGetPr(Aeq_mat);

	for (int i = 0; i < TOTAL_PATH_NUMBER; i++){
		for (int j = 0; j <pathRelation.size(); j++){
			pa[i*pathRelation.size() + j] = pathRelation[j][i];
		}
	}
	engPutVariable(ep, "Aeq", Aeq_mat);
	mxDestroyArray(Aeq_mat);

	//beq 等式列
	double beq_c[CROSS_POINT_NUMBER] = {1,0,0};
	mxArray *beq_mat = mxCreateDoubleMatrix(CROSS_POINT_NUMBER, 1, mxREAL);
	pa = mxGetPr(beq_mat);
	for (int i = 0; i < CROSS_POINT_NUMBER; i++){
		pa[i] = beq_c[i];
	}
	engPutVariable(ep, "beq", beq_mat);
	mxDestroyArray(beq_mat);


	//目标函数
	engEvalString(ep, "addpath('D:/Program Files (x86)/MATLAB/R2016a/bin');");

	engEvalString(ep, "[row,col]=size(A);");
	engEvalString(ep, "x0 = zeros(row,1);");//x0 初始值
	engEvalString(ep, "lb = zeros(row,1);");//lb  全0列
	engEvalString(ep, "ub = ones(row,1);");//ub  全1列

	engEvalString(ep, "[x_db, fval_db] = fmincon(@(x)testFun(x, A, M, P, R), x0, [], [], Aeq, beq, lb, ub);");
	

	double output[TOTAL_PATH_NUMBER];
	mxArray *output_matlab = mxCreateDoubleMatrix(1, TOTAL_PATH_NUMBER, mxREAL);
	output_matlab = engGetVariable(ep, "x_db");
	memcpy((void*)output, (void*)mxGetPr(output_matlab), sizeof(output));
	mxDestroyArray(output_matlab);

	for (int i = 0; i < TOTAL_PATH_NUMBER; i++){
		defenderBestStrategy[i]=output[i];
	}

	//在当前defender采取best response且attacker采取混合策略的情况下，defender收益为：
	mxArray *fval_matlab = mxCreateDoubleMatrix(1, 1, mxREAL);
	fval_matlab = engGetVariable(ep, "fval_db");
	pa = mxGetPr(fval_matlab);
	defenderUtility_Da = pa[0];
	mxDestroyArray(fval_matlab);
}

//删除双方纯策略集合中的较坏策略
//分为三种，每轮删除最差的策略，n轮平均删除最差的策略，n轮加折扣因子删除最差
void deleteBadStrategy_average(){

	if (totalRound < LEFT_NUMBER-1){  //20轮之前，累加
		for (int i = 0; i < totalRound; i++){
			attackerAverageMixedLog[i] = attackerAverageMixedLog[i] + attackerMixedStrategy[i];
			defenderAverageMixedLog[i] = defenderAverageMixedLog[i] + defenderMixedStrategy[i];
		}
	}
	else if (totalRound == LEFT_NUMBER-1){
		for (int i = 0; i < LEFT_NUMBER-1; i++){  //20轮时，进行平均处理
			attackerAverageMixedLog[i] = attackerAverageMixedLog[i] / (totalRound - i);
			defenderAverageMixedLog[i] = defenderAverageMixedLog[i] / (totalRound - i);
		}
	}
	else{
		for (int i = 0; i < LEFT_NUMBER; i++){  //20轮之后，每轮求平均，并删除一个纯策略
			attackerAverageMixedLog[i] = (attackerAverageMixedLog[i] + attackerMixedStrategy[i]) / (totalRound - LEFT_NUMBER + 2);
			defenderAverageMixedLog[i] = (defenderAverageMixedLog[i] + defenderMixedStrategy[i]) / (totalRound - LEFT_NUMBER + 2);
		}
		
		//利用冒泡思想，将最小的排到最后，删除最差的一个策略
		double tmpLog;
		vector<double> tmpStratery;
		for (int i = 0; i < LEFT_NUMBER-1; i++){
			if (attackerAverageMixedLog[i] < attackerAverageMixedLog[i + 1]){  
				//交换平均历史记录
				tmpLog = attackerAverageMixedLog[i];  
				attackerAverageMixedLog[i] = attackerAverageMixedLog[i+1];
				attackerAverageMixedLog[i + 1] = tmpLog;

				//交换纯策略集合中的位置
				tmpStratery = attackerPureStrategies[i];
				attackerPureStrategies[i] = attackerPureStrategies[i+1];
				attackerPureStrategies[i + 1] = tmpStratery;
			}
		}

		//删除最小的
		attackerAverageMixedLog[LEFT_NUMBER-1] = 0;
		attackerPureStrategies.pop_back();



		//利用冒泡思想，将最小的排到最后，删除最差的一个策略
		double tmpLog2;
		vector<double> tmpStratery2;
		for (int i = 0; i < LEFT_NUMBER-1; i++){
			if (defenderAverageMixedLog[i] < defenderAverageMixedLog[i + 1]){
				//交换平均历史记录
				tmpLog2 = defenderAverageMixedLog[i];
				defenderAverageMixedLog[i] = defenderAverageMixedLog[i + 1];
				defenderAverageMixedLog[i + 1] = tmpLog2;

				//交换纯策略集合中的位置
				tmpStratery2 = defenderPureStrategies[i];
				defenderPureStrategies[i] = defenderPureStrategies[i + 1];
				defenderPureStrategies[i + 1] = tmpStratery2;
			}
		}

		//删除最小的
		defenderAverageMixedLog[LEFT_NUMBER-1] = 0;
		defenderPureStrategies.pop_back();

	}
}

//判断收敛，若双方计算获得的纯策略均已在纯策略集合内，则收敛
bool strategy_convergence(){

	bool attackerFlag = false;
	bool defenderFlag = false;

	for (int i = 0; i < attackerPureStrategies.size(); i++){
		attackerFlag = true;
		for (int j = 0; j < attackerBestStrategy.size(); j++){
			if (abs(attackerBestStrategy[j] - attackerPureStrategies[i][j]) >0.001){  //若当前纯策略中有一个值与best策略中的值不同，则返回false,然后break，继续对比下一个纯策略
				attackerFlag = false;
				break;
			}
		}
		if (attackerFlag == true)  //若有一个纯策略与当前最优策略非常相近，则break，表示已找到
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

//根据payoff判断是否收敛
bool payoff_convergence(){
	//在对方策略固定的情况下（对方采取混合策略），我方产生的纯策略获得的收益不比混合策略获得的收益大时，即为收敛
	if ((attackerUtility_Ad <= attackerUtility_ad) && (defenderUtility_Da <= defenderUtility_da))
		return true;
	else
		return false;
}
