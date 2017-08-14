// bug detection.cpp : 定义控制台应用程序的入口点。
//

#include "computation.h"

int main()
{

	if (!(ep = engOpen(NULL))) //测试是否启动Matlab引擎成功。
	{
		cout << "Can't start Matlab engine!" << endl;
		exit(1);
	}
	engEvalString(ep, "clear;");
	engSetVisible(ep, true);

	fstream outstuf;
	outstuf.open("results.csv", ios::out);

	initPathAndPoint();  //初始化程序的路径和潜在插入点
	initPureStrategies();  //初始化双方的纯策略集合

	while (true){
		cout << ++totalRound << endl;

		computeMixedStrategy();  //计算双方在当前纯策略集下的混合策略
		
		computeAttackerBestResponse();  //计算最优纯策略
		computeDefenderBestResponse();

		if (payoff_convergence())  //判断收敛，若双方计算获得的纯策略均已在纯策略集合内，则收敛
			break;
		else{

			//deleteBadStrategy_average();  //删除双方纯策略集合中的较坏策略

			attackerPureStrategies.push_back(attackerBestStrategy);  //将最优纯策略放入纯策略集合
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

	//输出收敛时的best
	//defender
	for (int i = 0; i < defenderBestStrategy.size(); i++)
		outstuf << defenderBestStrategy[i] << ",";
	outstuf << "," << defenderUtility_da << ",";

	outstuf << "," << ",";//空两列

	//attacker
	for (int i = 0; i < attackerBestStrategy.size(); i++)
		outstuf << attackerBestStrategy[i] << ",";
	outstuf << "," << attackerUtility_ad << ",";
	outstuf << endl;
	outstuf << endl;
	outstuf << endl;

	outstuf << "保留的纯策略个数" << "," << "," << LEFT_NUMBER << endl;
	outstuf << "收敛轮数" << "," << "," << totalRound << endl;
	outstuf << "插入方收益" << "," << "," << attackerUtility_ad << endl;
	outstuf << endl;

	outstuf << "检测方混合策略" << "," << "," << "," << endl;
	for (int i = 0; i < defenderPureStrategies.size(); i++){
		outstuf << defenderMixedStrategy[i] << "," << ",";
		for (int j = 0; j < TOTAL_PATH_NUMBER; j++){
			outstuf << defenderPureStrategies[i][j] << ",";
		}
		outstuf << endl;
	}

	outstuf << "插入方混合策略" << "," << "," << "," << endl;
	for (int i = 0; i < attackerPureStrategies.size(); i++){
		outstuf << attackerMixedStrategy[i] << "," << ",";
		for (int j = 0; j < TOTAL_PATH_NUMBER; j++){
			outstuf << attackerPureStrategies[i][j] << ",";
		}
		outstuf << endl;
	}

	outstuf.close();//关闭结果输出流
	engClose(ep);//关闭matlab引擎

	return 0;
}


