/*    
 *    selectSignificance_asm.cpp		
 *    DiffSplice
 *
 *    Copyright (C) 2012 University of Kentucky and
 *                       Yin Hu
 *
 *    Authors: Yin Hu
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "selectSignificance_asm.h"


/************************************************************************/
/* TEST GENE CLASS                                                      */
/************************************************************************/


ASM::ASM(int pathnum)
{
	ASM_id = "ASM0";
	testID = 0;
	rangeLow = 0;
	rangeHigh = 0;
	pathNum = pathnum;
	ASMcategory = 0;

	
	long groupLoop, sampleLoop;

	individualExpression = new double* [GROUP_NUM+1];
	error_ratio = new double* [GROUP_NUM+1];
	meanExpression = new double [GROUP_NUM+1];
	individualPathExpression = new double** [GROUP_NUM+1];
	individualPathProportion = new double** [GROUP_NUM+1];
	meanPathExpression = new double* [GROUP_NUM+1];
	meanPathProportion = new double* [GROUP_NUM+1];
	for (groupLoop = 1; groupLoop <= GROUP_NUM; ++groupLoop)
	{
		individualExpression[groupLoop] = new double [SAMPLE_CNT_PER_GROUP+1];
		error_ratio[groupLoop] = new double [SAMPLE_CNT_PER_GROUP+1];
		individualPathExpression[groupLoop] = new double* [SAMPLE_CNT_PER_GROUP+1];
		individualPathProportion[groupLoop] = new double* [SAMPLE_CNT_PER_GROUP+1];
	}

	statistic_d_individual = new double [num_blocks+1];
	statistic_s_individual = new double [num_blocks+1];
	statistic_d_expected_individual = new double [num_blocks+1];
	significant_individual = new bool [num_blocks+1];
	
	for (groupLoop = 1; groupLoop <= GROUP_NUM; ++groupLoop)
	{
		for (sampleLoop = 0; sampleLoop <= SAMPLE_CNT_PER_GROUP; ++sampleLoop)
		{
			individualExpression[groupLoop][sampleLoop] = 0.0;
			error_ratio[groupLoop][sampleLoop] = 0.0;
		}
		meanExpression[groupLoop] = 0.0;
	}	
	minExpression = 0.0;
	minPathExpression = 0.0;

	for (groupLoop = 1; groupLoop <= GROUP_NUM; ++groupLoop)
	{
		for (sampleLoop = 1; sampleLoop <= SAMPLE_CNT_PER_GROUP; ++sampleLoop)
		{
			individualPathExpression[groupLoop][sampleLoop] = new double [pathnum+1];
			individualPathProportion[groupLoop][sampleLoop] = new double [pathnum+1];
		}
		meanPathExpression[groupLoop] = new double [pathnum+1];
		meanPathProportion[groupLoop] = new double [pathnum+1];
	}	

	statistic_d = 0.;
	statistic_s = 0.;
	statistic_d_expected = 0.;

	for (groupLoop = 0; groupLoop <= num_blocks; ++groupLoop)
	{
		statistic_d_individual[groupLoop] = 0.0;
		statistic_s_individual[groupLoop] = 0.0;
		statistic_d_expected_individual[groupLoop] = 0.0;
		significant_individual[groupLoop] = false;
	}

	significant = false;
}


void ASM::calculate_meanPathExprProp()
{
	int pathLoopCnt, groupLoopCnt, sampleLoopCnt;
	double groupSum_expr, groupSum_prop;

	for (pathLoopCnt = 1; pathLoopCnt <= pathNum; ++pathLoopCnt)
	{
		for (groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
		{
			groupSum_expr = 0.0;
			groupSum_prop = 0.0;
			for (sampleLoopCnt = 1; sampleLoopCnt <= SAMPLE_CNT_PER_GROUP; ++sampleLoopCnt)
			{
				groupSum_expr += individualPathExpression[groupLoopCnt][sampleLoopCnt][pathLoopCnt];
				groupSum_prop += individualPathProportion[groupLoopCnt][sampleLoopCnt][pathLoopCnt];
			}
			meanPathExpression[groupLoopCnt][pathLoopCnt] = groupSum_expr / SAMPLE_CNT_PER_GROUP;
			meanPathProportion[groupLoopCnt][pathLoopCnt] = groupSum_prop / SAMPLE_CNT_PER_GROUP;
		}
	}

	//normalize the proportions 
	for (groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
	{
		groupSum_prop = 0.0;
		for (pathLoopCnt = 1; pathLoopCnt <= pathNum; ++pathLoopCnt)
		{
			groupSum_prop += meanPathProportion[groupLoopCnt][pathLoopCnt];
		}		
		if (groupSum_prop == 0.0)
		{
			groupSum_prop = 1.0;
		}
		for (pathLoopCnt = 1; pathLoopCnt <= pathNum; ++pathLoopCnt)
		{
			meanPathProportion[groupLoopCnt][pathLoopCnt] = meanPathProportion[groupLoopCnt][pathLoopCnt] / groupSum_prop;
		}		
	}

	return;
}

ASM::~ASM()
{
	long groupLoop, sampleLoop;

	for (groupLoop = 1; groupLoop <= GROUP_NUM; ++groupLoop)
	{
		for (sampleLoop = 1; sampleLoop <= SAMPLE_CNT_PER_GROUP; ++sampleLoop)
		{
			delete [] individualPathExpression[groupLoop][sampleLoop];
			delete [] individualPathProportion[groupLoop][sampleLoop];
		}
		delete [] meanPathExpression[groupLoop];
		delete [] meanPathProportion[groupLoop];
	}	

	for (groupLoop = 1; groupLoop <= GROUP_NUM; ++groupLoop)
	{
		delete [] individualExpression[groupLoop];
		delete [] error_ratio[groupLoop];
		delete [] individualPathExpression[groupLoop];
		delete [] individualPathProportion[groupLoop];
	}
	delete [] individualExpression;
	delete [] error_ratio;
	delete [] meanExpression;
	delete [] individualPathExpression;
	delete [] individualPathProportion;
	delete [] meanPathExpression;
	delete [] meanPathProportion;

	delete [] statistic_d_individual;
	delete [] statistic_s_individual;
	delete [] statistic_d_expected_individual;
	delete [] significant_individual;
}



/************************************************************************/
/* BASIC FUNCTIONS                                                      */
/************************************************************************/

// 
// long partition(long p,long r, double *sortArray)
// {
// 	long i, j;
// 	double x, tmp;
// 
// 	//randomized partition
// 	i = p + (double)rand()/ (RAND_MAX) * (r - p);
// 
// 	if (sortArray[r] != sortArray[i])
// 	{
// 		tmp = sortArray[r];
// 		sortArray[r] = sortArray[i];
// 		sortArray[i] = tmp;
// 	}
// 
// 	x = sortArray[r];
// 	i = p - 1;
// 
// 	for(j = p; j < r; j++)
// 	{
// 		if (sortArray[j] <= x)
// 		{
// 			i++;
// 
// 			if (sortArray[i] != sortArray[j])
// 			{
// 				tmp = sortArray[j];
// 				sortArray[j] = sortArray[i];
// 				sortArray[i] = tmp;
// 			}
// 		}
// 	}
// 
// 	if (sortArray[r] != sortArray[i+1])
// 	{
// 		tmp = sortArray[r];
// 		sortArray[r] = sortArray[i+1];
// 		sortArray[i+1]=tmp;
// 	}
// 
// 	return i+1;
// }
// 
// 
// void quicksort(double *sortArray, long length)
// {
// 	long top = 0, p, r, q;
// 
// 	stack[top++] = 1;
// 	stack[top++] = length;
// 
// 	while (top != 0)
// 	{
// 		r = stack[--top];
// 		p = stack[--top];
// 
// 		if(p>=r)
// 			continue;
// 
// 		q = partition(p, r, sortArray);
// 
// 		stack[top++] = p;
// 		stack[top++] = q - 1;
// 
// 		stack[top++] = q + 1;
// 		stack[top++] = r;
// 	}
// 
// 	return;
// }


void merge_gene_sort(long p, long q, long r)
{
	long n1, n2, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	for (i = 1; i <= n1; i++)
	{
		mergeSort_Larray_ASM[i] = sortKey_ASM[p + i - 1];
		mergeSort_LorderedList_ASM[i] = sortList_ASM[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray_ASM[j] = sortKey_ASM[q + j];
		mergeSort_RorderedList_ASM[j] = sortList_ASM[q + j];
	}

	mergeSort_Larray_ASM[n1 + 1] = MAX_NUMBER;
	mergeSort_Rarray_ASM[n2 + 1] = MAX_NUMBER;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray_ASM[i] <= mergeSort_Rarray_ASM[j])
		{
			sortKey_ASM[k] = mergeSort_Larray_ASM[i];
			sortList_ASM[k] = mergeSort_LorderedList_ASM[i];

			i++;
		} 
		else
		{
			sortKey_ASM[k] = mergeSort_Rarray_ASM[j];
			sortList_ASM[k] = mergeSort_RorderedList_ASM[j];

			j++;
		}
	}

	return;
}


void mergeSort_gene_sort(long sortList_size)
{
	//non-recursive merge sort for sorting junctions
	long m, n, i, r;
	m = 1;
	n = sortList_size;

	while (m <= n)
	{
		i = 1;
		while (i <= n - m)
		{
			r = (i + 2 * m - 1) < n ? (i + 2 * m - 1) : n;
			merge_gene_sort(i, i + m - 1, r);
			i = i + 2 * m;
		}

		m = m * 2;
	}

	return;
}



/************************************************************************/
/* INPUT DATA                                                           */
/************************************************************************/


void input_ASMs(string resultPath)
{
	string filename, cur_chrNM, cur_ASMid;
	int groupLoopCnt, sampleLoopCnt, pathLoopCnt, pathCnt, badEstimationCnt;
	bool flag_filter_ASMexpr, flag_filter_pathexpr;
	double groupSum, minPathExpression, MSE_estimation, minASMExpression;
	long ASMrangeLow, ASMrangeHigh;
	ASM *newASM;
	string info;

	ifstream testASM_file;
	filename = resultPath + "stat/asm.txt";
	testASM_file.open(filename.c_str());

	ofstream filteredASM_file;
	filename = resultPath + "stat/filteredASM.txt";
	filteredASM_file.open(filename.c_str());


	while (testASM_file >> cur_ASMid)
	{
		testASM_file >> cur_chrNM;
		testASM_file >> ASMrangeLow;
		testASM_file >> ASMrangeHigh;
		testASM_file >> pathCnt;

		newASM = new ASM(pathCnt);
		flag_filter_ASMexpr = false;
		minASMExpression = MAX_NUMBER;

		newASM->ASM_id = cur_ASMid;
		newASM->testID = sortList_ASM_Num;
		newASM->chrNM = cur_chrNM;
		newASM->rangeLow = ASMrangeLow;
		newASM->rangeHigh = ASMrangeHigh;

		testASM_file >> newASM->ASMcategory;

		//input ASM expression
		for (groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
		{
			groupSum = 0.0;
			badEstimationCnt = 0;
			for (sampleLoopCnt = 1; sampleLoopCnt <= SAMPLE_CNT_PER_GROUP; ++sampleLoopCnt)
			{
				testASM_file >> newASM->individualExpression[groupLoopCnt][sampleLoopCnt];
				testASM_file >> MSE_estimation;

				newASM->error_ratio[groupLoopCnt][sampleLoopCnt] = MSE_estimation;

				groupSum += newASM->individualExpression[groupLoopCnt][sampleLoopCnt];
				if (MSE_estimation > THRESH_BAD_ESTIMATION)
					++badEstimationCnt;
			}
			newASM->meanExpression[groupLoopCnt] = groupSum / SAMPLE_CNT_PER_GROUP;
			if (newASM->meanExpression[groupLoopCnt] < minASMExpression)
				minASMExpression = newASM->meanExpression[groupLoopCnt];

			if (newASM->meanExpression[groupLoopCnt] < THRESHOLD_MIN_ASM_COVERAGE)// || badEstimationCnt > SAMPLE_CNT_PER_GROUP/2)
			{
				flag_filter_ASMexpr = true;
			}
		}	
		newASM->minExpression = minASMExpression;
		getline(testASM_file, info);

		if (COUNT_INTRON_RETENTION == false && newASM->ASMcategory == 3)
		{
			//throw intron retention
			flag_filter_ASMexpr = true;
		}

		//filtering by mean expression
		if (flag_filter_ASMexpr == true)
		{
			filteredASM_file << newASM->chrNM << "\t" << newASM->rangeLow << "\t" << newASM->rangeHigh << "\t" << newASM->pathNum << "\t" << newASM->ASMcategory << "\t"
				<< newASM->meanExpression[1] << "\t" << newASM->meanExpression[2] << endl;
			for (pathLoopCnt = 1; pathLoopCnt <= pathCnt; ++pathLoopCnt)
			{
				getline(testASM_file, info);
				filteredASM_file << info << endl;
			}
			delete newASM;
		}
		else
		{
			//input ASM paths
			flag_filter_pathexpr = false;
			minPathExpression = MAX_NUMBER;

			for (pathLoopCnt = 1; pathLoopCnt <= pathCnt; ++pathLoopCnt)
			{
				testASM_file >> ASMrangeLow;
				testASM_file >> ASMrangeHigh;
				for (groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
				{
					groupSum = 0.0;
					for (sampleLoopCnt = 1; sampleLoopCnt <= SAMPLE_CNT_PER_GROUP; ++sampleLoopCnt)
					{
						testASM_file >> newASM->individualPathExpression[groupLoopCnt][sampleLoopCnt][pathLoopCnt];
						testASM_file >> newASM->individualPathProportion[groupLoopCnt][sampleLoopCnt][pathLoopCnt];
						groupSum += newASM->individualPathExpression[groupLoopCnt][sampleLoopCnt][pathLoopCnt];
					}

					if (groupSum / SAMPLE_CNT_PER_GROUP < minPathExpression)
						minPathExpression = groupSum / SAMPLE_CNT_PER_GROUP;

					if (groupSum / SAMPLE_CNT_PER_GROUP < THRESHOLD_MIN_ASMPATH_COVERAGE)
					{
						flag_filter_pathexpr = true;
					}
				}		
				getline(testASM_file, info);
			}

			if (flag_filter_pathexpr == true)
			{
				filteredASM_file << newASM->chrNM << "\t" << newASM->rangeLow << "\t" << newASM->rangeHigh << "\t" << newASM->pathNum << "\t" << newASM->ASMcategory << "\t"
					<< newASM->meanExpression[1] << "\t" << newASM->meanExpression[2] << endl;
				for (pathLoopCnt = 1; pathLoopCnt <= pathCnt; ++pathLoopCnt)
				{
					filteredASM_file << ASMrangeLow << "\t" << ASMrangeHigh;
					for (groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
						for (sampleLoopCnt = 1; sampleLoopCnt <= SAMPLE_CNT_PER_GROUP; ++sampleLoopCnt)
							filteredASM_file << "\t" << newASM->individualPathExpression[groupLoopCnt][sampleLoopCnt][pathLoopCnt] << "\t" << newASM->individualPathProportion[groupLoopCnt][sampleLoopCnt][pathLoopCnt];
					filteredASM_file << endl;
				}
				delete newASM;
			}
			else
			{
				newASM->calculate_meanPathExprProp();
				newASM->minPathExpression = minPathExpression;
				++sortList_ASM_Num;
				sortList_ASM.push_back(newASM);
			}
		}
	}

	testASM_file.close();
	filteredASM_file.close();

	return;
}

/************************************************************************/
/* CALCULATE STATISTICS                                                 */
/************************************************************************/


void calculate_permutationCnt(string resultFileSuffix, int individualIndex_start, int individualIndex_end)
{
	permType = perm_full;
	permutationCnt = 1;
	permutation_limit_ratio = 1;

	for (int individualCnt = individualIndex_start; individualCnt <= individualIndex_end; ++individualCnt)
	{
		if (TECHNICAL_REP_NUM[individualCnt] > 10)
		{
			//too many combinations, take the random mode
			permType = perm_random;
			permutationCnt = MAX_PERMUTATION_NUMBER;
			break;
		}

		permutationCnt *= choose_non_recursive(GROUP_NUM * TECHNICAL_REP_NUM[individualCnt], TECHNICAL_REP_NUM[individualCnt]); 
	}

	if (permutationCnt > MAX_PERMUTATION_NUMBER * 10)
	{
		permType = perm_reducedfull;
		permutation_limit_ratio = int(permutationCnt / MAX_PERMUTATION_NUMBER);
		permutationCnt = MAX_PERMUTATION_NUMBER;
	}
	else if (permutationCnt > MAX_PERMUTATION_NUMBER)
	{
		permType = perm_random;
		permutationCnt = MAX_PERMUTATION_NUMBER;
	}

	allpermutedvalueCnt = permutationCnt * sortList_ASM_Num;

	cout << resultFileSuffix << ":\t#asm = " << sortList_ASM_Num << "\t#permutation = " << permutationCnt << "\t#allpermutedvalue = " << allpermutedvalueCnt << endl;

	return;
}


void calculate_q0_q50(string resultPath, string resultFileSuffix)
{
	//find the 0% point and 50%point
	long index_q0, index_q50, loopCnt;
	string filename;	
	ifstream permutationValue_file;
	string info;
	filename = resultPath + "stat/asm/d_permutation_sorted_" + resultFileSuffix + ".txt";
	permutationValue_file.open(filename.c_str());

	index_q0 = long(ceil((allpermutedvalueCnt - 1) * 0.0) + 1);
	index_q50 = long(ceil((allpermutedvalueCnt - 1) * 0.50) + 1);

	for (loopCnt = 1; loopCnt < index_q0; ++loopCnt)
	{
		getline(permutationValue_file, info); 
	}
	permutationValue_file >> q0;

	for (++loopCnt; loopCnt < index_q50; ++loopCnt)
	{
		getline(permutationValue_file, info); 
	}
	permutationValue_file >> q50;

	permutationValue_file.close();
	
	return;
}



void calculate_pi0(int index)
{
	long geneLoopCnt, diCnt = 0;

	if (index == 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
		{
			if (sortList_ASM[geneLoopCnt]->statistic_d > q0 && sortList_ASM[geneLoopCnt]->statistic_d < q50)
			{
				++diCnt;
			}
		}
	} 
	else if (index > 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
		{
			if (sortList_ASM[geneLoopCnt]->statistic_d_individual[index] > q0 && sortList_ASM[geneLoopCnt]->statistic_d_individual[index] < q50)
			{
				++diCnt;
			}
		}
	}
	else
	{
		cout << "error: unrecognized index" << endl;
		exit(1);
	}

	pi0 = double(diCnt) / (0.5 * sortList_ASM_Num);
	pi0 = pi0 > 1? 1 : pi0;

	pi0 = pi0 < 0.5? 0.5 : pi0;

	return;
}


void calculate_statistics(string resultPath, string resultFileSuffix, int index)
{
//	getCounts(resultPath);
		
	calculate_q0_q50(resultPath, resultFileSuffix);
#ifndef UNIX
	q0 = -0.5422;
	q50 = 0.5422;
#endif
	calculate_pi0(index);

	//cout << resultFileSuffix << ":\tq0 = " << q0 << "\tq50 = " << q50 << "\tpi0 = " << pi0 << endl;

	return;
}


/************************************************************************/
/* CALCULATE FDR                                                        */
/************************************************************************/

long count_falseGeneCnt_onePermutation(double delta, long permutationIndex)
{
	long falseGeneCnt = 0, geneLoopCnt;

	for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
	{
		if (fabs(statistics_permuted_d[geneLoopCnt][permutationIndex] - statistics_expected_d[geneLoopCnt]) > delta)
		{
			//falsely called gene
			++falseGeneCnt;
		}
	}

	return falseGeneCnt;
}

void calculate_falseGeneCnt(double delta, double &falseGeneCnt_median, double &falseGeneCnt_90percentile, double &falseGeneCnt_mean, string resultPath)
{
	long permutationLoopCnt;

	for (permutationLoopCnt = 1; permutationLoopCnt <= permutationCnt; ++permutationLoopCnt)
	{
		falseGeneCntList[permutationLoopCnt] = count_falseGeneCnt_onePermutation(delta, permutationLoopCnt);
	}

	quicksort(falseGeneCntList, permutationCnt);

	falseGeneCnt_median = calculate_percentile(falseGeneCntList, permutationCnt, 0.5);
	falseGeneCnt_90percentile = calculate_percentile(falseGeneCntList, permutationCnt, 0.9);
	falseGeneCnt_mean = calculate_mean(falseGeneCntList, permutationCnt);

//	cout << falseGeneCnt_mean << "\t" << falseGeneCnt_median << '\t' <<falseGeneCnt_90percentile << "\t" << falseGeneCntList[permutationCnt] << '\t';
//	for (permutationLoopCnt = 1; permutationLoopCnt <= permutationCnt; ++permutationLoopCnt)
//	{
//		cout << falseGeneCntList[permutationLoopCnt] << "\t";
//	}

	return;
}


long selectSignificance(double delta, int index)
{
	long geneLoopCnt, significanceCnt = 0;
	
	if (index == 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
		{
			if (fabs(sortList_ASM[geneLoopCnt]->statistic_d - statistics_expected_d[geneLoopCnt]) > delta)
			{
				//significant gene
				sortList_ASM[geneLoopCnt]->significant = true;
				++significanceCnt;
			}
		}
	}
	else if (index > 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
		{
			if (fabs(sortList_ASM[geneLoopCnt]->statistic_d_individual[index] - statistics_expected_d[geneLoopCnt]) > delta)
			{
				//significant gene
				sortList_ASM[geneLoopCnt]->significant_individual[index] = true;
				++significanceCnt;
			}
		}
	}
	
	return significanceCnt;
}

void calculate_FDR(double delta, long significanceCnt, double &FDR_median, double &FDR_90percentile, double &FDR_mean, string resultPath)
{
	double falseGeneCnt_median = 0.0, falseGeneCnt_90percentile = 0.0, falseGeneCnt_mean = 0.0;

	calculate_falseGeneCnt(delta, falseGeneCnt_median, falseGeneCnt_90percentile, falseGeneCnt_mean, resultPath);
	
	falseGeneCnt_median = falseGeneCnt_median * pi0;
	falseGeneCnt_90percentile = falseGeneCnt_90percentile * pi0;
	falseGeneCnt_mean = falseGeneCnt_mean * pi0;

	if (significanceCnt != 0)
	{
		FDR_median = falseGeneCnt_median / significanceCnt;
		FDR_90percentile = falseGeneCnt_90percentile / significanceCnt;
		FDR_mean = falseGeneCnt_mean / significanceCnt;
	}
	else
	{
		FDR_median = 0.0;
		FDR_90percentile = 0.0;
		FDR_mean = 0.0;
	}

	return;
}

/************************************************************************/
/* MAIN                                                                 */
/************************************************************************/

void initialization()
{
	sortKey_ASM = new double [sortList_ASM_Num+1];
	mergeSort_Larray_ASM = new double [sortList_ASM_Num+1];
	mergeSort_Rarray_ASM = new double [sortList_ASM_Num+1];
	mergeSort_LorderedList_ASM = new ASM* [sortList_ASM_Num+1];
	mergeSort_RorderedList_ASM = new ASM* [sortList_ASM_Num+1];

	statistics_distribution = new double*** [sortList_ASM_Num+1];
	for (long geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
	{
		statistics_distribution[geneLoopCnt] = new double** [GROUP_NUM+1];
		for (long groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
			statistics_distribution[geneLoopCnt][groupLoopCnt] = new double* [SAMPLE_CNT_PER_GROUP+1];
	}

	statistics_dimensionList = new long [sortList_ASM_Num+1];
	statistics_minPathExprList = new double [sortList_ASM_Num+1];
	statistics_expected_d = new double [sortList_ASM_Num+1];
	statistics_d = new double [sortList_ASM_Num+1];
	statistics_s = new double [sortList_ASM_Num+1];

	statistics_permuted_d = new double* [sortList_ASM_Num+1];

	individualList = new int [num_blocks+1];

	return;
}

void clearall()
{
	delete [] sortKey_ASM;
	delete [] mergeSort_Larray_ASM;
	delete [] mergeSort_Rarray_ASM;
	delete [] mergeSort_LorderedList_ASM;
	delete [] mergeSort_RorderedList_ASM;

	for (long geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
		delete sortList_ASM[geneLoopCnt];
	sortList_ASM.clear();
	
	for (long geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
	{
		for (long groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
			delete [] statistics_distribution[geneLoopCnt][groupLoopCnt];
		delete [] statistics_distribution[geneLoopCnt];
	}
	delete [] statistics_distribution;

	delete [] statistics_dimensionList;
	delete [] statistics_minPathExprList;
	delete [] statistics_expected_d;
	delete [] statistics_d;
	delete [] statistics_s;

	delete [] statistics_permuted_d;

	delete [] individualList;

	if (TECHNICAL_REP_NUM != NULL)
		delete [] TECHNICAL_REP_NUM;
	if (TECHNICAL_REP_ARRAY_BASE != NULL)
		delete [] TECHNICAL_REP_ARRAY_BASE;

	return;
}


void test(string resultPath, string resultFileSuffix, int individual_index)
{
	string filename, comd;
	double delta = 0.0, FDR_median = 0.0, FDR_90percentile = 0.0, FDR_mean = 0.0;
	long significanceCnt = 0, geneLoopCnt, group;
	int groupLoopCnt, sampleLoopCnt, sampleNum, individualNum, distLoopCnt, curDimension;


	/************************************************************************/
	/* Differential analysis                                                */
	/************************************************************************/
	if (individual_index == 0)
	{
		calculate_permutationCnt(resultFileSuffix, 1, num_blocks);
		sampleNum = SAMPLE_CNT_PER_GROUP;
		individualNum = num_blocks;
	}
	else if (individual_index > 0)
	{
		calculate_permutationCnt(resultFileSuffix, individual_index, individual_index);
		sampleNum = TECHNICAL_REP_NUM[individual_index];
		individualNum = 1;
	}

	for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
	{
		statistics_permuted_d[geneLoopCnt] = new double [permutationCnt + 1];
	}

	for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
	{
		curDimension = sortList_ASM[geneLoopCnt]->pathNum;
		statistics_dimensionList[geneLoopCnt] = curDimension;
		statistics_minPathExprList[geneLoopCnt] = sortList_ASM[geneLoopCnt]->minPathExpression;
		for (groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
		{
			for (sampleLoopCnt = 1; sampleLoopCnt <= sampleNum; ++sampleLoopCnt)
			{
				statistics_distribution[geneLoopCnt][groupLoopCnt][sampleLoopCnt] = new double [curDimension+1];
				for (distLoopCnt = 1; distLoopCnt <= curDimension; ++distLoopCnt)
				{
					statistics_distribution[geneLoopCnt][groupLoopCnt][sampleLoopCnt][distLoopCnt] = sortList_ASM[geneLoopCnt]->individualPathProportion[groupLoopCnt][sampleLoopCnt + TECHNICAL_REP_ARRAY_BASE[individual_index]][distLoopCnt];
				}
			}
		}
	}

	diff_expression_analysis_two_class(statistics_distribution, statistics_dimensionList, statistics_minPathExprList, statistics_d, statistics_s, statistics_expected_d, statistics_permuted_d, sortList_ASM_Num, GROUP_NUM, individualNum, individualList, sampleNum, resultPath, resultFileSuffix);
	
	if (individual_index == 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
		{
			sortList_ASM[geneLoopCnt]->statistic_d = statistics_d[geneLoopCnt];
			sortList_ASM[geneLoopCnt]->statistic_s = statistics_s[geneLoopCnt];
		}
	}
	else if (individual_index > 0) 
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
		{
			sortList_ASM[geneLoopCnt]->statistic_d_individual[individual_index] = statistics_d[geneLoopCnt];
			sortList_ASM[geneLoopCnt]->statistic_s_individual[individual_index] = statistics_s[geneLoopCnt];
		}
	}

#ifdef UNIX
	comd = "sort -n +0 -1 " + resultPath + "stat/asm/d_permutation_" + resultFileSuffix + ".txt > " + resultPath + "stat/asm/d_permutation_sorted_" + resultFileSuffix + ".txt";
	system(comd.c_str());
#endif

	calculate_statistics(resultPath, resultFileSuffix, individual_index);


	/************************************************************************/
	/* Sort All Genes                                                       */
	/************************************************************************/
	if (individual_index == 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
		{
			sortKey_ASM[geneLoopCnt] = sortList_ASM[geneLoopCnt]->statistic_d;
		}
	}
	else if (individual_index > 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
		{
			sortKey_ASM[geneLoopCnt] = sortList_ASM[geneLoopCnt]->statistic_d_individual[individual_index];
		}
	}
	mergeSort_gene_sort(sortList_ASM_Num);

	if (individual_index == 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
		{
			sortList_ASM[geneLoopCnt]->statistic_d_expected = statistics_expected_d[geneLoopCnt];
		}
	}
	else if (individual_index > 0) 
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
		{
			sortList_ASM[geneLoopCnt]->statistic_d_expected_individual[individual_index] = statistics_expected_d[geneLoopCnt];
		}
	}
	
	/************************************************************************/
	/* Compute False Discovery Rate                                         */
	/************************************************************************/
	//stack = new long [permutationCnt+1];
	falseGeneCntList = new double [permutationCnt+1];

	filename = resultPath + "FDR_transcription_" + resultFileSuffix + ".txt";
	ofstream FDRfile(filename.c_str());
	FDRfile << "delta is the significance cutoff of the differential transcription statistics (calculated as |stat_d_expected-stat_d|)" << endl; 
	FDRfile << "shifting delta will call different sets of significant ASMs and will result in different FDRs" << endl;
	FDRfile << "you may choose your desired FDR according to this list" << endl << endl;
	for (delta = 0; delta <= 10; delta += 0.01)
	{
		significanceCnt = selectSignificance(delta, individual_index);
		calculate_FDR(delta, significanceCnt, FDR_median, FDR_90percentile, FDR_mean, resultPath);

		FDRfile << "delta = " << delta << ": " << significanceCnt << " ASMs picked, with FDR(median) = "
			<< FDR_median << " and FDR(mean) = " << FDR_mean << " and FDR(90percentile) = " << FDR_90percentile << endl;

		//cout << "delta = " << delta << ": " << significanceCnt << " ASMs picked, with FDR(median) = "
		//	<< FDR_median << " and FDR(mean) = " << FDR_mean << " and FDR(90percentile) = " << FDR_90percentile << endl;
		
		if (significance_cutoff < 1e-5 && FDR_median <= false_discovery_rate)
		{
			significance_cutoff = delta;
		}
		if (significanceCnt <= 0)
			break;
	}
	FDRfile.close();

	if (significance_cutoff < 1e-5)
	{
		significance_cutoff = delta;
	}


	/************************************************************************/
	/* Clean up allocated arrays                                            */
	/************************************************************************/
	for (long geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
	{
		delete [] statistics_permuted_d[geneLoopCnt];
	}
	//delete [] stack;
	delete [] falseGeneCntList;
	for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
		for (groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
			for (sampleLoopCnt = 1; sampleLoopCnt <= sampleNum; ++sampleLoopCnt)
				delete [] statistics_distribution[geneLoopCnt][groupLoopCnt][sampleLoopCnt];

	return;
}

void output_statistics(string resultPath)
{
	string filename;
	long geneLoopCnt;
	int individualLoopCnt;
	bool sameDirection, sameDirectionVSexpected;

	/************************************************************************/
	/* Sort All Genes                                                       */
	/************************************************************************/
	for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
	{
		sortKey_ASM[geneLoopCnt] = sortList_ASM[geneLoopCnt]->statistic_d;
	}
	mergeSort_gene_sort(sortList_ASM_Num);


	/************************************************************************/
	/* Diagnostics                                                          */
	/************************************************************************/
//	sprintf(filename, "%sstat/compare_d_asm.txt", resultPath);
//	ofstream d_statistics_outputfile(filename);
//	ASM *curASM;
// 	for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
// 	{
// 		curASM = sortList_ASM[geneLoopCnt];
// 		sameDirection = true;
// 		sameDirectionVSexpected = true;
// 
// 		d_statistics_outputfile << curASM->ASM_id << "\t" << curASM->chrNM << "\t" << curASM->rangeLow << "\t" <<curASM->rangeHigh << "\t" << curASM->ASMcategory << "\t" << curASM->pathNum << "\t"
// 			<< curASM->statistic_d << '\t' << curASM->statistic_d_expected << "\t" << fabs(curASM->statistic_d - curASM->statistic_d_expected) << "\t" << curASM->statistic_s << "\t"
// 			<< sqrt(calculate_JSD(curASM->meanPathProportion[1], curASM->meanPathProportion[2], curASM->pathNum)) << "\t";
// 
// 		for (individualLoopCnt = 1; individualLoopCnt <= num_blocks; ++individualLoopCnt)
// 		{
// 			d_statistics_outputfile << curASM->statistic_d_individual[individualLoopCnt] << "\t" << curASM->statistic_d_expected_individual[individualLoopCnt] << "\t"
// 				<< fabs(curASM->statistic_d_individual[individualLoopCnt] - curASM->statistic_d_expected_individual[individualLoopCnt]) << "\t" << curASM->statistic_s_individual[individualLoopCnt] << "\t";
// 	
// 			if (curASM->statistic_d_individual[individualLoopCnt] * curASM->statistic_d < 0)
// 			{
// 				sameDirection = false;
// 			}
// 			if ((curASM->statistic_d_individual[individualLoopCnt] - curASM->statistic_d_expected_individual[individualLoopCnt]) * (curASM->statistic_d - curASM->statistic_d_expected) < 0)
// 			{
// 				sameDirectionVSexpected = false;
// 			}
// 		}
// 
// 		d_statistics_outputfile << sameDirection << "\t" << sameDirectionVSexpected << "\t" << curASM->meanExpression[1] << "\t" << curASM->meanExpression[2] << endl;
// 	}
// 	d_statistics_outputfile.close();
	filename = resultPath + "differential_transcription.txt";
	ofstream d_statistics_outputfile(filename.c_str());
	ASM *curASM;
	long significantCnt = 0;
	double stat_difftrans, jsd;
	d_statistics_outputfile << "chromosome\tposition_start\tposition_end\tcategory\tstat_diff_trans(|stat_d_expected-stat_d|)\tsqrtJSD\tcoverage_group1\tcoverage_group2\tsignificant" <<endl;
	for (geneLoopCnt = 1; geneLoopCnt <= sortList_ASM_Num; ++geneLoopCnt)
	{
		curASM = sortList_ASM[geneLoopCnt];
		stat_difftrans = fabs(curASM->statistic_d - curASM->statistic_d_expected);
		jsd = sqrt(calculate_JSD(curASM->meanPathProportion[1], curASM->meanPathProportion[2], curASM->pathNum));

		d_statistics_outputfile << curASM->chrNM << "\t" << curASM->rangeLow << "\t" <<curASM->rangeHigh << "\t" << alterSpliceCategory[curASM->ASMcategory] << "\t"
			<< stat_difftrans << "\t" << jsd << "\t";
		d_statistics_outputfile << curASM->meanExpression[1] << "\t" << curASM->meanExpression[2] << "\t";
		if (stat_difftrans >= significance_cutoff && jsd >= thresh_JSD)
		{
			d_statistics_outputfile << "yes";
			++significantCnt;
		}
		else
			d_statistics_outputfile << "no";

#ifdef OUTPUT_PAIRWISE_JSD
		for (int samplecnt = 1; samplecnt <= SAMPLE_CNT_PER_GROUP; ++samplecnt)
		{
			if (curASM->individualExpression[1][samplecnt] > THRESHOLD_MIN_ASM_COVERAGE*2 && curASM->individualExpression[2][samplecnt] > THRESHOLD_MIN_ASM_COVERAGE*2
				&& curASM->error_ratio[1][samplecnt] < 0.1 && curASM->error_ratio[2][samplecnt] < 0.1)
				d_statistics_outputfile << "\t" << sqrt(calculate_JSD(curASM->individualPathProportion[1][samplecnt], curASM->individualPathProportion[2][samplecnt], curASM->pathNum));
			else
				d_statistics_outputfile << "\t0";
		}
#endif

		d_statistics_outputfile << endl;
	}
	d_statistics_outputfile.close();

	cout << "under FDR(median)<=" << false_discovery_rate << ", requiring sqrtJSD >=" << thresh_JSD << ", select " << significantCnt << " ASMs with significant change on transcription" << endl;
}


void input_config(string filename)
{
	ifstream config_file;
	config_file.open(filename.c_str());

	string parameter, info;

	if (config_file.is_open() == true)
	{
		while (config_file >> parameter)
		{
			if (parameter.compare("GROUP_NUM") == 0)
				config_file >> GROUP_NUM;
			else if (parameter.compare("INDIVIDUAL_NUM") == 0)
				config_file >> INDIVIDUAL_NUM;
			else if (parameter.compare("SAMPLE_CNT_PER_GROUP") == 0)
				config_file >> SAMPLE_CNT_PER_GROUP;
			else if (parameter.compare("COUNT_INTRON_RETENTION") == 0)
			{
				config_file >> info;
				if (info.compare("true") == 0)
					COUNT_INTRON_RETENTION = true;
				else
					COUNT_INTRON_RETENTION = false;
			}
			else if (parameter.compare("TOTAL_SAMPLE_CNT") == 0)
				config_file >> TOTAL_SAMPLE_CNT;
			else if (parameter.compare("Num_of_blocks") == 0)
				config_file >> num_blocks;
			else if (parameter.compare("Num_of_samples_per_block") == 0)
			{
				TECHNICAL_REP_NUM = new int [1+num_blocks];
				TECHNICAL_REP_NUM[0] = 0;
				for (int tmpCnt = 1; tmpCnt <= num_blocks; ++tmpCnt)
					config_file >> TECHNICAL_REP_NUM[tmpCnt];
			}
			else if (parameter.compare("Array_base_of_blocks") == 0)
			{
				TECHNICAL_REP_ARRAY_BASE = new int [1+num_blocks];
				TECHNICAL_REP_ARRAY_BASE[0] = 0;
				for (int tmpCnt = 1; tmpCnt <= num_blocks; ++tmpCnt)
					config_file >> TECHNICAL_REP_ARRAY_BASE[tmpCnt];
			}
			else if (parameter.compare("false_discovery_rate") == 0)
				config_file >> false_discovery_rate;
			else if (parameter.compare("thresh_JSD") == 0)
				config_file >> thresh_JSD;
			else if (parameter.compare("THRESHOLD_MIN_ASM_COVERAGE") == 0)
				config_file >> THRESHOLD_MIN_ASM_COVERAGE;

			getline(config_file, info);
		}
	} 
	else
	{
		cout << "Error: fail to open the config file for differential transcription test" << endl;
		exit(1);
	}

	config_file.close();

	if (GROUP_NUM <= 0 || INDIVIDUAL_NUM <= 0 || SAMPLE_CNT_PER_GROUP <= 0 || TOTAL_SAMPLE_CNT <= 0 || num_blocks <= 0)
	{
		cout << "Error: incomplete config file for differential transcription test" << endl;
		exit(1);
	}

	return;
}



int main(int argc, char *argv[])
{
	string resultPath, resultFileSuffix, config_filename;
	int individualLoopCnt;

#ifdef UNIX
	if (argc != 3)
	{
		cout << argv[0] << "\t<result_path>\t<config_file>" << endl;
		exit(1);
	}
	resultPath = argv[1];
	config_filename = argv[2];
#else
	resultPath = "result\\";
	config_filename = "config_testtrans";
#endif

	input_config(config_filename);

	/* initialize random seed: */
	srand ( time(NULL) );

	sortList_ASM.reserve(DEFAULT_TESTGENE_NUM);
	sortList_ASM.push_back(NULL);

	input_ASMs(resultPath);

	initialization();


	//test for all individuals
	resultFileSuffix = "two_group_comparison";
	for (individualLoopCnt = 1; individualLoopCnt <= num_blocks; ++individualLoopCnt)
	{
		individualList[individualLoopCnt] = individualLoopCnt;
	}
	test(resultPath, resultFileSuffix, 0);


	//test for each individual
// 	for (int individualLoopCnt = 1; individualLoopCnt <= num_blocks; ++individualLoopCnt)
// 	{
// 		sprintf(resultFileSuffix, "id_%d", individualLoopCnt);
// 		individualList[1] = individualLoopCnt;
// 		test(resultPath, resultFileSuffix, individualLoopCnt);
// 	}

	output_statistics(resultPath);

	clearall();

	return 0;
}


