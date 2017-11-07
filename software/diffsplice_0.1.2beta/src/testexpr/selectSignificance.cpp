/*    
 *    selectSignificance.cpp		
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


#include "selectSignificance.h"


/************************************************************************/
/* TEST GENE CLASS                                                      */
/************************************************************************/

gene::gene()
{
	testID = 0;
	rangeLow = 0;
	rangeHigh = 0;

	long groupLoop, sampleLoop;

	individualExpression = new double* [GROUP_NUM+1];
	meanExpression = new double [GROUP_NUM+1];
	individualMedianExpression = new double* [GROUP_NUM+1];
	meanMedianExpression = new double [GROUP_NUM+1];
	noExpressedCnt = new int [GROUP_NUM+1];
	expressedCnt = new int [GROUP_NUM+1];
	for (groupLoop = 1; groupLoop <= GROUP_NUM; ++groupLoop)
	{
		individualExpression[groupLoop] = new double [SAMPLE_CNT_PER_GROUP+1];
		individualMedianExpression[groupLoop] = new double [SAMPLE_CNT_PER_GROUP+1];
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
			individualMedianExpression[groupLoop][sampleLoop] = 0.0;
		}
		meanExpression[groupLoop] = 0.0;
		meanMedianExpression[groupLoop] = 0.0;
		noExpressedCnt[groupLoop] = 0;
		expressedCnt[groupLoop] = 0;
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


gene::~gene()
{
	long groupLoop, sampleLoop;

	for (groupLoop = 1; groupLoop <= GROUP_NUM; ++groupLoop)
	{
		delete [] individualExpression[groupLoop];
		delete [] individualMedianExpression[groupLoop];
	}
	delete [] individualExpression;
	delete [] meanExpression;
	delete [] individualMedianExpression;
	delete [] meanMedianExpression;
	delete [] noExpressedCnt;
	delete [] expressedCnt;

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
		mergeSort_Larray_gene[i] = sortKey_gene[p + i - 1];
		mergeSort_LorderedList_gene[i] = sortList_gene[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray_gene[j] = sortKey_gene[q + j];
		mergeSort_RorderedList_gene[j] = sortList_gene[q + j];
	}

	mergeSort_Larray_gene[n1 + 1] = MAX_NUMBER;
	mergeSort_Rarray_gene[n2 + 1] = MAX_NUMBER;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray_gene[i] <= mergeSort_Rarray_gene[j])
		{
			sortKey_gene[k] = mergeSort_Larray_gene[i];
			sortList_gene[k] = mergeSort_LorderedList_gene[i];

			i++;
		} 
		else
		{
			sortKey_gene[k] = mergeSort_Rarray_gene[j];
			sortList_gene[k] = mergeSort_RorderedList_gene[j];

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


double calculate_foldchange(double expr_1, double expr_2)
{
	if (expr_1 < 1e-10)
	{
		return MAX_FOLDCHANGE;
	}
	return expr_2 / expr_1;
}

/************************************************************************/
/* INPUT DATA                                                           */
/************************************************************************/


void input_testGenes(string resultPath)
{
	string filename, cur_chrNM;
	int groupLoopCnt, sampleLoopCnt;
	bool flag_filter;
	double groupSum;
	gene *newTestGene;
	string info;

	ifstream testGene_file;
	filename = resultPath + "stat/expression.txt";
	testGene_file.open(filename.c_str());

	ofstream filteredGene_file;
	filename = resultPath + "stat/filteredGene.txt";
	filteredGene_file.open(filename.c_str());


	while (testGene_file >> cur_chrNM)
	{
		newTestGene = new gene;
		flag_filter = true;

		newTestGene->testID = sortList_gene_Num;
		newTestGene->chrNM = cur_chrNM;
		testGene_file >> newTestGene->rangeLow;
		testGene_file >> newTestGene->rangeHigh;

		for (groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
		{
			groupSum = 0.0;
			for (sampleLoopCnt = 1; sampleLoopCnt <= SAMPLE_CNT_PER_GROUP; ++sampleLoopCnt)
			{
				testGene_file >> newTestGene->individualExpression[groupLoopCnt][sampleLoopCnt];

				//normalize expression across samples
				newTestGene->individualExpression[groupLoopCnt][sampleLoopCnt] *= NORMALIZATION_FACTOR[SAMPLE_CNT_PER_GROUP*(groupLoopCnt-1) + sampleLoopCnt];

				if (newTestGene->individualExpression[groupLoopCnt][sampleLoopCnt] < 3)
					newTestGene->noExpressedCnt[groupLoopCnt] += 1;
				if (newTestGene->individualExpression[groupLoopCnt][sampleLoopCnt] > 20)
					newTestGene->expressedCnt[groupLoopCnt] += 1;
				groupSum += newTestGene->individualExpression[groupLoopCnt][sampleLoopCnt];
			}
			newTestGene->meanExpression[groupLoopCnt] = groupSum / SAMPLE_CNT_PER_GROUP;

			if (newTestGene->meanExpression[groupLoopCnt] >= THRESHOLD_MIN_EXPR_COVERAGE)
			{
				flag_filter = false;
			}
		}
		
		for (groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
		{
			groupSum = 0.0;
			for (sampleLoopCnt = 1; sampleLoopCnt <= SAMPLE_CNT_PER_GROUP; ++sampleLoopCnt)
			{
				testGene_file >> newTestGene->individualMedianExpression[groupLoopCnt][sampleLoopCnt];

				//normalize expression across samples
				newTestGene->individualMedianExpression[groupLoopCnt][sampleLoopCnt] *= NORMALIZATION_FACTOR[SAMPLE_CNT_PER_GROUP*(groupLoopCnt-1) + sampleLoopCnt];

				groupSum += newTestGene->individualMedianExpression[groupLoopCnt][sampleLoopCnt];
			}
			newTestGene->meanMedianExpression[groupLoopCnt] = groupSum / SAMPLE_CNT_PER_GROUP;
		}

		getline(testGene_file, info);

		//filtering by mean expression
		if (flag_filter == true)
		{
			filteredGene_file << newTestGene->chrNM << "\t" << newTestGene->rangeLow << "\t" << newTestGene->rangeHigh << "\t" << newTestGene->meanExpression[1] << "\t"
				<< newTestGene->meanExpression[2] << endl;
			delete newTestGene;
		}
		else
		{
			++sortList_gene_Num;
			sortList_gene.push_back(newTestGene);
		}
	}

	testGene_file.close();
	filteredGene_file.close();

	return;
}

/************************************************************************/
/* CALCULATE STATISTICS                                                 */
/************************************************************************/
// 
// void getCounts(char *resultPath)
// {
// 	long cur_permutedCnt = 0, cur_geneCnt = 0;
// 	char filename[1000];
// 	string info;
// 
// 	//get total count of the permuted d values
// 	ifstream permutationCnt_file;
// #ifdef UNIX
// 	sprintf(filename, "%sstat/d_permutation_cnt.txt", resultPath);
// #else
// 	sprintf(filename, "%sstat\\d_permutation_cnt.txt", resultPath);
// #endif
// 	permutationCnt_file.open(filename);
// 	while (permutationCnt_file >> cur_permutedCnt)
// 	{
// 		permutationCnt_file >> cur_geneCnt;
// 		permutationCnt_file >> permutationCnt;
// 		getline(permutationCnt_file, info);
// 		allpermutedvalueCnt += cur_permutedCnt;
// 		allgeneCnt += cur_geneCnt;
// 	}
// 	permutationCnt_file.close();
// 	
// 	cout << allgeneCnt << "\t" << permutationCnt << "\t";
// 	return;
// }

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

	allpermutedvalueCnt = permutationCnt * sortList_gene_Num;

	cout << resultFileSuffix << ":\t#gene = " << sortList_gene_Num << "\t#permutation = " << permutationCnt << "\t#allpermutedvalue = " << allpermutedvalueCnt << endl;

	return;
}


void calculate_q25_q75(string resultPath, string resultFileSuffix)
{
	//find the 25% point and 75%point
	long index_q25, index_q75, loopCnt;
	string filename;	
	ifstream permutationValue_file;
	string info;
	filename = resultPath + "stat/expr/d_permutation_sorted_" + resultFileSuffix + ".txt";
	permutationValue_file.open(filename.c_str());

	index_q25 = long(ceil((allpermutedvalueCnt - 1) * 0.25) + 1);
	index_q75 = long(ceil((allpermutedvalueCnt - 1) * 0.75) + 1);

	for (loopCnt = 1; loopCnt < index_q25; ++loopCnt)
	{
		getline(permutationValue_file, info); 
	}
	permutationValue_file >> q25;

	for (++loopCnt; loopCnt < index_q75; ++loopCnt)
	{
		getline(permutationValue_file, info); 
	}
	permutationValue_file >> q75;

	permutationValue_file.close();
	
	return;
}



void calculate_pi0(int index)
{
	long geneLoopCnt, diCnt = 0;

	if (index == 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
		{
			if (sortList_gene[geneLoopCnt]->statistic_d > q25 && sortList_gene[geneLoopCnt]->statistic_d < q75)
			{
				++diCnt;
			}
		}
	} 
	else if (index > 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
		{
			if (sortList_gene[geneLoopCnt]->statistic_d_individual[index] > q25 && sortList_gene[geneLoopCnt]->statistic_d_individual[index] < q75)
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

	pi0 = double(diCnt) / (0.5 * sortList_gene_Num);
	pi0 = pi0 > 1? 1 : pi0;

	pi0 = pi0 < 0.5? 0.5 : pi0;

	return;
}


void calculate_statistics(string resultPath, string resultFileSuffix, int index)
{
//	getCounts(resultPath);
		
	calculate_q25_q75(resultPath, resultFileSuffix);
#ifndef UNIX
	q25 = -0.5422;
	q75 = 0.5422;
#endif
	calculate_pi0(index);

	//cout << resultFileSuffix << ":\tq25 = " << q25 << "\tq75 = " << q75 << "\tpi0 = " << pi0 << endl;

	return;
}


/************************************************************************/
/* CALCULATE FDR                                                        */
/************************************************************************/

long count_falseGeneCnt_onePermutation(double delta, long permutationIndex)
{
	long falseGeneCnt = 0, geneLoopCnt;

	for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
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
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
		{
			if (fabs(sortList_gene[geneLoopCnt]->statistic_d - statistics_expected_d[geneLoopCnt]) > delta)
			{
				//significant gene
				sortList_gene[geneLoopCnt]->significant = true;
				++significanceCnt;
			}
		}
	}
	else if (index > 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
		{
			if (fabs(sortList_gene[geneLoopCnt]->statistic_d_individual[index] - statistics_expected_d[geneLoopCnt]) > delta)
			{
				//significant gene
				sortList_gene[geneLoopCnt]->significant_individual[index] = true;
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



void test(string resultPath, string resultFileSuffix, int individual_index)
{
	string filename, comd;
	double delta = 0.0, FDR_median = 0.0, FDR_90percentile = 0.0, FDR_mean = 0.0;
	long significanceCnt = 0, geneLoopCnt, group;
	int groupLoopCnt, sampleLoopCnt, sampleNum, individualNum;


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

	for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
	{
		statistics_permuted_d[geneLoopCnt] = new double [permutationCnt + 1];
	}

	for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
	{
		for (groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
		{
			for (sampleLoopCnt = 1; sampleLoopCnt <= sampleNum; ++sampleLoopCnt)
			{
				//changed on Aug26, use median expression to calculate differential expression level
				//statistics_expression[geneLoopCnt][groupLoopCnt][sampleLoopCnt] = sortList_gene[geneLoopCnt]->individualExpression[groupLoopCnt][sampleLoopCnt + TECHNICAL_REP_ARRAY_BASE[individual_index]];
				statistics_expression[geneLoopCnt][groupLoopCnt][sampleLoopCnt] = sortList_gene[geneLoopCnt]->individualMedianExpression[groupLoopCnt][sampleLoopCnt + TECHNICAL_REP_ARRAY_BASE[individual_index]];
			}
		}
	}

	diff_expression_analysis_two_class(statistics_expression, statistics_d, statistics_s, statistics_expected_d, statistics_permuted_d, sortList_gene_Num, GROUP_NUM, individualNum, individualList, sampleNum, resultPath, resultFileSuffix);
	
	if (individual_index == 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
		{
			sortList_gene[geneLoopCnt]->statistic_d = statistics_d[geneLoopCnt];
			sortList_gene[geneLoopCnt]->statistic_s = statistics_s[geneLoopCnt];
		}
	}
	else if (individual_index > 0) 
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
		{
			sortList_gene[geneLoopCnt]->statistic_d_individual[individual_index] = statistics_d[geneLoopCnt];
			sortList_gene[geneLoopCnt]->statistic_s_individual[individual_index] = statistics_s[geneLoopCnt];
		}
	}

#ifdef UNIX
	comd = "sort -n +0 -1 " + resultPath + "stat/expr/d_permutation_" + resultFileSuffix + ".txt > " + resultPath + "stat/expr/d_permutation_sorted_" + resultFileSuffix + ".txt";
	system(comd.c_str());
#endif

	calculate_statistics(resultPath, resultFileSuffix, individual_index);


	/************************************************************************/
	/* Sort All Genes                                                       */
	/************************************************************************/
	if (individual_index == 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
		{
			sortKey_gene[geneLoopCnt] = sortList_gene[geneLoopCnt]->statistic_d;
		}
	}
	else if (individual_index > 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
		{
			sortKey_gene[geneLoopCnt] = sortList_gene[geneLoopCnt]->statistic_d_individual[individual_index];
		}
	}
	mergeSort_gene_sort(sortList_gene_Num);

	if (individual_index == 0)
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
		{
			sortList_gene[geneLoopCnt]->statistic_d_expected = statistics_expected_d[geneLoopCnt];
		}
	}
	else if (individual_index > 0) 
	{
		for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
		{
			sortList_gene[geneLoopCnt]->statistic_d_expected_individual[individual_index] = statistics_expected_d[geneLoopCnt];
		}
	}
	
	/************************************************************************/
	/* Compute False Discovery Rate                                         */
	/************************************************************************/
	//stack = new long [permutationCnt+1];
	falseGeneCntList = new double [permutationCnt+1];

	filename = resultPath + "FDR_expression_" + resultFileSuffix + ".txt";
	ofstream FDRfile(filename.c_str());
	FDRfile << "delta is the significance cutoff of the differential expression statistics (calculated as |stat_d_expected-stat_d|)" << endl; 
	FDRfile << "shifting delta will call different sets of significant genes and will result in different FDRs" << endl;
	FDRfile << "you may choose your desired FDR according to this list" << endl << endl;
	for (delta = 0; delta <= 20; delta += 0.05)
	{
		significanceCnt = selectSignificance(delta, individual_index);
		calculate_FDR(delta, significanceCnt, FDR_median, FDR_90percentile, FDR_mean, resultPath);

		FDRfile << "delta = " << delta << ": " << significanceCnt << " genes picked, with FDR(median) = "
			<< FDR_median << " and FDR(mean) = " << FDR_mean << " and FDR(90percentile) = " << FDR_90percentile << endl;

		//cout << "delta = " << delta << ": " << significanceCnt << " genes picked, with FDR(median) = "
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
	for (long geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
	{
		delete [] statistics_permuted_d[geneLoopCnt];
	}
	//delete [] stack;
	delete [] falseGeneCntList;

	return;
}

void output_statistics(string resultPath)
{
	string filename;
	long geneLoopCnt;
	int individualLoopCnt, groupLoopCnt, sampleLoopCnt;
	bool sameDirection, sameDirectionVSexpected;

	/************************************************************************/
	/* Sort All Genes                                                       */
	/************************************************************************/
	for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
	{
		sortKey_gene[geneLoopCnt] = sortList_gene[geneLoopCnt]->statistic_d;
	}
	mergeSort_gene_sort(sortList_gene_Num);


	/************************************************************************/
	/* Diagnostics                                                          */
	/************************************************************************/
//	sprintf(filename, "%sstat/compare_d_expr.txt", resultPath);
//	ofstream d_statistics_outputfile(filename);
//	gene *curGene;
// 	for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
// 	{
// 		curGene = sortList_gene[geneLoopCnt];
// 		sameDirection = true;
// 		sameDirectionVSexpected = true;
// 
// 		d_statistics_outputfile << curGene->chrNM << "\t" << curGene->rangeLow << "\t" <<curGene->rangeHigh << "\t" << curGene->statistic_d << '\t'
// 			<< curGene->statistic_d_expected << "\t" << fabs(curGene->statistic_d - curGene->statistic_d_expected) << "\t" << curGene->statistic_s << "\t" << calculate_foldchange(curGene->meanExpression[1], curGene->meanExpression[2]) << "\t";
// 
// 		for (individualLoopCnt = 1; individualLoopCnt <= num_blocks; ++individualLoopCnt)
// 		{
// 			d_statistics_outputfile << curGene->statistic_d_individual[individualLoopCnt] << "\t" << curGene->statistic_d_expected_individual[individualLoopCnt] << "\t"
// 				<< fabs(curGene->statistic_d_individual[individualLoopCnt] - curGene->statistic_d_expected_individual[individualLoopCnt]) << "\t" << curGene->statistic_s_individual[individualLoopCnt] << "\t";
// 	
// 			if (curGene->statistic_d_individual[individualLoopCnt] * curGene->statistic_d < 0)
// 			{
// 				sameDirection = false;
// 			}
// 			if ((curGene->statistic_d_individual[individualLoopCnt] - curGene->statistic_d_expected_individual[individualLoopCnt]) * (curGene->statistic_d - curGene->statistic_d_expected) < 0)
// 			{
// 				sameDirectionVSexpected = false;
// 			}
// 		}
// 
// 		d_statistics_outputfile << sameDirection << "\t" << sameDirectionVSexpected << "\t" << curGene->meanExpression[1] << "\t" << curGene->meanExpression[2] << "\t" << curGene->meanMedianExpression[1] << "\t" << curGene->meanMedianExpression[2] << "\t" << curGene->noExpressedCnt[1] << "\t" << curGene->noExpressedCnt[2];// << endl;
// 
// 		//output all individual expression. temporary for Aug 26
// 		for (groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
// 			for (sampleLoopCnt = 1; sampleLoopCnt <= SAMPLE_CNT_PER_GROUP; ++sampleLoopCnt)
// 				d_statistics_outputfile << "\t" << curGene->individualMedianExpression[groupLoopCnt][sampleLoopCnt];
// 		d_statistics_outputfile << endl;
// 	}
// 	d_statistics_outputfile.close();
	filename = resultPath + "differential_expression.txt";
	ofstream d_statistics_outputfile(filename.c_str());
	gene *curGene;
	long significantCnt = 0;
	double stat_diffexpr, foldchange;
	d_statistics_outputfile << "chromosome\tposition_start\tposition_end\tstat_diff_expr(|stat_d_expected-stat_d|)\tfold_change\tcoverage_group1\tcoverage_group2\tsignificant" <<endl;
	for (geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
	{
		curGene = sortList_gene[geneLoopCnt];
		stat_diffexpr = fabs(curGene->statistic_d - curGene->statistic_d_expected);
		foldchange = calculate_foldchange(curGene->meanExpression[1], curGene->meanExpression[2]);
		
		d_statistics_outputfile << curGene->chrNM << "\t" << curGene->rangeLow << "\t" <<curGene->rangeHigh << "\t" << stat_diffexpr << "\t" << foldchange << "\t";
		d_statistics_outputfile << curGene->meanMedianExpression[1] << "\t" << curGene->meanMedianExpression[2] << "\t";
		if (stat_diffexpr >= significance_cutoff && (foldchange >= thresh_foldchange_up || foldchange <= thresh_foldchange_down))
		{
			d_statistics_outputfile << "yes";
			++significantCnt;
		}
		else
			d_statistics_outputfile << "no";

#ifdef OUTPUT_PAIRWISE_FOLDCHANGE
		for (int samplecnt = 1; samplecnt <= SAMPLE_CNT_PER_GROUP; ++samplecnt)
		{
			if (curGene->individualExpression[1][samplecnt] > THRESHOLD_MIN_EXPR_COVERAGE && curGene->individualExpression[2][samplecnt] > THRESHOLD_MIN_EXPR_COVERAGE
				|| curGene->individualExpression[1][samplecnt] > 2*THRESHOLD_MIN_EXPR_COVERAGE || curGene->individualExpression[2][samplecnt] > 2*THRESHOLD_MIN_EXPR_COVERAGE)
				d_statistics_outputfile << "\t" << calculate_foldchange(curGene->individualExpression[1][samplecnt], curGene->individualExpression[2][samplecnt]);
			else
				d_statistics_outputfile << "\t1";
		}
#endif

		d_statistics_outputfile << endl;
	}
	d_statistics_outputfile.close();

	cout << "under FDR(median)<=" << false_discovery_rate << ", requiring fold change >=" << thresh_foldchange_up << " or <=" << thresh_foldchange_down << ", select " << significantCnt << " loci with significant change on expression" << endl;
}


void input_normalization_factor(string filename)
{
	ifstream inputfile;
	inputfile.open(filename.c_str());

	for (int loopCnt = 1; loopCnt <= GROUP_NUM * SAMPLE_CNT_PER_GROUP; ++loopCnt)
	{
		inputfile >> NORMALIZATION_FACTOR[loopCnt];
		//cout << NORMALIZATION_FACTOR[loopCnt] << "\t";
	}
	//cout << endl;

	inputfile.close();
	return;
}

void input_normalization_factor()
{
	NORMALIZATION_FACTOR = new double [GROUP_NUM*SAMPLE_CNT_PER_GROUP+1];

	for (int loopCnt = 0; loopCnt <= GROUP_NUM * SAMPLE_CNT_PER_GROUP; ++loopCnt)
	{
		NORMALIZATION_FACTOR[loopCnt] = 1.0;
	}

//	string filename;
// 	filename = "normalization.txt";
// 	input_normalization_factor(filename);

	return;
}


void initialization()
{
	sortKey_gene = new double [sortList_gene_Num+1];
	mergeSort_Larray_gene = new double [sortList_gene_Num+1];
	mergeSort_Rarray_gene = new double [sortList_gene_Num+1];
	mergeSort_LorderedList_gene = new gene* [sortList_gene_Num+1];
	mergeSort_RorderedList_gene = new gene* [sortList_gene_Num+1];

	statistics_expression = new double** [sortList_gene_Num+1];
	for (long geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
	{
		statistics_expression[geneLoopCnt] = new double* [GROUP_NUM+1];
		for (long groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
			statistics_expression[geneLoopCnt][groupLoopCnt] = new double [SAMPLE_CNT_PER_GROUP+1];
	}

	statistics_expected_d = new double [sortList_gene_Num+1];
	statistics_d = new double [sortList_gene_Num+1];
	statistics_s = new double [sortList_gene_Num+1];

	statistics_permuted_d = new double* [sortList_gene_Num+1];

	individualList = new int [num_blocks+1];

	return;
}

void clearall()
{
	delete [] sortKey_gene;
	delete [] mergeSort_Larray_gene;
	delete [] mergeSort_Rarray_gene;
	delete [] mergeSort_LorderedList_gene;
	delete [] mergeSort_RorderedList_gene;

	for (long geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
		delete sortList_gene[geneLoopCnt];
	sortList_gene.clear();

	for (long geneLoopCnt = 1; geneLoopCnt <= sortList_gene_Num; ++geneLoopCnt)
	{
		for (long groupLoopCnt = 1; groupLoopCnt <= GROUP_NUM; ++groupLoopCnt)
			delete [] statistics_expression[geneLoopCnt][groupLoopCnt];
		delete [] statistics_expression[geneLoopCnt];
	}
	delete [] statistics_expression;

	delete [] statistics_expected_d;
	delete [] statistics_d;
	delete [] statistics_s;

	delete [] statistics_permuted_d;

	delete [] individualList;
	delete [] NORMALIZATION_FACTOR;
	
	if (TECHNICAL_REP_NUM != NULL)
		delete [] TECHNICAL_REP_NUM;
	if (TECHNICAL_REP_ARRAY_BASE != NULL)
		delete [] TECHNICAL_REP_ARRAY_BASE;

	return;
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
			else if (parameter.compare("thresh_foldchange_up") == 0)
				config_file >> thresh_foldchange_up;
			else if (parameter.compare("thresh_foldchange_down") == 0)
				config_file >> thresh_foldchange_down;
			else if (parameter.compare("THRESHOLD_MIN_EXPR_COVERAGE") == 0)
				config_file >> THRESHOLD_MIN_EXPR_COVERAGE;

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
	config_filename = "config_testexpr";
#endif

	input_config(config_filename);
	input_normalization_factor();

	/* initialize random seed: */
	srand ( time(NULL) );

	sortList_gene.reserve(DEFAULT_TESTGENE_NUM);
	sortList_gene.push_back(NULL);

	input_testGenes(resultPath);

	initialization();

	//test for all individuals
	resultFileSuffix = "two_group_comparison";
	for (individualLoopCnt = 1; individualLoopCnt <= num_blocks; ++individualLoopCnt)
	{
		individualList[individualLoopCnt] = individualLoopCnt;
	}
	test(resultPath, resultFileSuffix, 0);


	//test for each individual
//	for (int individualLoopCnt = 1; individualLoopCnt <= num_blocks; ++individualLoopCnt)
//	{
//		sprintf(resultFileSuffix, "id_%d", individualLoopCnt);
//		individualList[1] = individualLoopCnt;
//		test(resultPath, resultFileSuffix, individualLoopCnt);
//	}

	output_statistics(resultPath);

	clearall();

	return 0;
}


