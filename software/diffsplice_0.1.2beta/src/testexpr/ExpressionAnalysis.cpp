/*    
 *    ExpressionAnalysis.cpp		
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


#include "ExpressionAnalysis.h"



testGene::testGene(long permutationCnt, long ID, int groupNum, int groupSize)
{
	testID = ID;
	expression = new double* [groupNum+1];
	for (int iLoop = 1; iLoop <= groupNum; ++iLoop)
	{
		expression[iLoop] = new double [groupSize+1];
	}
	statistic_d = 0;
	statistic_s = 0;
//	statistic_d_permutation = new double [permutationCnt+1];
}

testGene::~testGene()
{
	for (int iLoop = 1; iLoop <= TEST_GROUP_NUM; ++iLoop)
	{
		delete [] expression[iLoop];
	}
	delete [] expression;
//	delete [] statistic_d_permutation;
}




/************************************************************************/
/* Basic functions                                                      */
/************************************************************************/


long partition(long p,long r, double *sortArray)
{
	long i, j;
	double x, tmp;

	//randomized partition
	i = long(p + (double)rand()/ (RAND_MAX) * (r - p));

	if (sortArray[r] != sortArray[i])
	{
		tmp = sortArray[r];
		sortArray[r] = sortArray[i];
		sortArray[i] = tmp;
	}

	x = sortArray[r];
	i = p - 1;

	for(j = p; j < r; j++)
	{
		if (sortArray[j] <= x)
		{
			i++;

			if (sortArray[i] != sortArray[j])
			{
				tmp = sortArray[j];
				sortArray[j] = sortArray[i];
				sortArray[i] = tmp;
			}
		}
	}

	if (sortArray[r] != sortArray[i+1])
	{
		tmp = sortArray[r];
		sortArray[r] = sortArray[i+1];
		sortArray[i+1]=tmp;
	}

	return i+1;
}


void quicksort(double *sortArray, long length)
{
	long top = 0, p, r, q;

	quicksortStack = new long [10+length*2];

	quicksortStack[top++] = 1;
	quicksortStack[top++] = length;

	while (top != 0)
	{
		r = quicksortStack[--top];
		p = quicksortStack[--top];

		if(p>=r)
			continue;

		q = partition(p, r, sortArray);

		quicksortStack[top++] = p;
		quicksortStack[top++] = q - 1;

		quicksortStack[top++] = q + 1;
		quicksortStack[top++] = r;
	}

	delete [] quicksortStack;

	return;
}


double calculate_mean(double *dataarray, long length)
{
	if (length < 1)
	{
		cout << "array length < 1 in calculate_mean";
		return 0.0;
	}

	double sum = 0.0;
	for (long tmpCnt = 1; tmpCnt <= length; ++tmpCnt)
	{
		sum += dataarray[tmpCnt];
	}
	return sum / length;
}

double calculate_percentile(double *dataarray, long length, double percentile_0to1)
{
	double percentileRank, rangeStart, rangeEnd, percentileValue = 0.0;

	if (percentile_0to1 < 0 || percentile_0to1 > 1)
	{
		cout << "ERROR: percentile should be 0 to 1 in calculate_percentile" << endl;
		exit(1);
	}

	percentileRank = percentile_0to1 * (length - 1) + 1;

	rangeStart = dataarray[long(floor(percentileRank))];
	rangeEnd   = dataarray[long(ceil(percentileRank))];

	percentileValue = rangeStart + (percentileRank - floor(percentileRank)) * (rangeEnd - rangeStart);

	return percentileValue;
}

double calculate_SSE(double *dataarray, long length)
{
	double SSE = 0.0, groupmean = calculate_mean(dataarray, length);
	for (long tmpCnt = 1; tmpCnt <= length; ++tmpCnt)
	{
		SSE += pow(dataarray[tmpCnt] - groupmean, 2);
	}
	return SSE;
}

long calculate_factorial(long n)
{
	if (n == 0)
		return 1;
	else
		return n * calculate_factorial(n-1);
}

long choose(long n, long k)
{
	return calculate_factorial(n) / calculate_factorial(k) / calculate_factorial(n-k);
}

long choose_non_recursive(long n, long k)
{
	if (k > n || k < 0)
		return 0;
	if (k > n - k)
		k = n - k;

	double comb = 1;
	for (long tmpLoop = 1; tmpLoop <= k; ++tmpLoop)
	{
		comb *= (n - (k - tmpLoop));
		comb = comb / tmpLoop;
	}

	return long(floor(comb + 0.5));
}

void merge_testGene_sort(long p, long q, long r)
{
	long n1, n2, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	for (i = 1; i <= n1; i++)
	{
		mergeSort_Larray_testGene[i] = sortKey_testGene[p + i - 1];
		mergeSort_LorderedList_testGene[i] = sortList_testGene[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray_testGene[j] = sortKey_testGene[q + j];
		mergeSort_RorderedList_testGene[j] = sortList_testGene[q + j];
	}

	mergeSort_Larray_testGene[n1 + 1] = MAX_NUMBER;
	mergeSort_Rarray_testGene[n2 + 1] = MAX_NUMBER;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray_testGene[i] <= mergeSort_Rarray_testGene[j])
		{
			sortKey_testGene[k] = mergeSort_Larray_testGene[i];
			sortList_testGene[k] = mergeSort_LorderedList_testGene[i];

			i++;
		} 
		else
		{
			sortKey_testGene[k] = mergeSort_Rarray_testGene[j];
			sortList_testGene[k] = mergeSort_RorderedList_testGene[j];

			j++;
		}
	}

	return;
}


void mergeSort_testGene_sort(long sortList_size)
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
			merge_testGene_sort(i, i + m - 1, r);
			i = i + 2 * m;
		}

		m = m * 2;
	}

	return;
}



/************************************************************************/
/* Statistic for a single gene                                          */
/************************************************************************/

double calculate_single_gene_r_two_class(double *expr_group1, int n1, double *expr_group2, int n2)
{
	return calculate_mean(expr_group2, n2) - calculate_mean(expr_group1, n1);	
}

double calculate_single_gene_s_two_class(double *expr_group1, int n1, double *expr_group2, int n2)
{
	if (n1 + n2 <= 2)
		return 0;
	return sqrt((1.0/n1 + 1.0/n2) * (calculate_SSE(expr_group1, n1) + calculate_SSE(expr_group2, n2)) / (n1+n2-2));
}

double calculate_single_gene_d(double *expr_group1, int n1, double *expr_group2, int n2, double s0)
{
	double value_s = calculate_single_gene_s_two_class(expr_group1, n1, expr_group2, n2);
	if (value_s + s0 < 1e-15)
	{
		//cout << "dividing a number ~= 0 in calculate_single_gene_d" << endl;
		return 1e-15;
	}
	return calculate_single_gene_r_two_class(expr_group1, n1, expr_group2, n2) / (value_s + s0);
}

double calculate_mad(double *dataArray, long arrayLength)
{
	//calculate the median absolute deviation from the median
	double median = 0.0, *abs_deviation, mad_value;
	abs_deviation = new double [arrayLength+1];

	quicksort(dataArray, arrayLength);
	//median = dataArray[long(ceil(arrayLength/2.))];
	median = calculate_percentile(dataArray, arrayLength, 0.5);

	for (long loopCnt = 1; loopCnt <= arrayLength; ++loopCnt)
	{
		abs_deviation[loopCnt] = fabs(dataArray[loopCnt] - median);
	}

	quicksort(abs_deviation, arrayLength);
	//mad_value = abs_deviation[long(ceil(arrayLength/2.))];
	mad_value = calculate_percentile(abs_deviation, arrayLength, 0.5);

	delete [] abs_deviation;
	return mad_value;
}

double calculate_coefficient_variation(double *dataArray, long arrayLength)
{
	//calculate the coefficient variation;

	return sqrt(calculate_SSE(dataArray, arrayLength)/(arrayLength-1)) / fabs(calculate_mean(dataArray, arrayLength));
}


/************************************************************************/
/* Statistic for all genes                                              */
/************************************************************************/

long permutationTest_within_individual(long groupNum, long sampleSize, double s0, string resultPath, double* permuted_d[], string resultFileSuffix, int individualNum, int *individualList)
{
	if (groupNum != 2)
	{
		//only consider permutation over two groups here
		return 0;
	}

	int iLoop, firstGroupIndex, secondGroupIndex, actualPermutationCnt = 0, enumerateRoundCnt = 0, **index, individualCnt, individualIndexBase;
	long geneLoopCnt, sampleLoopCnt;
	bool **firstGroup, fully_permutated = false;
	double *expr_group1, *expr_group2, *permuted_value_array = new double [1+sortList_testGene_Num];

//	char permutationTestFilename[1000];
//	sprintf(permutationTestFilename, "%sstat/permutationDiagnostics_s.txt", resultPath, resultFileSuffix);
//	ofstream permutationTestFile(permutationTestFilename);


	expr_group1 = new double [sampleSize+1];
	expr_group2 = new double [sampleSize+1];

	firstGroup = new bool* [individualNum+1];
	index = new int* [individualNum+1];
	for (individualCnt = 1; individualCnt <= individualNum; ++individualCnt)
	{
		firstGroup[individualCnt] = new bool [2*TECHNICAL_REP_NUM[individualList[individualCnt]]+1];
		index[individualCnt] = new int [TECHNICAL_REP_NUM[individualList[individualCnt]]+1];
	}

	for (individualCnt = 1; individualCnt <= individualNum; ++individualCnt)
	{
		for (iLoop = 1; iLoop <= TECHNICAL_REP_NUM[individualList[individualCnt]]; ++iLoop)
		{
			index[individualCnt][iLoop] = iLoop;
		}
	}

	while (actualPermutationCnt < permutationCnt && fully_permutated == false)
	{
		++enumerateRoundCnt;

		if (permType == perm_random || permType == perm_full || permType == perm_reducedfull && enumerateRoundCnt % permutation_limit_ratio == 0)
		{
			++actualPermutationCnt;

			//build the indicator matrix
			if (permType == perm_full || permType == perm_reducedfull)
			{
				for (individualCnt = 1; individualCnt <= individualNum; ++individualCnt)
				{
					for (iLoop = 1; iLoop <= 2 * TECHNICAL_REP_NUM[individualList[individualCnt]]; ++iLoop)
					{
						firstGroup[individualCnt][iLoop] = false;
					}
					for (iLoop = 1; iLoop <= TECHNICAL_REP_NUM[individualList[individualCnt]]; ++iLoop)
					{
						firstGroup[individualCnt][index[individualCnt][iLoop]] = true;
					}
				}
			} 
			else if (permType == perm_random)
			{
				for (individualCnt = 1; individualCnt <= individualNum; ++individualCnt)
				{
					for (iLoop = 1; iLoop <= 2 * TECHNICAL_REP_NUM[individualList[individualCnt]]; ++iLoop)
					{
						firstGroup[individualCnt][iLoop] = false;
					}
					int firstGroupCnt = 0, tmpindex;
					while (firstGroupCnt < TECHNICAL_REP_NUM[individualList[individualCnt]])
					{
						tmpindex = rand() % (2 * TECHNICAL_REP_NUM[individualList[individualCnt]]) + 1;
						if (firstGroup[individualCnt][tmpindex] == false)
						{
							firstGroup[individualCnt][tmpindex] = true;
							++firstGroupCnt;
						}						
					}					
				}
			}
			else
			{
				cout << "error: unrecognized permutation type" << endl;
				exit(1);
			}


			//calculate the statistic d
			for (geneLoopCnt = 1; geneLoopCnt <= sortList_testGene_Num; ++geneLoopCnt)
			{
				firstGroupIndex = 1; secondGroupIndex = 1; individualIndexBase = 0;

				for (individualCnt = 1; individualCnt <= individualNum; ++individualCnt)
				{
					for (iLoop = 1; iLoop <= TECHNICAL_REP_NUM[individualList[individualCnt]]; ++iLoop)
					{
						if (firstGroup[individualCnt][iLoop] == true)
						{
							expr_group1[firstGroupIndex] = sortList_testGene[geneLoopCnt]->expression[1][individualIndexBase + iLoop];
							++firstGroupIndex;
						}
						else
						{
							expr_group2[secondGroupIndex] = sortList_testGene[geneLoopCnt]->expression[1][individualIndexBase + iLoop];
							++secondGroupIndex;
						}
						if (firstGroup[individualCnt][iLoop + TECHNICAL_REP_NUM[individualList[individualCnt]]] == true)
						{
							expr_group1[firstGroupIndex] = sortList_testGene[geneLoopCnt]->expression[2][individualIndexBase + iLoop];
							++firstGroupIndex;
						}
						else
						{
							expr_group2[secondGroupIndex] = sortList_testGene[geneLoopCnt]->expression[2][individualIndexBase + iLoop];
							++secondGroupIndex;
						}
					}
					individualIndexBase += TECHNICAL_REP_NUM[individualList[individualCnt]];
				}

				if (firstGroupIndex != sampleSize + 1 || secondGroupIndex != sampleSize + 1)
				{
					cout << "ERROR: PERMUTATION" << endl;
				}

				permuted_value_array[geneLoopCnt] = calculate_single_gene_d(expr_group1, sampleSize, expr_group2, sampleSize, s0);
			}

			//sort and write the permuted values
			quicksort(permuted_value_array, sortList_testGene_Num);
			for (geneLoopCnt = 1; geneLoopCnt <= sortList_testGene_Num; ++geneLoopCnt)
			{
				permuted_d[geneLoopCnt][actualPermutationCnt] = permuted_value_array[geneLoopCnt];
				//			permutationTestFile << permuted_d[geneLoopCnt][roundCnt] << "\t";
			}
			//		permutationTestFile << endl;
		}

		if (permType == perm_full || permType == perm_reducedfull)
		{
			//update the index vector!!
			index[individualNum][TECHNICAL_REP_NUM[individualList[individualNum]]] += 1;
			if (index[individualNum][TECHNICAL_REP_NUM[individualList[individualNum]]] > 2*TECHNICAL_REP_NUM[individualList[individualNum]])
			{
				for (individualCnt = individualNum; individualCnt >= 1 && index[individualCnt][TECHNICAL_REP_NUM[individualList[individualCnt]]] > 2 * TECHNICAL_REP_NUM[individualList[individualCnt]]; --individualCnt)
				{
					for (iLoop = TECHNICAL_REP_NUM[individualList[individualCnt]] - 1; iLoop >= 1 && index[individualCnt][iLoop] >= TECHNICAL_REP_NUM[individualList[individualCnt]] + iLoop; --iLoop);
					if (iLoop >= 1)
					{
						break;
					}
					else
					{
						if (individualCnt > 1)
							index[individualCnt-1][TECHNICAL_REP_NUM[individualList[individualCnt]-1]] += 1;
					}
				}

				if (individualCnt >= 1)
				{
					//for this individual, set the following bits 
					if (iLoop >= 1)
					{
						index[individualCnt][iLoop] += 1;
						for (++iLoop; iLoop <= TECHNICAL_REP_NUM[individualList[individualCnt]]; ++iLoop)
						{
							index[individualCnt][iLoop] = index[individualCnt][iLoop-1] + 1;
						}
					}

					//for the following individuals, reset to 1,2,3,...
					for (++individualCnt; individualCnt <= individualNum; ++individualCnt)
					{
						for (iLoop = 1; iLoop <= TECHNICAL_REP_NUM[individualList[individualCnt]]; ++iLoop)
						{
							index[individualCnt][iLoop] = iLoop;
						}
					}
				} 
				else
				{
					//all combination enumerated
					fully_permutated = true;
				}
			}
		}
	}

	//for debug
	if (actualPermutationCnt != permutationCnt)
	{
		cout << "abnormal case, actualPermutationCnt < permutationCnt" << endl;
	}

	delete [] expr_group1;
	delete [] expr_group2;
	for (individualCnt = 1; individualCnt <= individualNum; ++individualCnt)
	{
		delete [] firstGroup[individualCnt];
		delete [] index[individualCnt];
	}
	delete [] firstGroup;
	delete [] index;

	delete [] permuted_value_array;

//	permutationTestFile.close();

// 	permutationCnt = actualPermutationCnt;
// 	allpermutedvalueCnt = permutationCnt * sortList_testGene_Num;
	return actualPermutationCnt;
}


double compute_s0(long sampleSize)
{
	//compute value of s0
	double s0 = 0.0;
	long geneLoopCnt;

	//compute s^\alpha for \alpha from 0 to 1.0 in step size = 0.05
	const int ALPHANUM = 21;
	int alphaLoopCnt;
	double alpha_values[50];
	
	for (geneLoopCnt = 1; geneLoopCnt <= sortList_testGene_Num; ++geneLoopCnt)
	{
		sortKey_testGene[geneLoopCnt] = sortList_testGene[geneLoopCnt]->statistic_s;
	}
	mergeSort_testGene_sort(sortList_testGene_Num);
	for (alphaLoopCnt = 1; alphaLoopCnt <= ALPHANUM; ++alphaLoopCnt)
	{
		//alpha_values[alphaLoopCnt] = sortKey_testGene[long(ceil((sortList_testGene_Num - 1) * (0.05 * (alphaLoopCnt - 1))) + 1)];
		alpha_values[alphaLoopCnt] = calculate_percentile(sortKey_testGene, sortList_testGene_Num, 0.05 * (alphaLoopCnt - 1));
	}


	//compute the 100 quantiles of the si values: q1, ... q100
	double quantiles[101];
	int quantileLoopCnt;

	for (quantileLoopCnt = 1; quantileLoopCnt <= 100; ++quantileLoopCnt)
	{
		//quantiles[quantileLoopCnt] = sortKey_testGene[long(ceil((sortList_testGene_Num - 1) * (0.01 * quantileLoopCnt)) + 1)];
		quantiles[quantileLoopCnt] = calculate_percentile(sortKey_testGene, sortList_testGene_Num, 0.01 * quantileLoopCnt);
	}


	//for every alpha, compute the coefficient of variation of the vj values
	double *di_alpha = new double [1+sortList_testGene_Num], vj_values[101], cv_values[ALPHANUM+1];
	long di_alpha_cnt;
	for (alphaLoopCnt = 1; alphaLoopCnt <= ALPHANUM; ++alphaLoopCnt)
	{
		geneLoopCnt = 1;

		//compute vj
		for (quantileLoopCnt = 1; quantileLoopCnt <= 100; ++quantileLoopCnt)
		{
			di_alpha_cnt = 0;
			//compute d_i^alpha
			for (; geneLoopCnt <= sortList_testGene_Num && sortList_testGene[geneLoopCnt]->statistic_s < quantiles[quantileLoopCnt]; ++geneLoopCnt)
			{
				di_alpha[++di_alpha_cnt] = calculate_single_gene_d(sortList_testGene[geneLoopCnt]->expression[1], sampleSize, sortList_testGene[geneLoopCnt]->expression[2], sampleSize, alpha_values[alphaLoopCnt]);
			}

			vj_values[quantileLoopCnt] = calculate_mad(di_alpha, di_alpha_cnt) / 0.64;
		}

		//compute coefficient of variation of the vj values
		cv_values[alphaLoopCnt] = calculate_coefficient_variation(vj_values, 100);
	}
	delete [] di_alpha;

	//find the min of the c_v's
	double min_cv = MAX_NUMBER;
	int min_cv_index = 1; 
	for (alphaLoopCnt = 1; alphaLoopCnt <= ALPHANUM; ++alphaLoopCnt)
	{
		if (cv_values[alphaLoopCnt] < min_cv)
		{
			min_cv = cv_values[alphaLoopCnt];
			min_cv_index = alphaLoopCnt;
		}
	}
	s0 = alpha_values[min_cv_index];

	return s0;
}



//output permuted d values
void output_permuted_d(string resultPath, long permutationCnt, double* permuted_d[], string resultFileSuffix)
{
	string outputfilename;
	ofstream outputfile;
	long totalCnt = 0; //should be testGeneNum * permutationCnt

	//output all permuted d values
	outputfilename = resultPath + "stat/expr/d_permutation_" + resultFileSuffix + ".txt";
	outputfile.open(outputfilename.c_str());
	for (long geneLoopCnt = 1; geneLoopCnt <= sortList_testGene_Num; ++geneLoopCnt)
	{
		for (long loopCnt = 1; loopCnt <= permutationCnt; ++loopCnt)
		{
			++totalCnt;
			outputfile << permuted_d[geneLoopCnt][loopCnt] << endl; 
		}
	}
	outputfile.close();
	
	//output count of permuted values, genes, and permutations
	outputfilename = resultPath + "stat/expr/d_permutation_cnt_" + resultFileSuffix + ".txt";
	outputfile.open(outputfilename.c_str());
	outputfile << totalCnt << "\t" << sortList_testGene_Num << "\t" << permutationCnt << endl;
	outputfile.close();

	return;
}


void initialize(long array_size)
{
	sortList_testGene = new testGene* [1+array_size];
	sortKey_testGene = new double [1+array_size];
	mergeSort_Larray_testGene = new double [1+array_size];
	mergeSort_Rarray_testGene = new double [1+array_size];
	mergeSort_LorderedList_testGene = new testGene* [1+array_size];
	mergeSort_RorderedList_testGene = new testGene* [1+array_size];
}

void clean_up()
{
	delete [] sortList_testGene;
	delete [] sortKey_testGene;
	delete [] mergeSort_Larray_testGene;
	delete [] mergeSort_Rarray_testGene;
	delete [] mergeSort_LorderedList_testGene;
	delete [] mergeSort_RorderedList_testGene;
}


//data format: 3-dim array, [geneID] [groupID] [sampleID]
//sampleSize is the number of samples in each group, which should be the sum of #tech-replicates of every individual; also, assume sampleSize is fixed for each group, for the convenience of permutation test
bool diff_expression_analysis_two_class(double*** expression, double *statistic_d, double *statistic_s, double *statistics_expected_d, double* statistics_permuted_d[], long testGeneNum, int groupNum, int individualNum, int *individualList, long sampleSize, string resultPath, string resultFileSuffix)
{
	if (groupNum != 2)
	{
		return false;
	}

	initialize(testGeneNum);

	TEST_GROUP_NUM = groupNum;

	long geneLoopCnt, sampleLoopCnt, groupLoopCnt, permutationCnt_actual;
	int individualCnt;
	sortList_testGene_Num = 0;

// #ifdef PERMUTATION_WITHIN_INDIVIDUAL
// 	for (individualCnt = 1; individualCnt <= INDIVIDUAL_NUM; ++individualCnt)
// 	{
// 		permutationCnt *= choose(groupNum * TECHNICAL_REP_NUM[individualCnt], TECHNICAL_REP_NUM[individualCnt]); 
// 	}
// #endif

	//generate all test genes
	testGene *newTestGene;
	for (geneLoopCnt = 1; geneLoopCnt <= testGeneNum; ++geneLoopCnt)
	{
		newTestGene = new testGene(permutationCnt, geneLoopCnt, groupNum, sampleSize);
		for (groupLoopCnt = 1; groupLoopCnt <= groupNum; ++groupLoopCnt)
		{
			for (sampleLoopCnt = 1; sampleLoopCnt <= sampleSize; ++sampleLoopCnt)
			{
				newTestGene->expression[groupLoopCnt][sampleLoopCnt] = expression[geneLoopCnt][groupLoopCnt][sampleLoopCnt];
			}
		}
		sortList_testGene[++sortList_testGene_Num] = newTestGene;
	}



	double s0;
	//compute s0
	for (geneLoopCnt = 1; geneLoopCnt <= sortList_testGene_Num; ++geneLoopCnt)
	{
		sortList_testGene[geneLoopCnt]->statistic_s = calculate_single_gene_s_two_class(sortList_testGene[geneLoopCnt]->expression[1], sampleSize, sortList_testGene[geneLoopCnt]->expression[2], sampleSize);
		if (sortList_testGene[geneLoopCnt]->statistic_s == 0)
		{
			//cout << "debug";
			sortList_testGene[geneLoopCnt]->statistic_s = MAX_NUMBER/2;
		}
	}
	s0 = compute_s0(sampleSize);
	//cout << "s0 = " << s0 << "\t";

	//compute di
	for (geneLoopCnt = 1; geneLoopCnt <= sortList_testGene_Num; ++geneLoopCnt)
	{
		sortList_testGene[geneLoopCnt]->statistic_d = calculate_single_gene_d(sortList_testGene[geneLoopCnt]->expression[1], sampleSize, sortList_testGene[geneLoopCnt]->expression[2], sampleSize, s0);
	}


	//permutation 
	permutationCnt_actual = permutationTest_within_individual(groupNum, sampleSize, s0, resultPath, statistics_permuted_d, resultFileSuffix, individualNum, individualList);

	
	//output permutation result
	output_permuted_d(resultPath, permutationCnt_actual, statistics_permuted_d, resultFileSuffix);

	//store the statistics
	for (geneLoopCnt = 1; geneLoopCnt <= sortList_testGene_Num; ++geneLoopCnt)
	{
		statistic_d[sortList_testGene[geneLoopCnt]->testID] = sortList_testGene[geneLoopCnt]->statistic_d;
		statistic_s[sortList_testGene[geneLoopCnt]->testID] = sortList_testGene[geneLoopCnt]->statistic_s;
		
		statistics_expected_d[geneLoopCnt] = calculate_mean(statistics_permuted_d[geneLoopCnt], permutationCnt_actual);
	}

	for (geneLoopCnt = 1; geneLoopCnt <= testGeneNum; ++geneLoopCnt)
	{
		delete sortList_testGene[geneLoopCnt];
	}

	clean_up();

	return true;
}







