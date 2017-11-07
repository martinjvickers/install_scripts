/*    
 *    TranscriptionAnalysis.cpp		
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


#include "TranscriptionAnalysis.h"


/************************************************************************/
/* TestASM                                                              */
/************************************************************************/


testASM::testASM(long ID, int dim, int groupNum, int groupSize)
{
	testID = ID;
	dimension = dim;
	
	individualDistribution = new double** [groupNum+1];
	meanDistribution = new double* [groupNum+1];
	minPathExpr = 0.0;

	for (int groupLoopCnt = 1; groupLoopCnt <= groupNum; ++groupLoopCnt)
	{
		individualDistribution[groupLoopCnt] = new double* [groupSize+1];
		for (int sampleLoopCnt = 1; sampleLoopCnt <= groupSize; ++sampleLoopCnt)
		{
			individualDistribution[groupLoopCnt][sampleLoopCnt] = new double [dimension+1];
		}
		meanDistribution[groupLoopCnt] = new double [dimension+1];
	}
	
	statistic_d = 0.;
	statistic_s = 0.;
}


void testASM::calculate_meanDistribution(int groupNum, int groupSize)
{
	int pathLoopCnt, groupLoopCnt, sampleLoopCnt;
	double groupSum_prop;

	for (pathLoopCnt = 1; pathLoopCnt <= dimension; ++pathLoopCnt)
	{
		for (groupLoopCnt = 1; groupLoopCnt <= groupNum; ++groupLoopCnt)
		{
			groupSum_prop = 0.0;
			for (sampleLoopCnt = 1; sampleLoopCnt <= groupSize; ++sampleLoopCnt)
			{
				groupSum_prop += individualDistribution[groupLoopCnt][sampleLoopCnt][pathLoopCnt];
			}
			meanDistribution[groupLoopCnt][pathLoopCnt] = groupSum_prop / groupSize;
		}
	}

	//normalize the proportions 
	for (groupLoopCnt = 1; groupLoopCnt <= groupNum; ++groupLoopCnt)
	{
		groupSum_prop = 0.0;
		for (pathLoopCnt = 1; pathLoopCnt <= dimension; ++pathLoopCnt)
		{
			groupSum_prop += meanDistribution[groupLoopCnt][pathLoopCnt];
		}		
		if (groupSum_prop == 0.0)
		{
			groupSum_prop = 1.0;
		}
		for (pathLoopCnt = 1; pathLoopCnt <= dimension; ++pathLoopCnt)
		{
			meanDistribution[groupLoopCnt][pathLoopCnt] = meanDistribution[groupLoopCnt][pathLoopCnt] / groupSum_prop;
		}		
	}

	return;
}

testASM::~testASM()
{
	for (int tmpLoop = 1; tmpLoop <= TEST_GROUP_NUM; ++tmpLoop)
	{
		for (int tmpLoop2 = 1; tmpLoop2 <= TEST_GROUP_SIZE; ++tmpLoop2)
		{
			delete [] individualDistribution[tmpLoop][tmpLoop2];
		}
		delete [] meanDistribution[tmpLoop];
		delete [] individualDistribution[tmpLoop];
	}	
	delete [] meanDistribution;
	delete [] individualDistribution;
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

void calculate_meanDistribution(double** distributions, double* meanDist, int distNum, int dimension)
{
	int pathLoopCnt, sampleLoopCnt;
	double sum;

	for (pathLoopCnt = 1; pathLoopCnt <= dimension; ++pathLoopCnt)
	{
		sum = 0.0;
		for (sampleLoopCnt = 1; sampleLoopCnt <= distNum; ++sampleLoopCnt)
		{
			sum += distributions[sampleLoopCnt][pathLoopCnt];
		}
		meanDist[pathLoopCnt] = sum / distNum;
	}

	//normalize the proportions 
	sum = 0.0;
	for (pathLoopCnt = 1; pathLoopCnt <= dimension; ++pathLoopCnt)
	{
		sum += meanDist[pathLoopCnt];
	}		
	if (sum == 0.0)
	{
		sum = 1.0;
	}
	for (pathLoopCnt = 1; pathLoopCnt <= dimension; ++pathLoopCnt)
	{
		meanDist[pathLoopCnt] = meanDist[pathLoopCnt] / sum;
	}	

	return;
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
		mergeSort_Larray_testASM[i] = sortKey_testASM[p + i - 1];
		mergeSort_LorderedList_testASM[i] = sortList_testASM[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray_testASM[j] = sortKey_testASM[q + j];
		mergeSort_RorderedList_testASM[j] = sortList_testASM[q + j];
	}

	mergeSort_Larray_testASM[n1 + 1] = MAX_NUMBER;
	mergeSort_Rarray_testASM[n2 + 1] = MAX_NUMBER;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray_testASM[i] <= mergeSort_Rarray_testASM[j])
		{
			sortKey_testASM[k] = mergeSort_Larray_testASM[i];
			sortList_testASM[k] = mergeSort_LorderedList_testASM[i];

			i++;
		} 
		else
		{
			sortKey_testASM[k] = mergeSort_Rarray_testASM[j];
			sortList_testASM[k] = mergeSort_RorderedList_testASM[j];

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


// KL divergence method for two numbers
double KL_divergence(double x1, double x2)
{
	double temp;
	if (x1 == 0 || x2 == 0)
		temp = 0.0;
	else
		temp = x1 * log(x1/x2) / log(2.0);

	return temp;
}

// array_P: probability distribution for P, sum should equal to 1
// array_Q:probability distribution for Q, sum should equal to 1
// length is the array length for both array_P and array_Q
double calculate_JSD(double *array_P, double *array_Q, long length)
{
	double temp, M, temp_P, temp_Q, propSUM1 = 0.0, propSUM2 = 0.0;
	long i_Loop;

	//normalization
	for (i_Loop = 1; i_Loop <= length; i_Loop++)
	{
		propSUM1 += array_P[i_Loop];
		propSUM2 += array_Q[i_Loop];
	}
	if (propSUM1 == 0)
		propSUM1 = 1.;
	if (propSUM2 == 0)
		propSUM2 = 1.;
	for (i_Loop = 1; i_Loop <= length; i_Loop++)
	{
		array_P[i_Loop] = array_P[i_Loop] / propSUM1;
		array_Q[i_Loop] = array_Q[i_Loop] / propSUM2;
	}

	temp = 0.0;
	for (i_Loop = 1; i_Loop <= length; i_Loop++)
	{
		M = (array_P[i_Loop]+array_Q[i_Loop])/2;
		temp_P = KL_divergence(array_P[i_Loop], M);
		temp_Q = KL_divergence(array_Q[i_Loop], M);
		temp += (temp_P + temp_Q)/2.;
	}

	if (temp < 1e-16)
	{
		temp = 0.0;
	}

	return temp; 
}

/************************************************************************/
/* Statistic for a single gene                                          */
/************************************************************************/

double calculate_single_gene_r_two_class(double *meanProp_group1, double *meanProp_group2, int dimension)
{
	return sqrt(calculate_JSD(meanProp_group1, meanProp_group2, dimension));	
}

double calculate_single_gene_s_two_class(double** individualProp_group1, double* meanProp_group1, int n1, double** individualProp_group2, double* meanProp_group2, int n2, int dimension)
{
	if (n1 + n2 <= 2)
		return 0;

	double SumJSD_group1 = 0.0, SumJSD_group2 = 0.0;
	int sampleLoopCnt;

	for (sampleLoopCnt = 1; sampleLoopCnt <= n1; ++sampleLoopCnt)
	{
		SumJSD_group1 += calculate_JSD(individualProp_group1[sampleLoopCnt], meanProp_group1, dimension);
	}
	for (sampleLoopCnt = 1; sampleLoopCnt <= n2; ++sampleLoopCnt)
	{
		SumJSD_group2 += calculate_JSD(individualProp_group2[sampleLoopCnt], meanProp_group2, dimension);
	}

	return sqrt((1.0/n1 + 1.0/n2) * (SumJSD_group1 + SumJSD_group2) / (n1+n2-2));
}

double calculate_single_gene_d(double** individualProp_group1, double* meanProp_group1, int n1, double** individualProp_group2, double* meanProp_group2, int n2, int dimension, double s0, double value_s)
{
	//if (value_s < 0)
		//value_s has not been calculated, calculate now
		value_s = calculate_single_gene_s_two_class(individualProp_group1, meanProp_group1, n1, individualProp_group2, meanProp_group2, n2, dimension);

	if (value_s + s0 < 1e-15)
	{
		//cout << "dividing a number ~= 0 in calculate_single_gene_d" << endl;
		return 1e-15;
	}
	return calculate_single_gene_r_two_class(meanProp_group1, meanProp_group2, dimension) / (value_s + s0);
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

long permutationTest_within_individual_full(long groupNum, long sampleSize, double s0, string resultPath, double* permuted_d[], string resultFileSuffix, int individualNum, int *individualList)
{
	if (groupNum != 2)
	{
		//only consider permutation over two groups here
		return 0;
	}

	int iLoop, firstGroupIndex, secondGroupIndex, actualPermutationCnt = 0, enumerateRoundCnt = 0, **index, individualCnt, individualIndexBase, sampleLoopCnt, distLoopCnt, curDimension;
	long geneLoopCnt;
	bool **firstGroup, fully_permutated = false;
	double **dist_group1, **dist_group2, *meanDist_1, *meanDist_2, *permuted_value_array = new double [1+sortList_testASM_Num];

//	char permutationTestFilename[1000];
//	sprintf(permutationTestFilename, "%sstat/permutationDiagnostics_s.txt", resultPath, resultFileSuffix);
//	ofstream permutationTestFile(permutationTestFilename);
	
	dist_group1 = new double* [sampleSize+1];
	dist_group2 = new double* [sampleSize+1];

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
			for (geneLoopCnt = 1; geneLoopCnt <= sortList_testASM_Num; ++geneLoopCnt)
			{
				firstGroupIndex = 1; secondGroupIndex = 1; individualIndexBase = 0; curDimension = sortList_testASM[geneLoopCnt]->dimension;
				for (sampleLoopCnt = 1; sampleLoopCnt <= sampleSize; ++sampleLoopCnt)
				{
					dist_group1[sampleLoopCnt] = new double [curDimension+1];
					dist_group2[sampleLoopCnt] = new double [curDimension+1];
				}
				meanDist_1 = new double [curDimension+1];
				meanDist_2 = new double [curDimension+1];

				for (individualCnt = 1; individualCnt <= individualNum; ++individualCnt)
				{
					for (iLoop = 1; iLoop <= TECHNICAL_REP_NUM[individualList[individualCnt]]; ++iLoop)
					{
						if (firstGroup[individualCnt][iLoop] == true)
						{
							for (distLoopCnt = 1; distLoopCnt <= curDimension; ++distLoopCnt)
							{
								dist_group1[firstGroupIndex][distLoopCnt] = sortList_testASM[geneLoopCnt]->individualDistribution[1][individualIndexBase + iLoop][distLoopCnt];
							}
							++firstGroupIndex;
						}
						else
						{
							for (distLoopCnt = 1; distLoopCnt <= curDimension; ++distLoopCnt)
							{
								dist_group2[secondGroupIndex][distLoopCnt] = sortList_testASM[geneLoopCnt]->individualDistribution[1][individualIndexBase + iLoop][distLoopCnt];
							}
							++secondGroupIndex;
						}
						if (firstGroup[individualCnt][iLoop + TECHNICAL_REP_NUM[individualList[individualCnt]]] == true)
						{
							for (distLoopCnt = 1; distLoopCnt <= curDimension; ++distLoopCnt)
							{
								dist_group1[firstGroupIndex][distLoopCnt] = sortList_testASM[geneLoopCnt]->individualDistribution[2][individualIndexBase + iLoop][distLoopCnt];
							}
							++firstGroupIndex;
						}
						else
						{
							for (distLoopCnt = 1; distLoopCnt <= curDimension; ++distLoopCnt)
							{
								dist_group2[secondGroupIndex][distLoopCnt] = sortList_testASM[geneLoopCnt]->individualDistribution[2][individualIndexBase + iLoop][distLoopCnt];
							}
							++secondGroupIndex;
						}
					}
					individualIndexBase += TECHNICAL_REP_NUM[individualList[individualCnt]];
				}

				if (firstGroupIndex != sampleSize + 1 || secondGroupIndex != sampleSize + 1)
				{
					cout << "ERROR: PERMUTATION" << endl;
				}

				calculate_meanDistribution(dist_group1, meanDist_1, sampleSize, curDimension);
				calculate_meanDistribution(dist_group2, meanDist_2, sampleSize, curDimension);

				permuted_value_array[geneLoopCnt] = calculate_single_gene_d(dist_group1, meanDist_1, sampleSize, dist_group2, meanDist_2, sampleSize, curDimension, s0, -1);

				for (sampleLoopCnt = 1; sampleLoopCnt <= sampleSize; ++sampleLoopCnt)
				{
					delete dist_group1[sampleLoopCnt];
					delete dist_group2[sampleLoopCnt];
				}
				delete meanDist_1;
				delete meanDist_2;
			}

			//sort and write the permuted values
			quicksort(permuted_value_array, sortList_testASM_Num);
			for (geneLoopCnt = 1; geneLoopCnt <= sortList_testASM_Num; ++geneLoopCnt)
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
	
	delete [] dist_group1;
	delete [] dist_group2;
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
// 	allpermutedvalueCnt = permutationCnt * sortList_testASM_Num;
	return actualPermutationCnt;
}

double compute_s_expr(double expr)
{
	//given an expression level, return the penalty parameter using a logistic function
	if (expr <= 1)
		return 1;
	//return 2.0 * (1.0 - 1.0 / (1.0 + exp(-1 * log10(expr)))); //P(-t)
	return 2.0 * (1.0 - 1.0 / (1.0 + exp(-1 * expr / 3))); //P(-t)
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
	
	for (geneLoopCnt = 1; geneLoopCnt <= sortList_testASM_Num; ++geneLoopCnt)
	{
		sortKey_testASM[geneLoopCnt] = sortList_testASM[geneLoopCnt]->statistic_s;
	}
	mergeSort_testGene_sort(sortList_testASM_Num);
	for (alphaLoopCnt = 1; alphaLoopCnt <= ALPHANUM; ++alphaLoopCnt)
	{
		//alpha_values[alphaLoopCnt] = sortKey_testGene[long(ceil((sortList_testGene_Num - 1) * (0.05 * (alphaLoopCnt - 1))) + 1)];
		alpha_values[alphaLoopCnt] = calculate_percentile(sortKey_testASM, sortList_testASM_Num, 0.05 * (alphaLoopCnt - 1));
	}


	//compute the 100 quantiles of the si values: q1, ... q100
	double quantiles[101];
	int quantileLoopCnt;

	for (quantileLoopCnt = 1; quantileLoopCnt <= 100; ++quantileLoopCnt)
	{
		//quantiles[quantileLoopCnt] = sortKey_testGene[long(ceil((sortList_testGene_Num - 1) * (0.01 * quantileLoopCnt)) + 1)];
		quantiles[quantileLoopCnt] = calculate_percentile(sortKey_testASM, sortList_testASM_Num, 0.01 * quantileLoopCnt);
	}


	//for every alpha, compute the coefficient of variation of the vj values
	double vj_values[101], cv_values[ALPHANUM+1];
	double* di_alpha = new double [1+sortList_testASM_Num];
	long di_alpha_cnt;
	for (alphaLoopCnt = 1; alphaLoopCnt <= ALPHANUM; ++alphaLoopCnt)
	{
		geneLoopCnt = 1;

		//compute vj
		for (quantileLoopCnt = 1; quantileLoopCnt <= 100; ++quantileLoopCnt)
		{
			di_alpha_cnt = 0;
			//compute d_i^alpha
			for (; geneLoopCnt <= sortList_testASM_Num && sortList_testASM[geneLoopCnt]->statistic_s < quantiles[quantileLoopCnt]; ++geneLoopCnt)
			{
				di_alpha[++di_alpha_cnt] = calculate_single_gene_d(sortList_testASM[geneLoopCnt]->individualDistribution[1], sortList_testASM[geneLoopCnt]->meanDistribution[1], sampleSize, 
					sortList_testASM[geneLoopCnt]->individualDistribution[2], sortList_testASM[geneLoopCnt]->meanDistribution[2], sampleSize, sortList_testASM[geneLoopCnt]->dimension, alpha_values[alphaLoopCnt], sortList_testASM[geneLoopCnt]->statistic_s);
			}

			if (di_alpha_cnt > 0)
			{
				vj_values[quantileLoopCnt] = calculate_mad(di_alpha, di_alpha_cnt) / 0.64;
			} 
			else
			{
				vj_values[quantileLoopCnt] = 0.0;
			}
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
	outputfilename = resultPath + "stat/asm/d_permutation_" + resultFileSuffix + ".txt";
	outputfile.open(outputfilename.c_str());
	for (long geneLoopCnt = 1; geneLoopCnt <= sortList_testASM_Num; ++geneLoopCnt)
	{
		for (long loopCnt = 1; loopCnt <= permutationCnt; ++loopCnt)
		{
			++totalCnt;
			outputfile << permuted_d[geneLoopCnt][loopCnt] << endl; 
		}
	}
	outputfile.close();
	
	//output count of permuted values, genes, and permutations
	outputfilename = resultPath + "stat/asm/d_permutation_cnt_" + resultFileSuffix + ".txt";
	outputfile.open(outputfilename.c_str());
	outputfile << totalCnt << "\t" << sortList_testASM_Num << "\t" << permutationCnt << endl;
	outputfile.close();

	return;
}

void initialize(long array_size)
{
	sortList_testASM = new testASM* [1+array_size];
	sortKey_testASM = new double [1+array_size];
	mergeSort_Larray_testASM = new double [1+array_size];
	mergeSort_Rarray_testASM = new double [1+array_size];
	mergeSort_LorderedList_testASM = new testASM* [1+array_size];
	mergeSort_RorderedList_testASM = new testASM* [1+array_size];
}

void clean_up()
{
	delete [] sortList_testASM;
	delete [] sortKey_testASM;
	delete [] mergeSort_Larray_testASM;
	delete [] mergeSort_Rarray_testASM;
	delete [] mergeSort_LorderedList_testASM;
	delete [] mergeSort_RorderedList_testASM;
}

//sampleSize is the number of samples in each group, which should be the sum of #tech-replicates of every individual; also, assume sampleSize is fixed for each group, for the convenience of permutation test
bool diff_expression_analysis_two_class(double**** pathDistribution, long *dimensionList, double *minPathExprList, double *statistic_d, double *statistic_s, double *statistics_expected_d, double* statistics_permuted_d[], long testASMNum, int groupNum, int individualNum, int *individualList, long sampleSize, string resultPath, string resultFileSuffix)
{
	if (groupNum != 2)
	{
		return false;
	}

	initialize(testASMNum);

	TEST_GROUP_NUM = groupNum; TEST_GROUP_SIZE = sampleSize;

	long ASMloopCnt, permutationCnt_actual, cnt_small_s = 0;;
	int sampleLoopCnt, groupLoopCnt, distLoopCnt;
	sortList_testASM_Num = 0;


	testASM *newTestASM;
	for (ASMloopCnt = 1; ASMloopCnt <= testASMNum; ++ASMloopCnt)
	{
		newTestASM = new testASM(ASMloopCnt, dimensionList[ASMloopCnt], groupNum, sampleSize);
		for (groupLoopCnt = 1; groupLoopCnt <= groupNum; ++groupLoopCnt)
		{
			for (sampleLoopCnt = 1; sampleLoopCnt <= sampleSize; ++sampleLoopCnt)
			{
				for (distLoopCnt = 1; distLoopCnt <= newTestASM->dimension; ++distLoopCnt)
				{
					newTestASM->individualDistribution[groupLoopCnt][sampleLoopCnt][distLoopCnt] = pathDistribution[ASMloopCnt][groupLoopCnt][sampleLoopCnt][distLoopCnt];
				}
			}
		}
		newTestASM->calculate_meanDistribution(groupNum, sampleSize);
		newTestASM->minPathExpr = minPathExprList[ASMloopCnt];
		sortList_testASM[++sortList_testASM_Num] = newTestASM;
	}


	double s0, s_max = 0.0;
	//compute s
	for (ASMloopCnt = 1; ASMloopCnt <= sortList_testASM_Num; ++ASMloopCnt)
	{
		sortList_testASM[ASMloopCnt]->statistic_s = calculate_single_gene_s_two_class(sortList_testASM[ASMloopCnt]->individualDistribution[1], sortList_testASM[ASMloopCnt]->meanDistribution[1],
			sampleSize, sortList_testASM[ASMloopCnt]->individualDistribution[2], sortList_testASM[ASMloopCnt]->meanDistribution[2], sampleSize, sortList_testASM[ASMloopCnt]->dimension);
		if (sortList_testASM[ASMloopCnt]->statistic_s > s_max)
			s_max = sortList_testASM[ASMloopCnt]->statistic_s;
		if (sortList_testASM[ASMloopCnt]->statistic_s < 1e-3)
		{
			//cout << "debug";
			sortList_testASM[ASMloopCnt]->statistic_s = 1e-3; //MAX_NUMBER/2;
			++cnt_small_s;
		}
	}
	//add expression term using min path expression
	for (ASMloopCnt = 1; ASMloopCnt <= sortList_testASM_Num; ++ASMloopCnt)
	{
		sortList_testASM[ASMloopCnt]->statistic_s += compute_s_expr(sortList_testASM[ASMloopCnt]->minPathExpr) * s_max * 1.3;
	}
	//compute s0
	s0 = compute_s0(sampleSize);
	//cout << "s0 = " << s0 << "\t";

	//compute di
	for (ASMloopCnt = 1; ASMloopCnt <= sortList_testASM_Num; ++ASMloopCnt)
	{
		sortList_testASM[ASMloopCnt]->statistic_d = calculate_single_gene_d(sortList_testASM[ASMloopCnt]->individualDistribution[1], sortList_testASM[ASMloopCnt]->meanDistribution[1],
			sampleSize, sortList_testASM[ASMloopCnt]->individualDistribution[2], sortList_testASM[ASMloopCnt]->meanDistribution[2], sampleSize, sortList_testASM[ASMloopCnt]->dimension, s0, sortList_testASM[ASMloopCnt]->statistic_s);
	}


	//permutation 
	permutationCnt_actual = permutationTest_within_individual_full(groupNum, sampleSize, s0, resultPath, statistics_permuted_d, resultFileSuffix, individualNum, individualList);

	
	//output permutation result
	output_permuted_d(resultPath, permutationCnt_actual, statistics_permuted_d, resultFileSuffix);

	//store the statistics
	for (ASMloopCnt = 1; ASMloopCnt <= sortList_testASM_Num; ++ASMloopCnt)
	{
		statistic_d[sortList_testASM[ASMloopCnt]->testID] = sortList_testASM[ASMloopCnt]->statistic_d;
		statistic_s[sortList_testASM[ASMloopCnt]->testID] = calculate_single_gene_s_two_class(sortList_testASM[ASMloopCnt]->individualDistribution[1], sortList_testASM[ASMloopCnt]->meanDistribution[1],
			sampleSize, sortList_testASM[ASMloopCnt]->individualDistribution[2], sortList_testASM[ASMloopCnt]->meanDistribution[2], sampleSize, sortList_testASM[ASMloopCnt]->dimension);
		
		statistics_expected_d[ASMloopCnt] = calculate_mean(statistics_permuted_d[ASMloopCnt], permutationCnt_actual);
	}

	for (ASMloopCnt = 1; ASMloopCnt <= testASMNum; ++ASMloopCnt)
		delete sortList_testASM[ASMloopCnt];
	
	clean_up();

	return true;
}







