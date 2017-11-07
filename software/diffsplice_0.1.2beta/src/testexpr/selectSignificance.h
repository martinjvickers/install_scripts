/*    
 *    selectSignificance.h		
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

#include "Configuration.h"

/************************************************************************/
/* GENERAL SETTINGS                                                     */
/************************************************************************/

int GROUP_NUM = 0;
int INDIVIDUAL_NUM = 0; //number of individuals
int SAMPLE_CNT_PER_GROUP = 0; //number of samples per group
int TOTAL_SAMPLE_CNT = 0; //total number of samples

int num_blocks = 0; //number of permutation blocks (individuals)
int* TECHNICAL_REP_NUM = NULL; //number of technical replicates of each individual
int* TECHNICAL_REP_ARRAY_BASE = NULL; //array base of technical replicates of each individual


double THRESHOLD_MIN_EXPR_COVERAGE = 10; //at least one group must have mean expression larger than this threshold to make differential expression
double false_discovery_rate = 0.01;
double thresh_foldchange_up = 2;
double thresh_foldchange_down = 0.5;

/************************************************************************/
/* FIXED SETTINGS	                                                    */
/************************************************************************/

const long DEFAULT_TESTGENE_NUM = 100000;



/************************************************************************/
/* DEFINITIONS & DECLARATIONS                                           */
/************************************************************************/



class gene
{
public:
	//basic information
	long testID;
	string chrNM;
	long rangeLow;
	long rangeHigh;
	double** individualExpression; //[GROUP_NUM+1][SAMPLE_CNT_PER_GROUP+1];
	double* meanExpression; //[GROUP_NUM + 1];
	double** individualMedianExpression; //[GROUP_NUM+1][SAMPLE_CNT_PER_GROUP+1];
	double* meanMedianExpression; //[GROUP_NUM + 1];
	int* noExpressedCnt; //[GROUP_NUM + 1];
	int* expressedCnt; //[GROUP_NUM + 1];

	//statistics of all-together
	double statistic_d;
	double statistic_s;
	double statistic_d_expected;
	
	//statistics of every individual
	double* statistic_d_individual; //[num_blocks + 1];
	double* statistic_s_individual; //[num_blocks + 1];
	double* statistic_d_expected_individual; //[num_blocks + 1];

	//whether significantly different or not
	bool significant;
	bool* significant_individual; //[num_blocks + 1];

	gene();
	~gene();
}; 

vector <gene*> sortList_gene; //fragment list to be sorted
double* sortKey_gene; //[MAX_TESTGENE_NUM];
long sortList_gene_Num = 0;

double* mergeSort_Larray_gene; //[MAX_TESTGENE_NUM];
double* mergeSort_Rarray_gene; //[MAX_TESTGENE_NUM];
gene* *mergeSort_LorderedList_gene; //[MAX_TESTGENE_NUM];
gene* *mergeSort_RorderedList_gene; //[MAX_TESTGENE_NUM];



double q25 = 0.0; //25% point of the permuted d values
double q75 = 0.0; //75% point of the permuted d values
double pi0 = 0.0; //pi_0, estimate for the proportion of true null genes

//all kinds of counts
long permutationCnt = 0; //count of permutations
long allpermutedvalueCnt = 0; //count for all permuted values, which is geneCnt * permutationCnt
int permutation_limit_ratio = 1;
permutationtype permType = perm_full;

//long* stack;
double* falseGeneCntList;

double*** statistics_expression; //[MAX_TESTGENE_NUM][GROUP_NUM+1][SAMPLE_CNT_PER_GROUP+1];
double* statistics_expected_d; //[MAX_TESTGENE_NUM];
double* statistics_d; //[MAX_TESTGENE_NUM];
double* statistics_s; //[MAX_TESTGENE_NUM];

double** statistics_permuted_d; //[MAX_TESTGENE_NUM]; //all permuted value

int* individualList; //[num_blocks+1];

bool diff_expression_analysis_two_class(double*** expression, double *statistic_d, double *statistic_s, double *statistics_expected_d, double* statistics_permuted_d[], long testGeneNum, int groupNum, int individualNum, int *individualList, long sampleSize, string resultPath, string resultFileSuffix);
extern long choose(long n, long k);
extern long choose_non_recursive(long n, long k);
extern void quicksort(double *sortArray, long length);
extern double calculate_mean(double *dataarray, long length);
extern double calculate_percentile(double *dataarray, long length, double percentile_0to1);

const double MAX_FOLDCHANGE = 10000;
double* NORMALIZATION_FACTOR; //[GROUP_NUM * SAMPLE_CNT_PER_GROUP + 1];

double significance_cutoff = 0; //cutoff for significant differences



