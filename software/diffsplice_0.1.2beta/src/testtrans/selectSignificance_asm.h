/*    
 *    selectSignificance_asm.h		
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


//#define OUTPUT_PAIRWISE_JSD


/************************************************************************/
/* GENERAL SETTINGS                                                     */
/************************************************************************/

int GROUP_NUM = 0;
int INDIVIDUAL_NUM = 0; //number of individuals, useless now, replaced by num_blocks
int SAMPLE_CNT_PER_GROUP = 0; //number of samples per group
bool COUNT_INTRON_RETENTION = false; //whether to take intron retention into account
int TOTAL_SAMPLE_CNT = 0; //total number of samples

int num_blocks = 0; //number of permutation blocks (individuals)
int* TECHNICAL_REP_NUM = NULL; //number of technical replicates of each individual
int* TECHNICAL_REP_ARRAY_BASE = NULL; //array base of technical replicates of each individual


double THRESHOLD_MIN_ASM_COVERAGE = 0; //all groups must have mean expression larger than this threshold to make differential transcription
double false_discovery_rate = 0.01;
double thresh_JSD = 0.1;

/************************************************************************/
/* FIXED SETTINGS	                                                    */
/************************************************************************/

const long DEFAULT_TESTGENE_NUM = 100000;

//for the leukemia data, do not apply any filter on the group, the individual cutoff will take charge
const double THRESH_BAD_ESTIMATION = MAX_NUMBER; //an ASM will not be considered if more than half samples have estimation error ratio larger than this cutoff

//this one does not make sense and has been disabled
const double THRESHOLD_MIN_ASMPATH_COVERAGE = -1; //all groups in all time points must have mean expression larger than this threshold to make differential transcription


/************************************************************************/
/* DEFINITIONS & DECLARATIONS                                           */
/************************************************************************/


char alterSpliceCategory[][30] = {"uncataloged", "exon_skipping", "mutual_exclusive", "intron_retention", "alter_splice_site", "alter_trans_start/end", "alter_trans_start/end", "diff_expression"};


class ASM
{
public:
	//basic information
	string ASM_id;
	long testID;
	string chrNM;
	long rangeLow;
	long rangeHigh;
	int pathNum;
	int ASMcategory;

	double** individualExpression; //[GROUP_NUM+1][SAMPLE_CNT_PER_GROUP+1];
	double** error_ratio; //[GROUP_NUM+1][SAMPLE_CNT_PER_GROUP+1]; //relative ratio of the estimation error against the expression
	double* meanExpression; //[GROUP_NUM + 1];
	double minExpression;
	double minPathExpression;

	double*** individualPathExpression; //[GROUP_NUM+1][SAMPLE_CNT_PER_GROUP+1];
	double*** individualPathProportion; //[GROUP_NUM+1][SAMPLE_CNT_PER_GROUP+1];
	double** meanPathExpression; //[GROUP_NUM+1];
	double** meanPathProportion; //[GROUP_NUM+1];

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

	ASM(int pathnum);
	void calculate_meanPathExprProp();
	~ASM();
}; 



vector <ASM*> sortList_ASM; //fragment list to be sorted
double* sortKey_ASM;
long sortList_ASM_Num = 0;

double* mergeSort_Larray_ASM;
double* mergeSort_Rarray_ASM;
ASM* *mergeSort_LorderedList_ASM;
ASM* *mergeSort_RorderedList_ASM;



double q0 = 0.0; //25% point of the permuted d values
double q50 = 0.0; //75% point of the permuted d values
double pi0 = 0.0; //pi_0, estimate for the proportion of true null genes

//all kinds of counts
long permutationCnt = 0; //count of permutations
long allpermutedvalueCnt = 0; //count for all permuted values, which is geneCnt * permutationCnt
int permutation_limit_ratio = 1;
permutationtype permType = perm_full;

//long* stack;
double* falseGeneCntList;

double**** statistics_distribution;
long* statistics_dimensionList;
double* statistics_minPathExprList;
double* statistics_expected_d;
double* statistics_d;
double* statistics_s;

double** statistics_permuted_d; //all permuted value

int* individualList;

bool diff_expression_analysis_two_class(double**** pathDistribution, long *dimensionList, double *minPathExprList, double *statistic_d, double *statistic_s, double *statistics_expected_d, double* statistics_permuted_d[], long testASMNum, int groupNum, int individualNum, int *individualList, long sampleSize, string resultPath, string resultFileSuffix);
extern long choose(long n, long k);
extern long choose_non_recursive(long n, long k);
extern void quicksort(double *sortArray, long length);
extern double calculate_mean(double *dataarray, long length);
extern double calculate_JSD(double *array_P, double *array_Q, long length);
extern double calculate_percentile(double *dataarray, long length, double percentile_0to1);


double significance_cutoff = 0; //cutoff for significant differences


