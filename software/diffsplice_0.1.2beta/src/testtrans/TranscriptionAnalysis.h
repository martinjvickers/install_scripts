/*    
 *    TranscriptionAnalysis.h		
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

//Array index starts from 1

#include "Configuration.h"

//permutation using block design: keep individuals balanced in each permutation. otherwise, a random permutation will be used (to be implemented later)
#define PERMUTATION_WITHIN_INDIVIDUAL


class testASM
{
public:
	long testID;
	
	int dimension;
	double*** individualDistribution;
	double** meanDistribution;
	double minPathExpr;

	double statistic_d;
	double statistic_s;
	
	testASM(long ID, int dim, int groupNum, int groupSize);
	void calculate_meanDistribution(int groupNum, int groupSize);
	~testASM();
};


testASM** sortList_testASM; //[MAX_TESTGENE_NUM]; //fragment list to be sorted
double* sortKey_testASM; //[MAX_TESTGENE_NUM];
long sortList_testASM_Num = 0;

double* mergeSort_Larray_testASM; //[MAX_TESTGENE_NUM];
double* mergeSort_Rarray_testASM; //[MAX_TESTGENE_NUM];
testASM** mergeSort_LorderedList_testASM; //[MAX_TESTGENE_NUM];
testASM** mergeSort_RorderedList_testASM; //[MAX_TESTGENE_NUM];

long* quicksortStack; //[MAX_TESTGENE_NUM];

int TEST_GROUP_NUM = 0;
int TEST_GROUP_SIZE = 0;




