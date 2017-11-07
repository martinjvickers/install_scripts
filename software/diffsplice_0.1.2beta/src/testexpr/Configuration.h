/*    
 *    ExpressionAnalysis.h		
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

#ifndef CONFIGURATION
#define CONFIGURATION

#define UNIX

#ifdef UNIX

#include <fstream>
#include <cstring>
#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <ctime>
#include <cstdio>
#include <vector>

#else

#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <math.h>
#include <stdlib.h> 
#include <time.h>
#include <vector>

#endif

using namespace std;

//#define OUTPUT_PAIRWISE_FOLDCHANGE

extern int* TECHNICAL_REP_NUM;
extern int* TECHNICAL_REP_ARRAY_BASE;

const double MAX_NUMBER = 2000000000;

//permutation count
extern long permutationCnt; //count of permutations
extern long allpermutedvalueCnt; //count for all permuted values, which is geneCnt * permutationCnt
extern int permutation_limit_ratio;
enum permutationtype {perm_full, perm_reducedfull, perm_random};
extern permutationtype permType;

//maximum permutation count
const long MAX_PERMUTATION_NUMBER = 10000; //MAX_NUMBER if not limited



#endif

