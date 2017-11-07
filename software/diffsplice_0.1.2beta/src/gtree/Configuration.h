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

/************************************************************************/
/* GENERAL SETTINGS                                                     */
/************************************************************************/

#define UNIX

#ifdef UNIX

#include <fstream>
#include <cstring>
#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <ctime>
#include <cstdio>
#include <sstream>
#include <vector>
#include <queue>

#else

#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <math.h>
#include <stdlib.h> 
#include <time.h>
#include <sstream>
#include <vector>
#include <queue>

#endif


using namespace std;


#define DO_ESTIMATION
#define NORMALIZE_COVERAGE
#define TRANSCRIPTION_DIRECTION
#define FILTER_JUNCTION
#define FILTER_FRAGMENTS
#define COVERAGE


//#define JUNCTIONONLY
//#define DO_MERGEALTERSITES
//#define COUTSTEPS

int SUPPORT_VECTOR_SIZE = 0;
double thresh_junction_filter_max_read_support = 0;
double thresh_junction_filter_mean_read_support = 0;
double thresh_junction_filter_num_samples_presence = 0;

bool COUNT_MAJOR_PATHS_ONLY = false;
double coverageThreshold_exon = 0; //throw-away := <= coverageThreshold * SUPPORT_VECTOR_SIZE
double coverageThreshold_intron = 1; //throw-away := <= coverageThreshold * SUPPORT_VECTOR_SIZE

#endif

