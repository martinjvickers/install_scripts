/*    
 *    diffsplice.h		
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


#define  UNIX

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

#endif


using namespace std;

bool SEP_SAM = false;
bool SEP_FRAG = false;

bool SEP_TREE = false;
bool SELECT_SIGNIFICANCE = false;

//bool SUMMARIZE = false;
//bool DELETE_OLD_RESULT = true;


vector <string> chrName;

void inputChrName(string filename);



class sample
{
public:
	string name;
	int groupID;
	int bioRepID;
	int techRepID;
	string datafilename;
	sample *next;

	sample(string, int, int, int);
	~sample();
};

class bioReplicate
{
public:
	string name;
	int bioRepID;
	int techRepNum;
	sample *techReps;
	sample *techRepsTail;
	bioReplicate *next;

	bioReplicate(string, int);
	~bioReplicate();
};

class sampleGroup
{
public:
	string name;
	int groupID;
	int bioRepNum;
	int sampleNum;
	bioReplicate *bioReps;
	bioReplicate *bioRepsTail;
	sampleGroup *next;

	sampleGroup(string, int);
	~sampleGroup();
};


int groupNum = 0;
int totalSampleNum = 0;
sampleGroup* allGroups = NULL;
sampleGroup* allGroupsTail = NULL;

double thresh_junction_filter_max_read_support = 0;
double thresh_junction_filter_mean_read_support = 0;
double thresh_junction_filter_num_samples_presence = 0;
bool COUNT_MAJOR_PATHS_ONLY = false;
double coverageThreshold_exon = 0; //throw-away := <= coverageThreshold * SUPPORT_VECTOR_SIZE
double coverageThreshold_intron = 1; //throw-away := <= coverageThreshold * SUPPORT_VECTOR_SIZE

bool config_bal_design = false; //balanced design for permutation test
double config_false_discovery_rate = 0.01;
double config_thresh_foldchange_up = 2.;
double config_thresh_foldchange_down = 0.5;
double config_thresh_sqrtJSD = 0.1;


void parse_inputfiles(char *filename); //extract the data files from the given file. 
void parse_settings(char *filename); //extract the settings from the given file. 
string itostr(int); //convert integer to string
void write_configfile(string target_path);
void check_sample_count();


//some other settings
double THRESHOLD_MIN_EXPR_COVERAGE = 10;//at least one group must have mean expression larger than this threshold to make differential expression;
double THRESHOLD_MIN_ASM_COVERAGE = 1; //all groups must have mean expression larger than this threshold to make differential transcription

