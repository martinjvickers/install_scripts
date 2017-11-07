/*    
 *    GTree.h		
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
/* SETTINGS		                                                        */
/************************************************************************/

//for the input of splice junctions
const int min_junction_intron_length = 10; //consider junctions whose end_position - start_position - 1 < min_intron_length as small deletion and do not include them in splice graph


const long MIN_EXON_LENGTH = 5;
const long MIN_ALTER_SPLICE_SITE_LENGTH = 1; //minimum exonic length for alternative splice sites. need a maximum??
const long MAX_NOCOVERAGE_LENGTH = 1; //maximal length for an exon to have no coverage

//filter exon in gene if its coverage is less than 0.05 of mean gene coverage in more than 0.95 samples
const double coverageThreshold_GeneExon = 0; //percentage of mean gene coverage

const double MIN_COVERAGE_FILTER_FRAGMENT = 0; //minimum coverage for fragments when filtering standalone fragments.

const int MAXJUNCNUMINDEPPATH = 15; //maximum number of junctions when separating dependent paths. if #junctions > this number, do not enumerate paths

const int COVERAGE_CHANGE_WINDOW = 10;
const double COVERAGE_CHANGE_THRESH = 5.0;

int thresh_path_expressed = 1;

/************************************************************************/
/* CLASS DECLARATION	                                                */
/************************************************************************/

class GTvertex;
class fragment;

class JuncGraphVertex;

/************************************************************************/
/* TYPE DEFINITION	                                                    */
/************************************************************************/

enum fragment_type {frag_exon, frag_intron, frag_retained_intron, frag_junction, frag_insertion, frag_deletion, frag_SNP};
enum trans_direction {undetermined, sense, antisense};
enum JuncGraphVertexType {virStart, virEnd, normal};

/************************************************************************/
/* CLASS DEFINITION - FRAGMENT                                          */
/************************************************************************/
class alter_junction
{
public:
	fragment *juncInfo;
//	int category; //0 for the merged target junction, -1 for 5', 1 for 3', 2 for exon
//	double proportion;
	alter_junction *next;

	alter_junction();
};

class fragment
{
public:
	string frag_name;
	unsigned long ID; //unique ID throughout the pipeline
	string chromosome_start;
	string chromosome_end;
	trans_direction transDirection;
	long start;
	long end;
	fragment_type type; //fragment type
	int altersite; //if the fragment is an exon, then this field indicates whether this fragment is an alternative splice site: 0 for exon, 1 for donor site, 2 for acceptor site

	double *support; //support from the fragment file

	alter_junction* alter;
	long alterFragCnt; //fragment count of alternative splice sites, note it's the number of fragments
//	alter_junction* fivePalterTail; //tail pointer of 5' alternative junction
	long *coverage;
	long start_real; //real start position in case the position is changed due to alternative splice sites
	long end_real; //real end position

	fragment();
	fragment* clone(); //make a clone, used for exon, so alter_junction not included
	~fragment();
};

class spliceSite
{
public:
	string chromosome;
	long position;
	bool directionOut; //true for out-going junction, false for in-coming junction

	spliceSite();
};


/************************************************************************/
/* CLASS DEFINITION - GenomeTree                                        */
/************************************************************************/
class rangeJunction
{
	//fragment within a range
public:
	fragment *junc;
	rangeJunction *next;

	rangeJunction();
};

class RangeJunctionList
{
	//fragment list within a range
public:
	long rangeLow;
	long rangeHigh;
	trans_direction transDirection;
	rangeJunction *list;
	rangeJunction *listtail;
	RangeJunctionList *nextList;

	RangeJunctionList();
	RangeJunctionList* clone();
	~RangeJunctionList();
};


class GTedge
{
	//edge of Genome Tree
public:
	GTvertex *linkedVertex;
	GTedge *next;

	GTedge();
};

class alternative_path
{
public:
	long path_start; //site that starts diverging
	long path_end; //site that ends diverging
	long whole_path_start; //the whole path including the starting exon
	long whole_path_end; //the whole path including the ending exon
	trans_direction transDirection;
	double *support; //support at this vertex
	double *proportion; //proportion at this vertex
	long junctionNum;
	long exonNum;
	GTvertex *pathVertex;
	alternative_path *next;

	alternative_path();
	~alternative_path();
};

class GTvertex
{
	//vertex of Genome Tree
public:
	long ID;
	long level;
	long rangeLow;
	long rangeHigh;
	GTedge *child;
	int childType; //0 for no child (i.e. leaf node), 1 for independent region, 2 for independent path, 3 for dependent path
	long childNum;
	RangeJunctionList *junctionInRange; //junctions within the vertex's range
	long junctionNum;
	long exonNum;
	GTvertex *prevSibling; //the sibling on the left (5') side
	GTvertex *nextSibling; //the sibling on the right (3') side

	GTedge *alterSpliceSite; //alternative splice sites involved in this vertex

	//for difference analysis
	double *support; //support at this vertex
	double *proportion; //proportion at this vertex (for genes, this array stores mean 25%-75% coverage)
	double *MSE_estimation; //mean squared error of estimated transcript abundance
	double *min_path_support;
	double *obs_support; //observed support directly calculated from raw support (with no estimation)

	double total_inflow;
	double total_outflow;

	bool estimated; //if true, the vertex has been estimated and will be treated as a whole
	fragment *representative;
	long estimate_exonNum; //number of blocks under estimation, equivalent to number of exons previously

	int ASMcategory;
	double ASMsupport_group1;
	double ASMsupport_group2;

	alternative_path* major_alter_paths;
	long major_alter_paths_num;

	GTvertex();
	~GTvertex();
};

class GenomeTree
{
public:
	GTvertex* root;

	GenomeTree();
	~GenomeTree();
};




/************************************************************************/
/* CLASS DEFINITION - JUNCTION GRAPH		                            */
/************************************************************************/

class JuncGraphEdge
{
	//edge of fragment graph
public:
	JuncGraphVertex *linkedVertex;
	JuncGraphEdge *next;

	JuncGraphEdge();
};

class JuncGraphVertex
{
	//vertex of fragment graph
public:
	rangeJunction* corresJunc;
	JuncGraphEdge* edges;
	JuncGraphVertex* next;
	bool traversed;
	bool hasInEdge; //if hasInEdge is false, consider it as a start vertex when enumerating transcripts
	bool hasOutEdge;

	JuncGraphVertexType vertexType;

	JuncGraphVertex();
	~JuncGraphVertex();
};

class JuncGraphPath
{
	//path of fragment graph
public:
	RangeJunctionList *pathJuncList;
	JuncGraphVertex *arrivedVertex;

	JuncGraphPath();
	~JuncGraphPath();
};

class JuncGraph
{
	//fragment graph
public:
	JuncGraphVertex *vertices;
	
	JuncGraph();
	~JuncGraph();
};


//class for filtering
class deletedSites
{
public:
	long sites;
	deletedSites *next;

	deletedSites();
};



/************************************************************************/
/* FUNCTION DECLARATION	                                                */
/************************************************************************/

void mergeSort_JunctionSort(long sortList_size);
RangeJunctionList* separateIndepRegion(RangeJunctionList* origList, bool &separable);
RangeJunctionList* separateIndepPath(RangeJunctionList* origList, bool &separable);
double abundance_estimation(GTvertex *rootVertex, int index_tissue, double &MSE);
bool countGTree(GTvertex *rootVertex);
double JSD_significance_level(double JSDvalue, long dimension, double N);
void calculate_ASM_group_meanExpression(GTvertex *rootVertex);
//bool alterSpliceReliability(GTvertex *curVertex);
void output_ASMpath_gtf(GTvertex *targetVertex, RangeJunctionList *pathList);
void filterFalseFragments(bool	);
void constructGTree(RangeJunctionList* origList);





/************************************************************************/
/* GLOBAL VARIABLES (NO NEED TO CHANGE)                                 */
/************************************************************************/

bool DOTRIMMING = true;
bool DORETAIN = true;

bool CLEAN_SPLICE_DIRECTION = false; //need to correct transcription strand for splice junctions (two splice junctions have same location but different strands)

string inputPath;
string resultPath;
string resultNamePrefix;

long GENEcount = 0;
long ASMcount = 0;
unsigned long fragmentID_Cnt = 0;


const long MAX_CHR_LENGTH = 1000000000; //maximal chromosome length
const long DEFAULT_JUNCTION_NUM = 1000000;

long CHROMOSOME_START = 0;
long CHROMOSOME_END   = MAX_CHR_LENGTH; 

const int unknown = 0;
const int exon_skipping = 1;
const int mutual_exclusive = 2;
const int intron_retention = 3;
const int alter_splice_site = 4;
const int alter_start = 5;
const int alter_end = 6;
const int diff_expression = 7;

long exon_skipping_cnt = 0;
long mutual_exclusive_cnt = 0;
long intron_retention_cnt = 0;


//the g-tree
GenomeTree *gTree;

//stack for DFS on the g-tree
vector <GTvertex*> stack_GTvertex;

//graph for the search of paths
JuncGraph *junctionGraph;

//for sorting junctions
vector <rangeJunction*> sortList_Junction; //for minimal change: check size every time putting in new elements, never clear; still use sortList_Junction_Num to control/identify # of elements in the vector
vector <long> sortKey_Junction;
long sortList_Junction_Num = 0;

vector <long> mergeSort_Larray;
vector <long> mergeSort_Rarray;
vector <rangeJunction*> mergeSort_LorderedList;
vector <rangeJunction*> mergeSort_RorderedList;


//for sorting splice sites
long sortList_spliceSite_Num;
vector <spliceSite*> sortList_spliceSite;
vector <long> sortKey_spliceSite;
vector <spliceSite*> mergeSort_LorderedList_spliceSite;
vector <spliceSite*> mergeSort_RorderedList_spliceSite;

//junction graph path queue for computing all paths
queue <JuncGraphPath*> juncPathQueue;

//main output
ofstream outfile_gtree;
ofstream outfile_not_enumerated;
ofstream outfile_stats_expr;
ofstream outfile_stat_asm;
ofstream outfile_gtf_asm_path;
ofstream outfile_gtf_splice_graph;
ofstream outfile_junction_all;
ofstream outfile_junction_filtered;
ofstream outfile_asm_composition; //output composition of ASMs

//For normalization
vector <long> DatasetReadCount_total;
double DatasetReadCount_normalization = 10000;
vector <double> DatasetNormalizationRatio;

//For collecting d statistics
vector <GTvertex*> vertexListForStatistics;




