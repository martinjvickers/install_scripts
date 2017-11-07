/*    
 *    diffsplice.cpp		
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


#include "diffsplice.h"

sample::sample(string samplename, int groupid, int biorepid, int techrepid)
{
	name = samplename;
	groupID = groupid;
	bioRepID = biorepid;
	techRepID = techrepid;
	datafilename = "";
	next = NULL;
}

sample::~sample(){}

bioReplicate::bioReplicate(string repname, int ID)
{
	name = repname;
	bioRepID = ID;
	techRepNum = 0;
	techReps = NULL;
	techRepsTail = NULL;
	next = NULL;
}

bioReplicate::~bioReplicate()
{
	sample *delSample = techReps;
	while (delSample != NULL)
	{
		techReps = techReps->next;
		delete delSample;
	}
}

sampleGroup::sampleGroup(string groupname, int ID)
{
	name = groupname;
	groupID = ID;
	bioRepNum = 0;
	sampleNum = 0;
	bioReps = NULL;
	bioRepsTail = NULL;
	next = NULL;
}

sampleGroup::~sampleGroup()
{
	bioReplicate *delBioRep = bioReps;
	while (delBioRep != NULL)
	{
		bioReps = bioReps->next;
		delete delBioRep;
	}
}



int main(int argc, char* argv[])
{
	int timepointLoopCnt, individualLoopCnt, duplicateLoopCnt, chrLoopCnt, i;
	long permutationLoopCnt, permutationCnt, tmpLong;
	string target_path, fname, comd, runName = "diffsplice", filepath, resultpath, filename;
	sampleGroup *curGroup;
	bioReplicate *curBioRep;
	sample *curSample;
	time_t t;
	string cur_time;


	time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
	cout << "[" << cur_time << "]\tDiffSplice 0.1.1" << endl;

#ifdef UNIX
	if (argc < 4)
	{
		cout << argv[0] << "\t[options]\tsettings.cfg\tdatafile.cfg\tresult_path" << endl;
		return 1;
	}
	else
	{
		int cmdLoop;
		for (cmdLoop = 1; cmdLoop < argc && argv[cmdLoop][0] == '-'; ++cmdLoop)
		{
			for (int tmpLoop = 1; argv[cmdLoop][tmpLoop] != '\0'; ++tmpLoop)
			{
				if (argv[cmdLoop][tmpLoop] == 'p')
				{
					SEP_SAM = true;
					SEP_FRAG = true;
				}
				else if (argv[cmdLoop][tmpLoop] == 's')
					SEP_TREE = true;
				else if (argv[cmdLoop][tmpLoop] == 'x')
					SELECT_SIGNIFICANCE = true;
				else
				{
					cout << "Error: unrecognized option in command line." << endl;
					return 1;
				}
			}
		}

		if (cmdLoop == 1)
		{
			//all steps as default 
			SEP_SAM = true;
			SEP_FRAG = true;
			SEP_TREE = true;
			SELECT_SIGNIFICANCE = true;
		}
		if (cmdLoop + 3 > argc)
		{
			cout << argv[0] << "\t[options]\tsettings.cfg\tdatafile.cfg\tresult_path" << endl;
			return 1;
		}

		parse_settings(argv[cmdLoop++]);
		parse_inputfiles(argv[cmdLoop++]);

		target_path = argv[cmdLoop];
		if (target_path[target_path.length()-1] != '/')
			target_path += "/";
	}
#else
	parse_settings("settings.cfg");
	parse_inputfiles("datafile.cfg");
	target_path = "";
#endif

	comd = "mkdir -p " + target_path + "data"; system(comd.c_str());

	check_sample_count();
	write_configfile(target_path + "data/");


	time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
	cout << "[" << cur_time << "]\tData set: " << groupNum << " groups with " << allGroups->sampleNum << " samples in each group" << endl;

	
	if (SEP_SAM == true)
	{
		time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
		cout << "[" << cur_time << "]\tBegin to parse data files (in SAM format)" << endl;

		comd = "rm -f -r " + target_path + "data/sam"; system(comd.c_str());
		comd = "mkdir -p " + target_path + "data/sam"; system(comd.c_str());

		curGroup = allGroups;
		while (curGroup != NULL)
		{
			comd = "mkdir -p " + target_path + "data/sam/group" + itostr(curGroup->groupID); system(comd.c_str());

			curBioRep = curGroup->bioReps;
			while (curBioRep != NULL)
			{
				curSample = curBioRep->techReps;
				while (curSample != NULL)
				{
					comd = "mkdir -p " + target_path + "data/sam/group" + itostr(curGroup->groupID) + "/ID" + itostr(curBioRep->bioRepID) + "_rep" + itostr(curSample->techRepID); system(comd.c_str());
					comd = "sepSAM " + curSample->datafilename + " " + target_path + "data/sam/group" + itostr(curGroup->groupID) + "/ID" + itostr(curBioRep->bioRepID) + "_rep" + itostr(curSample->techRepID) + "/"; system(comd.c_str());

					curSample = curSample->next;
				}

				curBioRep = curBioRep->next;
			}

			curGroup = curGroup->next;
		}
	}

	filepath = target_path + "data/sam/group" + itostr(allGroups->bioReps->techReps->groupID) + "/ID" + itostr(allGroups->bioReps->techReps->bioRepID) + "_rep" + itostr(allGroups->bioReps->techReps->techRepID) + "/";
	inputChrName(filepath);

// 	time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
// 	cout << "[" << cur_time << "]\tTotal # of chromosomes = " << chrName.size() << endl;

	if (SEP_FRAG == true)
	{
		time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
		cout << "[" << cur_time << "]\tBegin to prepare input files for DiffSplice" << endl;

		comd = "rm -f -r " + target_path + "data/frag"; system(comd.c_str());
		comd = "mkdir -p " + target_path + "data/frag"; system(comd.c_str());

		for (chrLoopCnt = 0; chrLoopCnt < chrName.size(); ++chrLoopCnt)
		{
			comd = "mkdir -p " + target_path + "data/frag/" + chrName[chrLoopCnt];
			system(comd.c_str());
		}

		int sampleCnt = 0;

		curGroup = allGroups;
		while (curGroup != NULL)
		{
			curBioRep = curGroup->bioReps;
			while (curBioRep != NULL)
			{
				curSample = curBioRep->techReps;
				while (curSample != NULL)
				{
					++sampleCnt;

					filepath = target_path + "data/sam/group" + itostr(curGroup->groupID) + "/ID" + itostr(curBioRep->bioRepID) + "_rep" + itostr(curSample->techRepID) + "/";
					for (chrLoopCnt = 0; chrLoopCnt < chrName.size(); ++chrLoopCnt)
					{
						comd = "fragment " + filepath + " " + chrName[chrLoopCnt] + " " + itostr(sampleCnt) + " noDB " + target_path + "data/frag/" + chrName[chrLoopCnt] + "/"; 
						system(comd.c_str());
					}

					curSample = curSample->next;
				}

				curBioRep = curBioRep->next;
			}

			curGroup = curGroup->next;
		}

		for (chrLoopCnt = 0; chrLoopCnt < chrName.size(); ++chrLoopCnt)
		{
			comd = "sort -n +1 -4 " + target_path + "data/frag/" + chrName[chrLoopCnt] + "/allJunction.txt > " + target_path + "data/frag/" + chrName[chrLoopCnt] + "/junction_all.txt";
			system(comd.c_str());	
			comd = "rm -f " + target_path + "data/frag/" + chrName[chrLoopCnt] + "/allJunction.txt";
			system(comd.c_str());
			comd = "sort -n +1 -4 " + target_path + "data/frag/" + chrName[chrLoopCnt] + "/allExon.txt > " + target_path + "data/frag/" + chrName[chrLoopCnt] + "/exon_all.txt";
			system(comd.c_str());	
			comd = "rm -f " + target_path + "data/frag/" + chrName[chrLoopCnt] + "/allExon.txt";
			system(comd.c_str());
		}
	}

// 	if (DELETE_OLD_RESULT == true)
// 	{
// 		sprintf(comd, "rm -r %sresult", path.c_str());
// 		system(comd);
// 	}
	
	if (SEP_TREE == true)
	{
		time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
		cout << "[" << cur_time << "]\tBegin to run DiffSplice" << endl;


		filepath = target_path + "data/frag/";
		resultpath = target_path + "result/";
		comd = "rm -f -r " + target_path + "result"; system(comd.c_str());
		comd = "mkdir -p " + target_path + "result"; system(comd.c_str());
		comd = "mkdir -p " + target_path + "result/stat"; system(comd.c_str());
		comd = "mkdir -p " + target_path + "result/detail"; system(comd.c_str());
		comd = "mkdir -p " + target_path + "result/stat/expr"; system(comd.c_str());
		comd = "mkdir -p " + target_path + "result/stat/asm"; system(comd.c_str());


		ofstream asm_path_file, non_asm_file, splice_all_file, splice_filter_file;
		filename = target_path + "result/asm_path.gtf";
		asm_path_file.open(filename.c_str());
		asm_path_file << "browser full " << runName << "_ASM_path\ntrack name=\"" << runName << "_ASM_path\" description=\"" << runName << "_ASM_path\" visibility=2 useScore=1 itemRgb=on" << endl;
		asm_path_file.close();

		filename = target_path + "result/splice_graph.gtf";
		non_asm_file.open(filename.c_str());
		non_asm_file << "browser dense " << runName << "_splice_graph\ntrack name=\"" << runName << "_splice_graph\" description=\"" << runName << "_splice_graph\" visibility=1 useScore=1 itemRgb=on" << endl;
		non_asm_file.close();

		filename = target_path + "result/stat/splice_all.bed";
		splice_all_file.open(filename.c_str());
		splice_all_file << "browser full " << runName << "_junction_all\ntrack name=\"" << runName << "_junction_all\" description=\"" << runName << "_junction_all\" visibility=2 useScore=1 itemRgb=on" << endl;
		splice_all_file.close();

		filename = target_path + "result/stat/splice_filtered.bed";
		splice_filter_file.open(filename.c_str());
		splice_filter_file << "browser full " << runName << "_junction_filtered\ntrack name=\"" << runName << "_junction_filtered\" description=\"" << runName << "_junction_filtered\" visibility=2 useScore=1 itemRgb=on" << endl;
		splice_filter_file.close();


		for (i = 0; i < chrName.size(); i++)
		{
			time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
			cout << "[" << cur_time << "]\t" << flush;
			comd = "tree " + filepath + chrName[i] + "/ " + resultpath + " " + chrName[i] + " " + target_path + "data/config_gtree";
			system(comd.c_str());
			cout << endl;
		}
	}

	if (SELECT_SIGNIFICANCE == true)
	{
		resultpath = target_path + "result/";

		time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
		cout << "[" << cur_time << "]\tBegin the differential expression analysis on genomic loci" << endl;

		comd = "diff_expr_analysis " + resultpath + " " + target_path + "data/config_testexpr"; 
		system(comd.c_str());

		time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
		cout << "[" << cur_time << "]\tBegin the differential transcription analysis on ASMs" << endl;

		comd = "diff_asm_analysis " + resultpath + " " + target_path + "data/config_testtrans"; 
		system(comd.c_str());
	}

	time(&t); cur_time = ctime(&t); cur_time.erase(cur_time.size() - 1);
	cout << "[" << cur_time << "]\tFinish DiffSplice run" << endl;

	
// 	if (SUMMARIZE == false)
// 		return 0;
// 
// 	sprintf(comd, "rm -r ./Analysis/%s", runName);
// 	system(comd);
// 	sprintf(comd, "mkdir ./Analysis/plot");
// 	system(comd);
// 
// 	sprintf(filename, "./Analysis/plot/comparison/gene_summary.txt");
// 	ofstream gene_summary_file;
// 	gene_summary_file.open(filename);
// 	gene_summary_file << "chrNm\tstart\tend\t#db_trans\t#db_ASM\tdb_ASM_size\tdb_ASM_size_overlap\tdb_ASM_level\t#data_ASM\tdata_ASM_size\tdata_ASM_size_overlap\tdata_ASM_level\n";
// 	gene_summary_file.close();
// 
// 	sprintf(filename, "./Analysis/plot/ASM_summary.txt");
// 	ofstream ASM_summary_file;
// 	ASM_summary_file.open(filename);
// 	ASM_summary_file << "chrNm\tstart\tend\tcategory\tannotated\t#path_annotation\tlevel_annotation\tinData\t#path_data\tlevel_data\tdata_expr_G1\tdata_expr_G2\n";
// 	ASM_summary_file.close();
// 
// 	sprintf(filename, "./Analysis/plot/geneList_validation.txt");
// 	ofstream validationList;
// 	validationList.open(filename);
// 	validationList << "Entrez_ID\tGene\tGene_ID\tmiR_targets\tChrNm\tStrand\tStart\tEnd\tSig_diff_expr\tStat_expr (|d - Ed|)\tFold_change\tSig_diff_trans\tStat_trans (|d - Ed|)\tASM category\tMean JSD\tCuffdiff_ref\tJSD_ref\tpvalue_ref\tCuffdiff_noref\tJSD_noref\tpvalue_noref\tEstTransNum\t\tmeanExpr_1\tmeanExpr_2\tdb_isoform_num\tdata_ASM_num\n";
// 	validationList.close();
// 
// 
// 	sprintf(comd, "./summarize_%d %s ./data/UCSC_KnownGene_Dec2009.gff ./result/%s/ ./Analysis/ ./ChromosomeName.txt ./result/annotationASM/ ", runID, 
// 		runName, runName);
// 	system(comd);
// 
// 	sprintf(comd, "./ASMcompare_%d ./result/%s/ ./Analysis/plot/", runID, runName);
// 	system(comd);
// 
// 
// 	cout << "done" << endl;
// 
// 	sprintf(comd, "mv ./Analysis/plot ./Analysis/%s", runName);
// 	system(comd);

	return 0;
}



void inputChrName(string filename)
{
	ifstream chrFile;
	string cur_chr, info;
	filename += "ChromosomeName.txt";
	chrFile.open(filename.c_str());
	
	while (chrFile >> cur_chr)
	{
		getline(chrFile, info);

		if (cur_chr.empty() == false && cur_chr[0] != '#' && cur_chr.compare("*") != 0)
			chrName.push_back(cur_chr);
	}
	chrFile.close();

	return;
}


void parse_inputfiles(char *filename)
{
	ifstream infile;
	infile.open(filename);

	string info;
	sampleGroup *curGroup;
	bioReplicate *curBioRep;
	sample *curSample;

	if (infile.is_open())
	{
		while (infile >> info)
		{
			if (info[0] == '#' || info[0] == '\0')
			{
				//comment line
				getline(infile, info);
			}
			else
			{
				//data line

				//find the group
				curGroup = allGroups;
				while (curGroup != NULL)
				{
					if (curGroup->name.compare(info) == 0)
						break;
					curGroup = curGroup->next;
				}

				if (curGroup == NULL)
				{
					curGroup = new sampleGroup(info, ++groupNum);
					if (allGroupsTail == NULL)
					{
						allGroups = curGroup;
						allGroupsTail = curGroup;
					}
					else
					{
						allGroupsTail->next = curGroup;
						allGroupsTail = curGroup;
					}
				}


				//find the bioReplicate
				infile >> info;
				curBioRep = curGroup->bioReps;
				while (curBioRep != NULL)
				{
					if (curBioRep->name.compare(info) == 0)
						break;
					curBioRep = curBioRep->next;
				}

				if (curBioRep == NULL)
				{
					curBioRep = new bioReplicate(info, ++(curGroup->bioRepNum));
					if (curGroup->bioRepsTail == NULL)
					{
						curGroup->bioReps = curBioRep;
						curGroup->bioRepsTail = curBioRep;
					}
					else
					{
						curGroup->bioRepsTail->next = curBioRep;
						curGroup->bioRepsTail = curBioRep;
					}
				}

				//find the sample
				infile >> info;
				curSample = curBioRep->techReps;
				while (curSample != NULL)
				{
					if (curSample->name.compare(info) == 0)
					{
						cout << "Warning: duplicated sample file." << endl;
						break;
					}
					curSample = curSample->next;
				}

				if (curSample == NULL)
				{
					++totalSampleNum;
					++curGroup->sampleNum;
					curSample = new sample(info, curGroup->groupID, curBioRep->bioRepID, ++(curBioRep->techRepNum));
					if (curBioRep->techRepsTail == NULL)
					{
						curBioRep->techReps = curSample;
						curBioRep->techRepsTail = curSample;
					}
					else
					{
						curBioRep->techRepsTail->next = curSample;
						curBioRep->techRepsTail = curSample;
					}
				}

				infile >> curSample->datafilename;
				getline(infile, info);
			}
		}
	}
	else
	{
		cout << "fail to open the data file" << endl;
		exit(1);
	}

	if (allGroups == NULL)
	{
		cout << "fail to import datafile.cfg" << endl;
		exit(1);
	}

	if (config_bal_design == true)
	{
		//check and correct the order of the biological replicates
		sampleGroup *firstgroup = allGroups;
		bioReplicate *allReps, *prevRep, *curRep;

		curGroup = firstgroup->next;
		while (curGroup != NULL)
		{
			allReps = curGroup->bioReps;
			curGroup->bioReps = NULL;
			curGroup->bioRepsTail = NULL;

			curBioRep = firstgroup->bioReps;
			while (curBioRep != NULL)
			{
				curRep = allReps;
				prevRep = NULL;
				while (curRep != NULL && curRep->name.compare(curBioRep->name) != 0)
				{
					prevRep = curRep;
					curRep = curRep->next;
				}

				if (curRep == NULL)
				{
					cout << "Error: missing biological replicate for balanced permutation test in group " << curGroup->name << "." << endl;
					exit(1);
				}
				else
				{
					if (curRep->techRepNum != curBioRep->techRepNum)
					{
						cout << "Error: please have the same number of technical replicates for biological replicate " << curRep->name << " in every group.\n";
						exit(1);
					}

					if (prevRep == NULL)
						allReps = curRep->next;
					else
						prevRep->next = curRep->next;

					if (curGroup->bioReps == NULL)
					{
						curGroup->bioReps = curRep;
						curGroup->bioRepsTail = curRep;
					}
					else
					{
						curGroup->bioRepsTail->next = curRep;
						curGroup->bioRepsTail = curRep;
					}
				}

				curBioRep = curBioRep->next;
			}

			if (allReps != NULL)
			{
				cout << "Error: extra biological replicate for balanced permutation test in group " << curGroup->name << "." << endl;
				exit(1);
			}

			curGroup = curGroup->next;
		}
	}

	return;
}


void parse_settings(char *filename)
{
	ifstream infile;
	infile.open(filename);

	string info;

	if (infile.is_open())
	{
		while (infile >> info)
		{
			if (info[0] == '#' || info[0] == '\0')
			{
				//comment line
				getline(infile, info);
			}
			else
			{
				//data line
				if (info.compare("balanced_design_for_permutation_test") == 0)
				{
					infile >> info;
					if (info.compare("yes") == 0)
						config_bal_design = true;
					else
						config_bal_design = false;
				}
				else if (info.compare("ignore_minor_alternative_splicing_variants") == 0)
				{
					infile >> info;
					if (info.compare("yes") == 0)
						COUNT_MAJOR_PATHS_ONLY = true;
					else
						COUNT_MAJOR_PATHS_ONLY = false;
				}
				else if (info.compare("thresh_junction_filter_max_read_support") == 0)
				{
					infile >> thresh_junction_filter_max_read_support;
				}
				else if (info.compare("thresh_junction_filter_mean_read_support") == 0)
				{
					infile >> thresh_junction_filter_mean_read_support;
				}
				else if (info.compare("thresh_junction_filter_num_samples_presence") == 0)
				{
					infile >> thresh_junction_filter_num_samples_presence;
				}
				else if (info.compare("thresh_average_read_coverage_exon") == 0)
				{
					infile >> coverageThreshold_exon;
				}
				else if (info.compare("thresh_average_read_coverage_intron") == 0)
				{
					infile >> coverageThreshold_intron;
				}
				else if (info.compare("false_discovery_rate") == 0)
				{
					infile >> config_false_discovery_rate;
				}
				else if (info.compare("thresh_foldchange_up") == 0)
				{
					infile >> config_thresh_foldchange_up;
				}
				else if (info.compare("thresh_foldchange_down") == 0)
				{
					infile >> config_thresh_foldchange_down;
				}
				else if (info.compare("thresh_sqrtJSD") == 0)
				{
					infile >> config_thresh_sqrtJSD;
				}
				else
				{
					cout << "Warning: unrecognized option in settings.cfg - " << info << endl;
				}
				
				getline(infile, info);
			}
		}
	}
	else
	{
		cout << "Fail to open the setting file. Will try with default settings." << endl;
	}

	return;
}

//check whether the samples provided are good to use
//now only allow two groups, and require same sample size for each group
void check_sample_count()
{
	if (groupNum != 2)
	{
		cout << "Please check datafile.cfg and make sure 2 sample groups are listed. DiffSplice now can only perform two-group comparison." << endl;
		exit(1);
	}
	if (allGroups->sampleNum != allGroupsTail->sampleNum)
	{
		cout << "Please check datafile.cfg and make sure 2 groups have the same number of samples." << endl;
		exit(1);
	}

	return;
}


string itostr(int t)
{
	ostringstream oss;
	oss << t;
	return oss.str();
}


void write_configfile_gtree(string target_path)
{
	ofstream conf_file;
	string filename;

	filename = target_path + "config_gtree";
	conf_file.open(filename.c_str());

	conf_file << "SUPPORT_VECTOR_SIZE\t" << totalSampleNum << endl; //total number of samples, for all groups
	conf_file << "thresh_junction_filter_max_read_support\t" << thresh_junction_filter_max_read_support << endl;
	conf_file << "thresh_junction_filter_mean_read_support\t" << thresh_junction_filter_mean_read_support << endl;
	conf_file << "thresh_junction_filter_num_samples_presence\t" << thresh_junction_filter_num_samples_presence << endl;


	conf_file << "COUNT_MAJOR_PATHS_ONLY\t";
	if (COUNT_MAJOR_PATHS_ONLY == true)
		conf_file << "yes" << endl;
	else
		conf_file << "no" << endl;
	conf_file << "coverageThreshold_exon\t" << coverageThreshold_exon << endl;
	conf_file << "coverageThreshold_intron\t" << coverageThreshold_intron << endl;

	conf_file.close();

	return;
}

void write_configfile_testexpr(string target_path)
{
	ofstream conf_file;
	string filename;

	filename = target_path + "config_testexpr";
	conf_file.open(filename.c_str());

	conf_file << "GROUP_NUM " << groupNum << endl;
	conf_file << "INDIVIDUAL_NUM " << allGroups->bioRepNum << endl;
	conf_file << "SAMPLE_CNT_PER_GROUP " << totalSampleNum / 2 << endl;
	conf_file << "TOTAL_SAMPLE_CNT " << totalSampleNum << endl;
	
	if (config_bal_design == true)
	{
		conf_file << "Num_of_blocks " << allGroups->bioRepNum << endl;
		conf_file << "Num_of_samples_per_block ";
		bioReplicate *curRep = allGroups->bioReps;
		while (curRep != NULL)
		{
			conf_file << curRep->techRepNum << " ";
			curRep = curRep->next;
		}
		conf_file << endl; //number of technical replicates of each individual

		conf_file << "Array_base_of_blocks ";
		curRep = allGroups->bioReps; int tmpsum = 0;
		while (curRep != NULL)
		{
			conf_file << tmpsum << " ";
			tmpsum += curRep->techRepNum;
			curRep = curRep->next;
		}
		conf_file << endl; //array base of technical replicates of each individual
	} 
	else
	{
		conf_file << "Num_of_blocks 1" << endl;
		conf_file << "Num_of_samples_per_block " << totalSampleNum / 2 << endl; //number of technical replicates of each individual
		conf_file << "Array_base_of_blocks 0" << endl; //array base of technical replicates of each individual
	}

	conf_file << "false_discovery_rate " << config_false_discovery_rate << endl;
	conf_file << "thresh_foldchange_up " << config_thresh_foldchange_up << endl;
	conf_file << "thresh_foldchange_down " << config_thresh_foldchange_down << endl;
	conf_file << "THRESHOLD_MIN_EXPR_COVERAGE " << THRESHOLD_MIN_EXPR_COVERAGE << endl; //at least one group must have mean expression larger than this threshold to make differential expression;

	conf_file.close();

	return;
}


void write_configfile_testtrans(string target_path)
{
	ofstream conf_file;
	string filename;

	filename = target_path + "config_testtrans";
	conf_file.open(filename.c_str());

	conf_file << "GROUP_NUM " << groupNum << endl;
	conf_file << "INDIVIDUAL_NUM " << allGroups->bioRepNum << endl;
	conf_file << "SAMPLE_CNT_PER_GROUP " << totalSampleNum / 2 << endl; //number of samples per group
	conf_file << "COUNT_INTRON_RETENTION true" << endl; //whether to take intron retention into account
	conf_file << "TOTAL_SAMPLE_CNT " << totalSampleNum << endl;

	if (config_bal_design == true)
	{
		conf_file << "Num_of_blocks " << allGroups->bioRepNum << endl;
		conf_file << "Num_of_samples_per_block ";
		bioReplicate *curRep = allGroups->bioReps;
		while (curRep != NULL)
		{
			conf_file << curRep->techRepNum << " ";
			curRep = curRep->next;
		}
		conf_file << endl; //number of technical replicates of each individual

		conf_file << "Array_base_of_blocks ";
		curRep = allGroups->bioReps; int tmpsum = 0;
		while (curRep != NULL)
		{
			conf_file << tmpsum << " ";
			tmpsum += curRep->techRepNum;
			curRep = curRep->next;
		}
		conf_file << endl; //array base of technical replicates of each individual
	}
	else
	{
		conf_file << "Num_of_blocks 1" << endl;
		conf_file << "Num_of_samples_per_block " << totalSampleNum / 2 << endl; //number of technical replicates of each individual
		conf_file << "Array_base_of_blocks 0" << endl; //array base of technical replicates of each individual
	}
	
	
	conf_file << "false_discovery_rate " << config_false_discovery_rate << endl;
	conf_file << "thresh_JSD " << config_thresh_sqrtJSD << endl;
	conf_file << "THRESHOLD_MIN_ASM_COVERAGE " << THRESHOLD_MIN_ASM_COVERAGE << endl; //all groups must have mean expression larger than this threshold to make differential transcription
	
	conf_file.close();

	return;
}

void write_configfile(string target_path)
{
	write_configfile_gtree(target_path);
	write_configfile_testexpr(target_path);
	write_configfile_testtrans(target_path);

	return;
}

