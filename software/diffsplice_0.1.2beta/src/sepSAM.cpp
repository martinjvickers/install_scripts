/*    
 *    sepSAM.cpp		
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


#define UNIX

#ifdef UNIX
#include <fstream>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <vector>
#else
#include <fstream>
#include <stdio.h>
#include <string>
#include <iostream>
#include <stdlib.h> 
#include <sstream>
#include <vector>
#endif

using namespace std;

class reads
{
public:
	string readname;
	string field1;
	string field2;
	string field3;
	string field4;
	string field5;
	string XS;

	reads();
};

class chromosome
{
public:
	string name; //should replace targetChromosome
	long total_num_reads; //should replace readAlignIndex
	vector <reads*> cur_readlist; //current list of reads held in the memory

	chromosome(string chr_name) {name = chr_name; total_num_reads = 0;};
	int add_read(reads *new_read);
	int flush_cur_readlist();
};

long totalReadNum;

vector <chromosome> list_chr;

long global_num_reads_in_memory = 0;
const long LIMIT_NUM_READS_IN_MEMORY = 100000;

string lastTimeName;
string dirPrefix;

ofstream abnormal_reads;
ofstream dropped_reads;

reads::reads()
{
	XS = "+";
}

//add one read to the current read list
int chromosome::add_read(reads *new_read)
{
	++total_num_reads;
	cur_readlist.push_back(new_read);

	return 0;
}

//write current reads held in the memory to output files
int chromosome::flush_cur_readlist()
{
	if (cur_readlist.size() <= 0)
		return 0;

	string outputfilename = dirPrefix + name + ".txt";
	ofstream outputfile;
	
	outputfile.open(outputfilename.c_str(), fstream::app);
	if (!outputfile.is_open())
	{
		cout << "warning: cannot open output file " << outputfilename << "." << endl;
		return 1;
	}
	
	for (unsigned long iLoop = 0; iLoop < cur_readlist.size(); ++iLoop)
	{
		reads *curRead = cur_readlist[iLoop];
		outputfile << curRead->readname << "\t" << curRead->field1 << "\t" << curRead->field2 << "\t" << curRead->field3 << "\t" << curRead->field4 << "\t" << curRead->field5 << "\t" << curRead->XS << endl; 
		delete curRead;
	}

	cur_readlist.clear();
	outputfile.close();
	
	return 0;
}



void process(reads* curRead)
{
	//output a single read
	string chr_name;
	unsigned int tmpLoop;

	chr_name = curRead->field2;

	for (tmpLoop = 0; tmpLoop < list_chr.size(); ++tmpLoop)
	{
		if (chr_name.compare(list_chr[tmpLoop].name) == 0)
			break;
	}

	if (tmpLoop >= list_chr.size())
	{
		//chromosome not found
		chromosome new_chr(chr_name);
		list_chr.push_back(new_chr);
	}

	list_chr[tmpLoop].add_read(curRead);

	++global_num_reads_in_memory;
	if (lastTimeName.compare(curRead->readname) != 0)
	{
		lastTimeName = curRead->readname;
		++totalReadNum;
	}

	if (global_num_reads_in_memory > LIMIT_NUM_READS_IN_MEMORY)
	{
		for (tmpLoop = 0; tmpLoop < list_chr.size(); ++tmpLoop)
		{
			list_chr[tmpLoop].flush_cur_readlist();
		}
		global_num_reads_in_memory = 0;		
	}	

	return;
}


void outputRead()
{
	unsigned int tmpLoop;
	string outputfilename;
	ofstream outputChrname, outputDataStat;
	long totalAlignNum = 0;

	for (tmpLoop = 0; tmpLoop < list_chr.size(); ++tmpLoop)
	{
		list_chr[tmpLoop].flush_cur_readlist();
	}

	outputfilename = dirPrefix + "ChromosomeName.txt";
	outputChrname.open(outputfilename.c_str());
	outputfilename = dirPrefix + "DatasetStat.txt";
	outputDataStat.open(outputfilename.c_str(), fstream::app);

	for (tmpLoop = 0; tmpLoop < list_chr.size(); ++tmpLoop)
	{
		outputChrname << list_chr[tmpLoop].name << '\t' << list_chr[tmpLoop].total_num_reads << endl; 
		totalAlignNum += list_chr[tmpLoop].total_num_reads;
	}

	outputChrname << "#" << totalReadNum << '\t' << totalAlignNum << endl;

	outputDataStat << totalAlignNum << endl;

	outputChrname.close();
	dropped_reads.close();
	abnormal_reads.close();
	outputDataStat.close();

	return;
}


void parse(string inputfilename)
{
	reads *ptr = NULL;

	ifstream inputfile;
	inputfile.open(inputfilename.c_str());

	string name, curName, notUsedField;
	char optField[10000];
	int tmp;
	stringstream iss;
	string info;
	
	while (inputfile >> name)
	{
		if (name[0] == '@')
		{
			getline(inputfile, info);
			continue;
		}

// 		for (tmp = 0; tmp < 99; tmp++)
// 		{
// 			if (name[tmp] == '/' || name[tmp] == '\0')
// 			{
// 				break;
// 			}
// 		}
// 		name[tmp] = '\0'; 

		ptr = new reads;

		ptr->readname = name;
		inputfile >> ptr->field1;
		inputfile >> ptr->field2;
		inputfile >> ptr->field3;
		inputfile >> ptr->field4;
		inputfile >> ptr->field5;

		for (tmp = 0; tmp < 5; ++tmp)
		{
			inputfile >> notUsedField;
		}
		
		getline(inputfile, info);
		
		if (info.empty() == false && info.size() < 10000)
		{
			iss << info;
			while (iss.getline(optField, 10000, '\t'))
			{
				if (optField[0] == 'X' && optField[1] == 'S')
				{
					ptr->XS = optField[5];
				}
			}
			iss.clear();
		}

		process(ptr);
	}


	inputfile.close();

	return;
}


void initialization()
{
	totalReadNum = 0;
	lastTimeName = "XXXX";

	return;
}


int main(int argc, char* argv[])
{
#ifdef UNIX
	if (argc != 3)
	{
		cout << argv[0] << "\t<filename>\t<target_path>" << endl;
		return 1;
	}

	dirPrefix = argv[2];

	if (dirPrefix.compare("./") == 0)
		dirPrefix = "./data/tmp/";
#else
	dirPrefix = "";
#endif

	string outpurfilename;
	outpurfilename = dirPrefix + "DroppedReads.txt";
	dropped_reads.open(outpurfilename.c_str());
	outpurfilename = dirPrefix + "AbnormalReads.txt";
	abnormal_reads.open(outpurfilename.c_str());

	string inputfilename;
#ifdef UNIX
	inputfilename = argv[1];
#else
	inputfilename = "test.sam";
#endif

	initialization();

	parse(inputfilename);

	outputRead();

	return 0;
}


