/*    
 *    parseFragment.cpp		
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



/************************************************************************/
/* Get Fragments From SAM File                                          */
/************************************************************************/

//Should be unique alignments

#define UNIX

#ifdef UNIX
#include <fstream>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <cstdlib> 
#include <unistd.h>
#else
#include <fstream>
#include <stdio.h>
#include <string>
#include <iostream>
#include <stdlib.h> 
#endif

using namespace std;

//#define DEBUG

#define CIGAR_EXPRESSION_FULL
//otherwise only accept M and N

const long MAX = 2000000000;

ofstream exonFile;
ofstream junctionFile;
//ofstream readFile;
ofstream statFile;
ofstream datastatFile;
ofstream crossGeneFile;

ofstream allExonFile;
ofstream allJunctionFile;


class transcript
{
public:
	string transID;
	long start;
	long end;

	transcript *next;

	transcript();
};

class gene
{
public:
	string geneID;
	string chrNm;
	long start;
	long end;

	long db_trans_cnt;

	transcript *transList;

	gene();
};

// gene* geneList[1000000];
// long geneListNum = 0;
// 
// long sortKey[1000000];
// long mergeSort_Larray[1000000];
// long mergeSort_Rarray[1000000];
// gene* mergeSort_LorderedList[1000000];
// gene* mergeSort_RorderedList[1000000];

string dirPrefix;


enum read_type {read_ingene, read_crossgene, read_betweengene, read_crossbetweengene, read_crossIngeneBetweengene};


transcript::transcript()
{
	start = MAX;
	end = 0;

	next = NULL;
}

gene::gene()
{
	start = MAX;
	end = 0;

	db_trans_cnt = 0;

	transList = NULL;
}


long compute_endpoint(long startPosition, char *end)
{
	long i, tmp = 0, endPosition;

	endPosition = startPosition - 1;

	for (i = 0; end[i] != '\0'; i++)
	{
		if (end[i] == 'M')
		{
			endPosition = endPosition + tmp;
			tmp = 0;
		} 
		else if (end[i] == 'N')
		{
			endPosition = endPosition + tmp;
			tmp = 0;
		}
		else if (end[i] >= '0' && end[i] <= '9')
		{
			tmp = tmp * 10 + end[i] - 48;
		}
		else
		{
			tmp = 0;
		}
	}

	return endPosition;
}
// 
// //sort fragment
// void merge(long p, long q, long r)
// {
// 	long n1, n2, i, j, k;
// 
// 	n1 = q - p + 1;
// 	n2 = r - q;
// 
// 	for (i = 1; i <= n1; i++)
// 	{
// 		mergeSort_Larray[i] = sortKey[p + i - 1];
// 		mergeSort_LorderedList[i] = geneList[p + i - 1];
// 	}
// 	for (j = 1; j <= n2; j++)
// 	{
// 		mergeSort_Rarray[j] = sortKey[q + j];
// 		mergeSort_RorderedList[j] = geneList[q + j];
// 	}
// 
// 	mergeSort_Larray[n1 + 1] = MAX;
// 	mergeSort_Rarray[n2 + 1] = MAX;
// 
// 	i = 1;
// 	j = 1;
// 
// 	for (k = p; k <= r; k++)
// 	{
// 		if (mergeSort_Larray[i] <= mergeSort_Rarray[j])
// 		{
// 			sortKey[k] = mergeSort_Larray[i];
// 			geneList[k] = mergeSort_LorderedList[i];
// 
// 			i++;
// 		} 
// 		else
// 		{
// 			sortKey[k] = mergeSort_Rarray[j];
// 			geneList[k] = mergeSort_RorderedList[j];
// 
// 			j++;
// 		}
// 	}
// 
// 	return;
// }
// 
// 
// void mergeSort(long sortList_size)
// {
// 	//non-recursive merge sort for sorting junctions
// 	long m, n, i, r;
// 	m = 1;
// 	n = sortList_size;
// 
// 	while (m <= n)
// 	{
// 		i = 1;
// 		while (i <= n - m)
// 		{
// 			r = (i + 2 * m - 1) < n ? (i + 2 * m - 1) : n;
// 			merge(i, i + m - 1, r);
// 			i = i + 2 * m;
// 		}
// 
// 		m = m * 2;
// 	}
// 
// 	return;
// }
// 
// void parseDB(char* inputfilename)
// {
// 	gene *newGene;
// 	newGene = NULL;
// 
// 	ifstream inputfile;
// 	inputfile.open(inputfilename);
// 
// 	char chrNm[100], lineCategory[100], IDinfo[500], tmpChar[100];
// 	int tmp;
// 	long start, end, iLoop;
// 	string otherInfo;
// 	transcript *newTrans;
// 
// 	for (tmp = 0; tmp < 100; tmp++)
// 	{
// 		chrNm[tmp] = '\0';
// 	}
// 
// 	inputfile >> chrNm;
// 	while (chrNm[0] != '\0')
// 	{
// 		inputfile >> tmpChar;
// 		inputfile >> lineCategory;
// 		inputfile >> start;
// 		inputfile >> end;
// 		inputfile >> tmpChar;
// 		inputfile >> tmpChar;
// 		inputfile >> tmpChar;
// 		inputfile >> IDinfo;
// 		getline(inputfile, otherInfo);
// 
// 		if (strcmp(lineCategory, "gene") == 0)
// 		{
// 			newGene = new gene;
// 			strcpy(newGene->chrNm, chrNm);
// 			newGene->start = start;
// 			newGene->end = end;
// 			newGene->db_trans_cnt = 0;
// 
// 			if (IDinfo[0] == 'I' && IDinfo[1] == 'D' && IDinfo[2] == '=')
// 			{
// 				for (tmp = 0; IDinfo[tmp + 3] != '|' && IDinfo[tmp + 3] != '\0'; tmp++)
// 				{
// 					newGene->geneID[tmp] = IDinfo[tmp + 3];
// 				}
// 				newGene->geneID[tmp] = '\0';
// 			} 
// 			else
// 			{
// 				cout << "Error: abnormal gene ID. Please confirm... ";
// 				cin >> tmpChar;
// 				exit(1);
// 			}
// 
// 			geneListNum++;
// 			geneList[geneListNum] = newGene;
// 		}
// 		else if (strcmp(lineCategory, "transcript") == 0)
// 		{
// 			newTrans = new transcript;
// 
// 			newGene->db_trans_cnt += 1;
// 
// 			newTrans->start = start;
// 			newTrans->end = end;
// 			newTrans->next = newGene->transList;
// 			newGene->transList = newTrans; 
// 
// 			if (IDinfo[0] == 'I' && IDinfo[1] == 'D' && IDinfo[2] == '=')
// 			{
// 				for (tmp = 0; IDinfo[tmp + 3] != ';' && IDinfo[tmp + 3] != '|' && IDinfo[tmp + 3] != '\0'; tmp++)
// 				{
// 					newTrans->transID[tmp] = IDinfo[tmp + 3];
// 				}
// 				newTrans->transID[tmp] = '\0';
// 			} 
// 			else
// 			{
// 				cout << "Error: abnormal gene ID. Please confirm... ";
// 				cin >> tmpChar;
// 				exit(1);
// 			}
// 		}
// 		else
// 		{
// 			//do nothing
// 		}
// 
// 		chrNm[0] = '\0';
// 		inputfile >> chrNm;		
// 	}
// 
// 	inputfile.close();
// 
// 	for (iLoop = 1; iLoop <= geneListNum; iLoop++)
// 	{
// 		sortKey[iLoop] = geneList[iLoop]->end;
// 	}
// 	mergeSort(geneListNum);
// 
// 	for (iLoop = 1; iLoop <= geneListNum; iLoop++)
// 	{
// 		sortKey[iLoop] = geneList[iLoop]->start;
// 	}
// 	mergeSort(geneListNum);
// 
// 	return;
// }
// 
// double getGeneIndex(long position)
// {
// 	double posIndex = 0;
// 	long iLoop;
// 
// 	if (position < geneList[1]->start)
// 	{
// 		posIndex = 0.5;
// 	}
// 	else if (position > geneList[geneListNum]->end)
// 	{
// 		posIndex = geneListNum + 0.5;
// 	}
// 	else
// 	{
// 		for (iLoop = 1; iLoop < geneListNum; iLoop++)
// 		{
// 			if (geneList[iLoop]->start <= position && position <= geneList[iLoop]->end)
// 			{
// 				posIndex = iLoop;
// 			}
// 			else if (geneList[iLoop]->end < position && position < geneList[iLoop + 1]->start)
// 			{
// 				posIndex = iLoop + 0.5;
// 			}
// 		}
// 		if (geneList[geneListNum]->start <= position && position <= geneList[geneListNum]->end)
// 		{
// 			posIndex = geneListNum;
// 		}
// 	}
// 
// 	return posIndex;
// }
// 
// read_type crossGene(long startPosition, long endPosition)
// {
// 	return read_ingene;
// 
// 	//return true if the read is cross gene
// 	long iLoop;
// 	double startIndex = 0, endIndex = 0;
// 
// 	startIndex = getGeneIndex(startPosition);
// 	endIndex = getGeneIndex(endPosition);
// 
// 	if (startIndex == endIndex)
// 	{
// 		if (startIndex == int(startIndex))
// 		{
// 			return read_ingene;
// 		} 
// 		else
// 		{
// 			return read_betweengene;
// 		}
// 	}
// 	else
// 	{
// 		if (startIndex == int(startIndex) && endIndex == int(endIndex))
// 		{
// 			return read_crossgene;
// 		}
// 		else if (startIndex != int(startIndex) && endIndex != int(endIndex))
// 		{
// 			return read_crossbetweengene;
// 		}
// 		else
// 		{
// 			return read_crossIngeneBetweengene;
// 		}
// 	}
// }

void parse(string inputfilename)
{
	ifstream inputfile;
	inputfile.open(inputfilename.c_str());
		
	string name;
	string info, curLine;

	string field1, field2, field3, field4, end, XS;
	long startPoint, endPoint, tmp, strLoop, sign;
	long totalReadNum = 0, crossGeneReadNum = 0, cnt_read_ingene = 0, cnt_read_crossgene = 0, cnt_read_betweengene = 0, cnt_read_crossbetweengene = 0, cnt_read_crossIngeneBetweengene = 0;
	bool spliced, transDirection, flag;
	read_type curReadType;

	long stat_MinStart = MAX, stat_MaxEnd = 0;
	
	while (inputfile >> name)
	{
		inputfile >> field1;
		inputfile >> field2;
		inputfile >> startPoint;
		inputfile >> field4;
		inputfile >> end;
		inputfile >> XS;
		getline(inputfile, info);

#ifndef CIGAR_EXPRESSION_FULL
		flag = false;
		for (strLoop = 0; strLoop < end.size(); ++strLoop)
		{
			if (end[strLoop] != 'M' && end[strLoop] != 'N' && (end[strLoop] < '0' || end[strLoop] > '9'))
			{
				flag = true;
				break;
			}
		}
		if (flag == true)
		{
			end = "*";
		}
#endif

		if (end.compare("*") != 0)
		{
			curReadType = read_ingene; //crossGene(startPoint, compute_endpoint(startPoint, end));
			
			totalReadNum++;

			if (curReadType == read_ingene)
			{
				cnt_read_ingene++;
			}
			else if (curReadType == read_crossgene)
			{
				cnt_read_crossgene++;
			}
			else if (curReadType == read_betweengene)
			{
				cnt_read_betweengene++;
			}
			else if (curReadType == read_crossbetweengene)
			{
				cnt_read_crossbetweengene++;
			}
			else if (curReadType == read_crossIngeneBetweengene)
			{
				cnt_read_crossIngeneBetweengene++;
			}

			if (curReadType == read_crossgene || curReadType == read_crossbetweengene || curReadType == read_crossIngeneBetweengene)
			{
				crossGeneReadNum++;
				crossGeneFile << name << "\t" << field1 << "\t" << field2 << "\t" << startPoint << "\t" << field4 << "\t" << end << endl;
			}
			else if (curReadType == read_ingene || curReadType == read_betweengene)
			{
				tmp = 0;
				sign = 1;
				spliced = false;
				if (XS.compare("-") == 0)
					transDirection = false;
				else
					transDirection = true;
				
				if (startPoint < stat_MinStart)
					stat_MinStart = startPoint;

				for (strLoop = 0; strLoop < end.size(); ++strLoop)
				{
					if (end[strLoop] == 'N')
					{
						spliced = true;
						break;
					}
				}

				//readFile << field2 << '\t' << startPoint << '\t' << startPoint << '\t' << spliced << endl;

				for (strLoop = 0; strLoop < end.size(); ++strLoop)
				{
					if (end[strLoop] == 'M')
					{
						endPoint = startPoint + tmp * sign;
						if (startPoint < 0)
						{
							exit(1);
						}
						exonFile << field2 << '\t' << startPoint << '\t' << endPoint << '\t' << spliced << endl;
						allExonFile << field2 << '\t' << startPoint << '\t' << endPoint << '\t' << spliced << endl;
						startPoint = endPoint;
						tmp = 0;
						sign = 1;
					} 
					else if (end[strLoop] == 'N')
					{
						endPoint = startPoint + tmp * sign;
						if (startPoint < 0)
						{
							exit(1);
						}
						junctionFile << field2 << '\t' << startPoint - 1 << '\t' << endPoint << "\t" << transDirection << endl;
						allJunctionFile << field2 << '\t' << startPoint - 1 << '\t' << endPoint << "\t" << transDirection << endl;
						startPoint = endPoint;
						tmp = 0;
						sign = 1;
					}
					else if (end[strLoop] == 'I')
					{
						tmp = 0;
						sign = 1;
					}
					else if (end[strLoop] == 'D' || end[strLoop] == 'S')
					{
						endPoint = startPoint + tmp * sign;
						startPoint = endPoint;
						tmp = 0;
						sign = 1;
					}
					else if (end[strLoop] == '-')
					{
						exit(1);
					}
					else if (end[strLoop] >= '0' && end[strLoop] <= '9')
					{
						tmp = tmp * 10 + end[strLoop] - 48;
					}
					else
					{
						//					endPoint = startPoint + tmp * sign;
						//					if (startPoint < 0)
						//					{
						//						exit(1);
						//					}
						//					startPoint = endPoint;
						tmp = 0;
						sign = 1;
					}
				}

				if (endPoint > stat_MaxEnd)
					stat_MaxEnd = endPoint;
			}
		}
	}
	
	statFile << stat_MinStart << '\t' << stat_MaxEnd << endl;
	datastatFile << totalReadNum << "\t" << crossGeneReadNum << "\t" << cnt_read_ingene << "\t" << cnt_read_crossgene << "\t" << cnt_read_betweengene << "\t" << cnt_read_crossbetweengene << "\t" << cnt_read_crossIngeneBetweengene << endl;

	return;
}



int main(int argc, char* argv[])
{
	if (argc == 5)
	{
		dirPrefix = "./tmp/";
	}
	else if (argc == 6)
	{
		dirPrefix = argv[5];
	}
	else
	{
		cout << argv[0] << "\t<readfilepath>\t<chr>\t<outputfile_suffix>\t<geneDB_filename>\t<output_folder>" << endl;
		return 1;
	}

	string inputfilename, outputfilename, comd;
	
	outputfilename = dirPrefix + "junction.txt";
	junctionFile.open(outputfilename.c_str());
	outputfilename = dirPrefix + "exon.txt";
	exonFile.open(outputfilename.c_str());
	//sprintf(outputfilename, "%sread.txt", dirPrefix);
	//readFile.open(outputfilename);
	outputfilename = dirPrefix + "stat.txt";
	statFile.open(outputfilename.c_str(), fstream::app);
	outputfilename = dirPrefix + "datastat.txt";
	datastatFile.open(outputfilename.c_str(), fstream::app);
	outputfilename = dirPrefix + "crossgenedata.txt";
	crossGeneFile.open(outputfilename.c_str(), fstream::app);
	
	outputfilename = dirPrefix + "allJunction.txt";
	allJunctionFile.open(outputfilename.c_str(), fstream::app);
	outputfilename = dirPrefix + "allExon.txt";
	allExonFile.open(outputfilename.c_str(), fstream::app);

// 	strcpy(inputfilename, argv[4]);
// 	parseDB(inputfilename);

	datastatFile << argv[3] << "\t";

	inputfilename = argv[1];
	inputfilename += argv[2];
	inputfilename += ".txt";
#ifdef DEBUG
	cout << inputfilename << "... ";
#endif
	parse(inputfilename);

	junctionFile.close();
	exonFile.close();
	//readFile.close();
	statFile.close();
	datastatFile.close();
	crossGeneFile.close();

	allJunctionFile.close();
	allExonFile.close();

	comd = "sort -n +1 -4 " + dirPrefix + "junction.txt > " + dirPrefix + "junction_" + argv[3] + ".txt";
	system(comd.c_str());
	comd = "sort -n +1 -4 " + dirPrefix + "exon.txt > " + dirPrefix + "exon_" + argv[3] + ".txt";
	system(comd.c_str());
// 	sprintf(comd, "sort -n +1 -4 %sread.txt > %sread_%s.txt", dirPrefix, dirPrefix, argv[3]);
// 	system(comd);

	comd = "rm -f -r " + dirPrefix + "junction.txt";
	system(comd.c_str());
	comd = "rm -f -r " + dirPrefix + "exon.txt";
	system(comd.c_str());
// 	sprintf(comd, "rm -r %sread.txt", dirPrefix);
// 	system(comd);

	string readCount_total;
	ifstream inputDataStat;
	inputfilename = argv[1]; inputfilename += "DatasetStat.txt";
	inputDataStat.open(inputfilename.c_str());
	getline(inputDataStat, readCount_total);
	inputDataStat.close();

	ofstream outputDataStat;
	outputfilename = dirPrefix + "DatasetStat.txt";
	outputDataStat.open(outputfilename.c_str(), fstream::app);
	outputDataStat << readCount_total << endl;
	outputDataStat.close();


	return 0;
}



