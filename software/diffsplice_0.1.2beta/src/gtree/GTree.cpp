/*    
 *    GTree.cpp		
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

#include "GTree.h"


/************************************************************************/
/* functions                                                            */
/************************************************************************/

//clear a vector and release its memory
template <class T>
void free_vector(T &t)
{
	T tmp;
	tmp.swap(t);
}

string itostr(long t)
{
	ostringstream oss;
	oss << t;
	return oss.str();
}


//*/
//sort fragment
void merge_JunctionSort(long p, long q, long r)
{
	long n1, n2, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	for (i = 1; i <= n1; i++)
	{
		mergeSort_Larray[i] = sortKey_Junction[p + i - 1];
		mergeSort_LorderedList[i] = sortList_Junction[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray[j] = sortKey_Junction[q + j];
		mergeSort_RorderedList[j] = sortList_Junction[q + j];
	}

	mergeSort_Larray[n1 + 1] = MAX_CHR_LENGTH * 2;
	mergeSort_Rarray[n2 + 1] = MAX_CHR_LENGTH * 2;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray[i] <= mergeSort_Rarray[j])
		{
			sortKey_Junction[k] = mergeSort_Larray[i];
			sortList_Junction[k] = mergeSort_LorderedList[i];

			i++;
		} 
		else
		{
			sortKey_Junction[k] = mergeSort_Rarray[j];
			sortList_Junction[k] = mergeSort_RorderedList[j];

			j++;
		}
	}

	return;
}


void mergeSort_JunctionSort(long sortList_size)
{
	//non-recursive merge sort for sorting junctions
	long m, n, i, r;
	m = 1;
	n = sortList_size;

	while (m <= n)
	{
		i = 1;
		while (i <= n - m)
		{
			r = (i + 2 * m - 1) < n ? (i + 2 * m - 1) : n;
			merge_JunctionSort(i, i + m - 1, r);
			i = i + 2 * m;
		}

		m = m * 2;
	}

	return;
}


long partition(long p,long r, double *sortArray)
{
	long i, j;
	double x, tmp;

	//randomized partition
	i = long(p + (double)rand()/ (RAND_MAX) * (r - p));

	if (sortArray[r] != sortArray[i])
	{
		tmp = sortArray[r];
		sortArray[r] = sortArray[i];
		sortArray[i] = tmp;
	}

	x = sortArray[r];
	i = p - 1;

	for(j = p; j < r; j++)
	{
		if (sortArray[j] <= x)
		{
			i++;

			if (sortArray[i] != sortArray[j])
			{
				tmp = sortArray[j];
				sortArray[j] = sortArray[i];
				sortArray[i] = tmp;
			}
		}
	}

	if (sortArray[r] != sortArray[i+1])
	{
		tmp = sortArray[r];
		sortArray[r] = sortArray[i+1];
		sortArray[i+1]=tmp;
	}

	return i+1;
}


void quicksort(double *sortArray, long length, long *quicksortStack)
{
	long top = 0, p, r, q;

	quicksortStack[top++] = 1;
	quicksortStack[top++] = length;

	while (top != 0)
	{
		r = quicksortStack[--top];
		p = quicksortStack[--top];

		if(p>=r)
			continue;

		q = partition(p, r, sortArray);

		quicksortStack[top++] = p;
		quicksortStack[top++] = q - 1;

		quicksortStack[top++] = q + 1;
		quicksortStack[top++] = r;
	}

	return;
}



//RangeJunctionList Queue
//for separating dependent paths
void juncPathQueue_initialization()
{
	while (juncPathQueue.empty() == false)
	{
		juncPathQueue.pop();
	}

	return;
}

bool juncPathQueue_enqueue(JuncGraphPath *x)
{
	juncPathQueue.push(x);

	return true;
}

JuncGraphPath* juncPathQueue_dequeue()
{
	if (juncPathQueue.empty() == true)
	{
		return NULL;
	} 
	else
	{
		JuncGraphPath *next_ele = juncPathQueue.front();
		juncPathQueue.pop();
		return next_ele;
	}
}


/************************************************************************/
/* JUNCTION GRAPH                                                       */
/************************************************************************/


JuncGraphEdge::JuncGraphEdge()
{
	linkedVertex = NULL;
	next = NULL;
}

JuncGraphVertex::JuncGraphVertex()
{
	corresJunc = NULL;
	edges = NULL;
	next = NULL;
	traversed = false;
	hasInEdge = false;
	hasOutEdge = false;
	vertexType = normal;
}

JuncGraphVertex::~JuncGraphVertex()
{
	if (corresJunc != NULL)
	{
		delete corresJunc;
	}

	JuncGraphEdge *delEdge;

	while (edges != NULL)
	{
		delEdge = edges;
		edges = delEdge->next;
		delete delEdge;
	}
}

JuncGraphPath::JuncGraphPath()
{
	pathJuncList = NULL;
	arrivedVertex = NULL;
}

JuncGraphPath::~JuncGraphPath()
{
	arrivedVertex = NULL;
	delete pathJuncList;
}

JuncGraph::JuncGraph()
{
	vertices = NULL;
}

JuncGraph::~JuncGraph()
{
	JuncGraphVertex *delVertex;

	while (vertices != NULL)
	{
		delVertex = vertices;
		vertices = delVertex->next;
		delete delVertex;
	}
}


//fragment graph
//assume junctions have been sorted based on position

//separate genes 
RangeJunctionList* separateGene(RangeJunctionList* origList, bool &separable)
{
	//separate genes
	//return a set of lists, each of which corresponds to a gene
	//IMPORTANT: RangeJunctionList is head-inserting, for the convenience of stack

	if (origList == NULL)
	{
		separable = false;
		return NULL;
	}

	RangeJunctionList *resultList, *curList;
	resultList = NULL;
	curList = NULL;

	rangeJunction *curJunc, *curListTail;
	curJunc = NULL;
	curListTail = NULL;

	long endBoard = -2;
	int numIndep = 0; //number of independent regions

	curJunc = origList->list;
	while (curJunc != NULL)
	{
		if (curJunc->junc->start > endBoard+1)
		{
			//get an independent region
			numIndep++;

			//omit the first list, it is empty
			if (curList == NULL)
			{
				//do nothing
			} 
			else
			{
				curList->rangeHigh = endBoard;
			}

			curList = new RangeJunctionList;
			curList->nextList = resultList;

			// 			if (resultList == NULL)
			// 			{
			// 				//first list
			// 				curList->rangeLow = origList->rangeLow;
			// 			} 
			// 			else
			// 			{
			// 				curList->rangeLow = endBoard;
			// 			}
			curList->rangeLow = curJunc->junc->start;

			resultList = curList;

			curListTail = NULL;
		}

		//process current fragment
		if (curListTail == NULL)
		{
			curList->list = curJunc;
			curListTail = curJunc;
		}
		else
		{
			curListTail->next = curJunc;
			curListTail = curJunc;
		}

		if (curJunc->junc->end > endBoard)
		{
			endBoard = curJunc->junc->end;
		}

		curJunc = curJunc->next;
		curListTail->next = NULL;
	}

	if (curList == NULL)
	{
		cout << "curList == NULL" << endl;
		exit(1);
	} 
	else
	{
		//curList->rangeHigh = origList->rangeHigh;
		curList->rangeHigh = endBoard;
	}

	if (numIndep == 1)
		separable = false; //cannot be separated into multiple regions
	else if (numIndep > 1)
		separable = true; //can be separated
	else
	{
		cout << "numGene == 0" << endl;
		exit(1);
	}

	origList->list = NULL;
	delete origList;

	return resultList;
}




//separate independent paths

// bool compatibleJunctions(rangeJunction* junctionA, rangeJunction* junctionB)
// {
// 	//check the compatibility of two junctions
// 	if (junctionA->junc->end < junctionB->junc->start)
// 	{
// 		return true;
// 	}
// 	else if (junctionA->junc->start > junctionB->junc->end)
// 	{
// 		return true;
// 	}
// 	else
// 		return false;
// }
bool compatibleJunctions(rangeJunction* junctionA, rangeJunction* junctionB)
{
	//check the compatibility of two junctions
// 	if (fabs(junctionB->junc->start - junctionA->junc->end) <= 1)
// 	{
// 		return true;
// 	}
// 	else if (fabs(junctionB->junc->end - junctionA->junc->start) <= 1)
// 	{
// 		return true;
// 	}
// 	else
// 		return false;

	if (junctionA->junc->end <= junctionB->junc->start && junctionB->junc->start - junctionA->junc->end <= 1)
	{
		return true;
	}
	else if (junctionA->junc->start >= junctionB->junc->end && junctionA->junc->start - junctionB->junc->end <= 1)
	{
		return true;
	}
	else
		return false;
}

void constructJuncGraph_undirected(RangeJunctionList* origList, bool virtualSE)
{
	//construct fragment graph based on given junctions
	//add virtual start and virtual end if virtualSE is true 

	JuncGraphVertex *curVertex, *tailVertex, *edgeVertex;
	JuncGraphEdge *newEdge;
	rangeJunction *curJunc;

	//build vertices

	curJunc = origList->list;
	while (curJunc != NULL)
	{
		curVertex = new JuncGraphVertex;
		curVertex->corresJunc = curJunc;
		curJunc = curJunc->next;
		curVertex->corresJunc->next = NULL;

		if (junctionGraph->vertices == NULL)
		{
			junctionGraph->vertices = curVertex;
			tailVertex = curVertex;
		} 
		else
		{
			tailVertex->next = curVertex;
			tailVertex = curVertex;
		}
	}

	if (virtualSE == true)
	{
		curVertex = new JuncGraphVertex;
		curVertex->vertexType = virStart;
		curVertex->next = junctionGraph->vertices;
		junctionGraph->vertices = curVertex;

		curVertex = new JuncGraphVertex;
		curVertex->vertexType = virEnd;
		tailVertex->next = curVertex;
		tailVertex = curVertex;
	}

	//build edges

	curVertex = junctionGraph->vertices;
	while (curVertex != NULL)
	{
		if (curVertex->vertexType == normal)
		{
			edgeVertex = curVertex->next;
			while (edgeVertex != NULL)
			{
				if (edgeVertex->vertexType == normal)
				{
					if (compatibleJunctions((curVertex->corresJunc), (edgeVertex->corresJunc)) == true)
					{
						//create an edge for curVertex
						newEdge = new JuncGraphEdge;
						newEdge->linkedVertex = edgeVertex;

						newEdge->next = curVertex->edges;
						curVertex->edges = newEdge;

						//create an edge for edgeVertex
						newEdge = new JuncGraphEdge;
						newEdge->linkedVertex = curVertex;

						newEdge->next = edgeVertex->edges;
						edgeVertex->edges = newEdge;

						edgeVertex->hasInEdge = true;
						curVertex->hasOutEdge = true;
					} 
					else
					{
						//break; 
					}
				}

				edgeVertex = edgeVertex->next;
			}
		}

		curVertex = curVertex->next;
	}

	if (virtualSE == true)
	{
		curVertex = junctionGraph->vertices;
		edgeVertex = curVertex->next;
		while (edgeVertex != NULL)
		{
			if (edgeVertex->vertexType == normal && edgeVertex->hasInEdge == false)
			{
				//create an edge for curVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = edgeVertex;

				newEdge->next = curVertex->edges;
				curVertex->edges = newEdge;

				//create an edge for edgeVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = curVertex;

				newEdge->next = edgeVertex->edges;
				edgeVertex->edges = newEdge;

			}

			edgeVertex = edgeVertex->next;
		}

		curVertex = tailVertex;
		edgeVertex = junctionGraph->vertices->next;
		while (edgeVertex != NULL)
		{
			if (edgeVertex->vertexType == normal && edgeVertex->hasOutEdge == false)
			{
				//create an edge for curVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = edgeVertex;

				newEdge->next = curVertex->edges;
				curVertex->edges = newEdge;

				//create an edge for edgeVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = curVertex;

				newEdge->next = edgeVertex->edges;
				edgeVertex->edges = newEdge;

			}

			edgeVertex = edgeVertex->next;
		}
	}


	return;
}

void constructJuncGraph_directed(RangeJunctionList* origList)
{
	//construct fragment graph based on given junctions

	JuncGraphVertex *curVertex, *tailVertex, *edgeVertex;
	JuncGraphEdge *newEdge;
	rangeJunction *curJunc;

	//build vertices

	curJunc = origList->list;
	while (curJunc != NULL)
	{
		curVertex = new JuncGraphVertex;
		curVertex->corresJunc = curJunc;
		curJunc = curJunc->next;
		curVertex->corresJunc->next = NULL;

		if (junctionGraph->vertices == NULL)
		{
			junctionGraph->vertices = curVertex;
			tailVertex = curVertex;
		} 
		else
		{
			tailVertex->next = curVertex;
			tailVertex = curVertex;
		}
	}

	//build edges

	curVertex = junctionGraph->vertices;
	while (curVertex != NULL)
	{
		edgeVertex = curVertex->next;
		while (edgeVertex != NULL)
		{
			if (compatibleJunctions((curVertex->corresJunc), (edgeVertex->corresJunc)) == true)
			{
				//create an edge for curVertex
				newEdge = new JuncGraphEdge;
				newEdge->linkedVertex = edgeVertex;

				newEdge->next = curVertex->edges;
				curVertex->edges = newEdge;

				// 				//create an edge for edgeVertex
				// 				newEdge = new JuncGraphEdge;
				// 				newEdge->linkedVertex = curVertex;
				// 
				// 				newEdge->next = edgeVertex->edges;
				// 				edgeVertex->edges = newEdge;

				edgeVertex->hasInEdge = true;
				curVertex->hasOutEdge = true;
			} 
			else
			{
				//break; 
			}
			edgeVertex = edgeVertex->next;
		}

		curVertex = curVertex->next;
	}


	return;
}


void DFS_visit(JuncGraphVertex *u)
{
	//visit a vertex u during DFS
	u->traversed = true;

	JuncGraphVertex *v;
	JuncGraphEdge *curEdge;
	curEdge = u->edges;
	while (curEdge != NULL)
	{
		v = curEdge->linkedVertex;
		if (v->traversed == false)
		{
			DFS_visit(v);
		}
		curEdge = curEdge->next;
	}

	//add u into sortList_Junction
	sortList_Junction_Num++;
	if (sortList_Junction_Num >= sortList_Junction.size())
	{
		sortList_Junction.resize(sortList_Junction.size() + DEFAULT_JUNCTION_NUM, NULL);
		sortKey_Junction.resize(sortKey_Junction.size() + DEFAULT_JUNCTION_NUM, 0);
		mergeSort_Larray.resize(mergeSort_Larray.size() + DEFAULT_JUNCTION_NUM, 0);
		mergeSort_Rarray.resize(mergeSort_Rarray.size() + DEFAULT_JUNCTION_NUM, 0);
		mergeSort_LorderedList.resize(mergeSort_LorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
		mergeSort_RorderedList.resize(mergeSort_RorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
	}
	sortList_Junction[sortList_Junction_Num] = u->corresJunc;
	u->corresJunc = NULL;

	return;
}

RangeJunctionList* search_conn_comp(int &numCC, long rangeLow, long rangeHigh)
{
	//search all connected component within junctionGraph
	JuncGraphVertex *curVertex;
	RangeJunctionList *resultList, *curList;
	long i;

	numCC = 0; //number of connected component

	resultList = NULL;
	sortList_Junction_Num = 0;
	curVertex = junctionGraph->vertices;
	while (curVertex != NULL)
	{
		if (curVertex->traversed == false)
		{
			DFS_visit(curVertex);

			numCC++;

			curList = new RangeJunctionList;
			// 			curList->rangeLow = rangeLow;
			// 			curList->rangeHigh = rangeHigh;


			//get a connected component
			for (i = 1; i <= sortList_Junction_Num; i++)
			{
				sortKey_Junction[i] = sortList_Junction[i]->junc->end;
			}
			mergeSort_JunctionSort(sortList_Junction_Num);

			curList->rangeHigh = sortList_Junction[sortList_Junction_Num]->junc->end;

			for (i = 1; i <= sortList_Junction_Num; i++)
			{
				sortKey_Junction[i] = sortList_Junction[i]->junc->start;
			}
			mergeSort_JunctionSort(sortList_Junction_Num);

			curList->rangeLow = sortList_Junction[1]->junc->start;

			curList->nextList = resultList;
			resultList = curList;

			for (i = sortList_Junction_Num; i >= 1; i--)
			{
				sortList_Junction[i]->next = curList->list;
				curList->list = sortList_Junction[i];
			}			

			sortList_Junction_Num = 0;
		}

		curVertex = curVertex->next;
	}

	return resultList;
}

//separate independent regions
RangeJunctionList* separateIndepRegion(RangeJunctionList* origList, bool &separable)
{
	//separate independent regions
	//return a set of lists, each of which corresponds to an independent region
	//IMPORTANT: RangeJunctionList is head-inserting, for the convenience of stack

	if (origList == NULL)
	{
		separable = false;
		return NULL;
	}

	junctionGraph = new JuncGraph;
	constructJuncGraph_undirected(origList, true);

	RangeJunctionList *resultList, *curList;
	resultList = NULL;
	curList = NULL;

	rangeJunction *curJunc, *curListTail;
	curJunc = NULL;
	curListTail = NULL;
	
	JuncGraphVertex *curVertex;
	JuncGraphEdge *curEdge;

	long endBoard = 0, geneEnd = 0;
	int numIndep = 0; //number of independent regions

	//handling alternative start
	curEdge = junctionGraph->vertices->edges;
	while (curEdge != NULL)
	{
		if (curEdge->linkedVertex->corresJunc->junc->start > endBoard)
		{
			endBoard = curEdge->linkedVertex->corresJunc->junc->start;
		}
		curEdge = curEdge->next;
	}

	curVertex = junctionGraph->vertices->next;
	while (curVertex != NULL)
	{
		if (curVertex->vertexType == normal)
		{
			curJunc = curVertex->corresJunc;
			curVertex->corresJunc = NULL;

			if (curJunc->junc->start >= endBoard)
			{
				//omit the first list, it is empty
				if (curList != NULL)
				{
					if (endBoard < MAX_CHR_LENGTH)
						curList->rangeHigh = endBoard;
					else
						curList->rangeHigh = geneEnd;
				}

				curList = NULL;
			}

			if (curList == NULL)
			{
				//get an independent region
				numIndep++;

				curList = new RangeJunctionList;
				curList->nextList = resultList;
				curList->rangeLow = curJunc->junc->start;
				resultList = curList;

				curListTail = NULL;
			}

			//process current fragment
			if (curListTail == NULL)
			{
				curList->list = curJunc;
				curListTail = curJunc;
			}
			else
			{
				curListTail->next = curJunc;
				curListTail = curJunc;
			}
			curListTail->next = NULL;

			
			if (curVertex->edges->linkedVertex->vertexType == virEnd)
			{
				endBoard = MAX_CHR_LENGTH;
				if (curJunc->junc->end > geneEnd)
				{
					geneEnd = curJunc->junc->end;
				}
			}
			else if (curJunc->junc->end > endBoard)
			{
				endBoard = curJunc->junc->end;
			}
		}

		curVertex = curVertex->next;
	}

	if (curList == NULL)
	{
		cout << "curList == NULL" << endl;
		exit(1);
	} 
	else
	{
		//curList->rangeHigh = origList->rangeHigh;
		if (endBoard < MAX_CHR_LENGTH)
			curList->rangeHigh = endBoard;
		else
			curList->rangeHigh = geneEnd;
	}


	delete junctionGraph;

	if (numIndep == 1)
		separable = false; //cannot be separated into multiple regions
	else if (numIndep > 1)
		separable = true; //can be separated
	else
	{
		cout << "numIndep == 0" << endl;
		exit(1);
	}

	origList->list = NULL;
	delete origList;

	return resultList;
}

////separate independent regions
//RangeJunctionList* separateIndepRegion(RangeJunctionList* origList, bool &separable)
//{
//	//separate independent regions
//	//return a set of lists, each of which corresponds to an independent region
//	//IMPORTANT: RangeJunctionList is head-inserting, for the convenience of stack
//
//	if (origList == NULL)
//	{
//		separable = false;
//		return NULL;
//	}
//
//	RangeJunctionList *resultList, *curList;
//	resultList = NULL;
//	curList = NULL;
//
//	rangeJunction *curJunc, *curListTail;
//	curJunc = NULL;
//	curListTail = NULL;
//
//	long endBoard = 0;
//	int numIndep = 0; //number of independent regions
//
//	curJunc = origList->list;
//	while (curJunc != NULL)
//	{
//		if (curJunc->junc->start >= endBoard)
//		{
//			//get an independent region
//			numIndep++;
//
//			//omit the first list, it is empty
//			if (curList == NULL)
//			{
//				//do nothing
//			} 
//			else
//			{
//				curList->rangeHigh = endBoard;
//			}
//
//			curList = new RangeJunctionList;
//			curList->nextList = resultList;
//
//			// 			if (resultList == NULL)
//			// 			{
//			// 				//first list
//			// 				curList->rangeLow = origList->rangeLow;
//			// 			} 
//			// 			else
//			// 			{
//			// 				curList->rangeLow = endBoard;
//			// 			}
//			curList->rangeLow = curJunc->junc->start;
//
//			resultList = curList;
//
//			curListTail = NULL;
//		}
//
//		//process current fragment
//		if (curListTail == NULL)
//		{
//			curList->list = curJunc;
//			curListTail = curJunc;
//		}
//		else
//		{
//			curListTail->next = curJunc;
//			curListTail = curJunc;
//		}
//
//		if (curJunc->junc->end > endBoard)
//		{
//			endBoard = curJunc->junc->end;
//		}
//
//		curJunc = curJunc->next;
//		curListTail->next = NULL;
//	}
//
//	if (curList == NULL)
//	{
//		cout << "curList == NULL" << endl;
//		exit(1);
//	} 
//	else
//	{
//		//curList->rangeHigh = origList->rangeHigh;
//		curList->rangeHigh = endBoard;
//	}
//
//	if (numIndep == 1)
//		separable = false; //cannot be separated into multiple regions
//	else if (numIndep > 1)
//		separable = true; //can be separated
//	else
//	{
//		cout << "numIndep == 0" << endl;
//		exit(1);
//	}
//
//	origList->list = NULL;
//	delete origList;
//
//	return resultList;
//}

RangeJunctionList* separateIndepPath(RangeJunctionList* origList, bool &separable)
{
	//separate independent paths 
	//return a set of lists, each of which corresponds to an independent path
	//return decomposable
	if (origList == NULL)
	{
		separable = false;
		return NULL;
	}

	junctionGraph = new JuncGraph;
	constructJuncGraph_undirected(origList, false);
	origList->list = NULL;

	RangeJunctionList *resultList;
	int numCC = 0;
	resultList = search_conn_comp(numCC, origList->rangeLow, origList->rangeHigh);

	delete junctionGraph;

	if (numCC == 1)
		separable = false; //cannot be separated into multiple paths
	else if (numCC > 1)
		separable = true; //can be separated
	else
	{
		cout << "numCC == 0";
		exit(1);
	}

	origList->list = NULL;
	delete origList;

	return resultList;
}


//separate dependent paths

RangeJunctionList* separateDepPath(RangeJunctionList* origList, bool &separable)
{
	//separate dependent paths 
	//return a set of lists, each of which corresponds to a dependent path
	//return decomposable
	if (origList == NULL)
	{
		separable = false;
		return origList;
	}

	rangeJunction *countJunc;
	long junctionCntinList = 0;
	countJunc = origList->list;
	while (countJunc != NULL)
	{
		if (countJunc->junc->type == frag_junction)
		{
			junctionCntinList++;
		}
		countJunc = countJunc->next;
	}
	if (junctionCntinList > MAXJUNCNUMINDEPPATH)
	{
		outfile_not_enumerated << origList->rangeLow << '\t' << origList->rangeHigh << "\t" << junctionCntinList << endl;
		countJunc = origList->list;
		while (countJunc != NULL)
		{
			if (countJunc->junc->type == frag_exon)
			{
				outfile_not_enumerated << "1\t";
			} 
			else
			{
				outfile_not_enumerated << "0\t";
			}
			outfile_not_enumerated << countJunc->junc->start << "\t" << countJunc->junc->end;
			for (int tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
			{
				outfile_not_enumerated << "\t" << countJunc->junc->support[tmp];
			}
			outfile_not_enumerated << endl;
			
			countJunc = countJunc->next;
		}
		outfile_not_enumerated << endl;

		separable = false;
		return origList;
	}

	RangeJunctionList *backupList;
	backupList = origList->clone();
	junctionGraph = new JuncGraph;
	constructJuncGraph_directed(backupList);
	backupList->list = NULL;

	RangeJunctionList *resultList, *newList;
	rangeJunction *newJunc, *curJunc;
	int numPath = 0;
	long startPosition, endPosition;
	bool pathExtended;

	JuncGraphVertex *curVertex;
	JuncGraphEdge *curEdge;
	JuncGraphPath *curPath, *newPath;
	//find start vertex and end vertex

	startPosition = origList->rangeLow;
	endPosition = origList->rangeHigh;

	resultList = NULL;
	juncPathQueue_initialization();


	//build initial queue
	curVertex = junctionGraph->vertices;
	while (curVertex != NULL)
	{
		//		if (curVertex->corresJunc->junc->start == startPosition)
		if (curVertex->hasInEdge == false)
		{
			newList = new RangeJunctionList;
			newList->rangeLow = curVertex->corresJunc->junc->start;
			newList->rangeHigh = curVertex->corresJunc->junc->end;
			newList->transDirection = curVertex->corresJunc->junc->transDirection;

			newJunc = new rangeJunction;
			newJunc->junc = curVertex->corresJunc->junc;
			newList->list = newJunc;
			newList->listtail = newJunc;

			newPath = new JuncGraphPath;
			newPath->arrivedVertex = curVertex;
			newPath->pathJuncList = newList;

			juncPathQueue_enqueue(newPath);
		}
		curVertex = curVertex->next;
	}

	//main part
	curPath = juncPathQueue_dequeue();
	while (curPath != NULL)
	{
		//		if (curPath->arrivedVertex->corresJunc->junc->end == endPosition)
		pathExtended = false;

		if (curPath->arrivedVertex->edges != NULL)
		{
			curEdge = curPath->arrivedVertex->edges;
			while (curEdge != NULL)
			{
				if (curEdge->linkedVertex->corresJunc->junc->end >= curPath->arrivedVertex->corresJunc->junc->end)
				{
					if (curPath->pathJuncList->transDirection == undetermined || curEdge->linkedVertex->corresJunc->junc->transDirection == undetermined
						|| curEdge->linkedVertex->corresJunc->junc->transDirection == curPath->pathJuncList->transDirection)
					{
						pathExtended = true;
						newList = curPath->pathJuncList->clone();

						newJunc = new rangeJunction;
						newJunc->junc = curEdge->linkedVertex->corresJunc->junc;
						//newJunc->next = newList->list;
						//newList->list = newJunc;
						newList->listtail->next = newJunc;
						newList->listtail = newJunc;

						if (newList->transDirection == undetermined && (newJunc->junc->transDirection == sense || newJunc->junc->transDirection == antisense))
						{
							newList->transDirection = newJunc->junc->transDirection;
						}

						if (newJunc->junc->start < newList->rangeLow)
						{
							newList->rangeLow = newJunc->junc->start;
						}
						if (newJunc->junc->end > newList->rangeHigh)
						{
							newList->rangeHigh = newJunc->junc->end;
						}

						newPath = new JuncGraphPath;
						newPath->arrivedVertex = curEdge->linkedVertex;
						newPath->pathJuncList = newList;

						juncPathQueue_enqueue(newPath);
					}
				}

				curEdge = curEdge->next;
			}			
		}

		if (pathExtended == false)
		{
			//current path arrives destination
			numPath++;

			curPath->pathJuncList->nextList = resultList;
			resultList = curPath->pathJuncList;			
		}

		curPath->pathJuncList = NULL;
		delete curPath;

		curPath = juncPathQueue_dequeue();
	}


	delete junctionGraph;

	if (numPath == 1)
		separable = false; //cannot be separated into multiple paths
	else if (numPath > 1)
		separable = true; //can be separated
	else
	{
		cout << origList->rangeLow << " - " << origList->rangeHigh << "  numPath == 0" << endl;
		exit(1);
	}

	delete backupList;

	if (separable == true)
	{
		return resultList;
	} 
	else
	{
		return origList;
	}
}






/************************************************************************/
/* GENOME TREE	                                                        */
/************************************************************************/


alter_junction::alter_junction()
{
	juncInfo = NULL;
//	category = 0;
//	proportion = 0.0;
	next = NULL;
}


fragment::fragment()
{
	ID = 0;
	transDirection = undetermined;

	start = 0;
	end = 0;

	type = frag_junction;
	altersite = 0;

	support = new double [SUPPORT_VECTOR_SIZE+1];
	for (int tmpCnt = 0; tmpCnt <= SUPPORT_VECTOR_SIZE; ++tmpCnt)
		support[tmpCnt] = 0.0;

	alter = NULL;
//	fivePalterTail = NULL;
	coverage = NULL;

	start_real = 0;
	end_real = 0;
}

fragment* fragment::clone()
{
	fragment *newFrag = new fragment;

	newFrag->frag_name = frag_name;
	newFrag->ID = ID;
	newFrag->chromosome_start = chromosome_start;
	newFrag->chromosome_end = chromosome_end;
	newFrag->start = start;
	newFrag->end = end;
	newFrag->type = type;

	for (int iLoop = 0; iLoop <= SUPPORT_VECTOR_SIZE; iLoop++)
		newFrag->support[iLoop] = support[iLoop];
	
	newFrag->alter = NULL;
//	newFrag->fivePalterTail = NULL;
	newFrag->coverage = NULL;
	newFrag->start_real = start_real;
	newFrag->end_real = end_real;

	return newFrag;
}

fragment::~fragment()
{
	alter_junction *delJunc;
	delJunc = alter;
	while (delJunc != NULL)
	{
		alter = delJunc->next;
		delete delJunc;
		delJunc = alter;
	}

	if (coverage != NULL)
		delete [] coverage;
	
	delete [] support;
}

spliceSite::spliceSite()
{
	position = 0;
	directionOut = true;
}

rangeJunction::rangeJunction()
{
	junc = NULL;
	next = NULL;
}

RangeJunctionList::RangeJunctionList()
{
	rangeLow = 0;
	rangeHigh = 0;
	transDirection = undetermined;
	list = NULL;
	listtail = NULL;
	nextList = NULL;
}

RangeJunctionList* RangeJunctionList::clone()
{
	//clone a same RangeJunctionList
	RangeJunctionList *resultList;
	resultList = new RangeJunctionList;

	resultList->rangeLow = rangeLow;
	resultList->rangeHigh = rangeHigh;
	resultList->transDirection = transDirection;
	resultList->list = NULL;
	resultList->listtail = NULL;
	resultList->nextList = NULL;

	rangeJunction *curList, *newList;
	curList = list;
	while (curList != NULL)
	{
		newList = new rangeJunction;
		newList->junc = curList->junc;
		newList->next = NULL;

		if (resultList->list == NULL)
		{
			resultList->list = newList;
			resultList->listtail = newList;
		}
		else
		{
			resultList->listtail->next = newList;
			resultList->listtail = newList;
		}

		curList = curList->next;
	}

	return resultList;
}

RangeJunctionList::~RangeJunctionList()
{
	nextList = NULL;

	rangeJunction *delList;

	while (list != NULL)
	{
		delList = list;
		list = delList->next;
		delete delList;
	}
}

GTedge::GTedge()
{
	linkedVertex = NULL;
	next = NULL;
}

alternative_path::alternative_path()
{
	path_start = 0;
	path_end = 0;
	whole_path_start = 0;
	whole_path_end = 0;
	transDirection = undetermined;
	support = new double [SUPPORT_VECTOR_SIZE];
	proportion = new double [SUPPORT_VECTOR_SIZE];
	for (int i = 0; i < SUPPORT_VECTOR_SIZE; i++)
	{
		support[i] = 0.0;
		proportion[i] = 1.0;
	}
	junctionNum = 0;
	exonNum = 0;
	pathVertex = NULL;
	next = NULL;
}

alternative_path::~alternative_path()
{
	delete [] support;
	delete [] proportion;
}

GTvertex::GTvertex()
{
	ID = 0;
	level = 0;
	rangeLow = 0;
	rangeHigh = 0;
	child = NULL;
	childType = 0;
	childNum = 0;
	junctionInRange = NULL;
	junctionNum = 0;
	exonNum = 0;
	prevSibling = NULL;
	nextSibling = NULL;

	alterSpliceSite = NULL;

	support = new double [SUPPORT_VECTOR_SIZE];
	proportion = new double [SUPPORT_VECTOR_SIZE];
	MSE_estimation = new double [SUPPORT_VECTOR_SIZE];
	min_path_support = new double [SUPPORT_VECTOR_SIZE];
	obs_support = new double [SUPPORT_VECTOR_SIZE];
	for (int i = 0; i < SUPPORT_VECTOR_SIZE; i++)
	{
		support[i] = 0.0;
		proportion[i] = 1.0;
		MSE_estimation[i] = 0.0;
		min_path_support[i] = 0.0;
		obs_support[i] = 0.0;
	}

	total_inflow = 0.0;
	total_outflow = 0.0;

	estimated = false;
	representative = NULL;
	estimate_exonNum = 0;

	ASMcategory = -1;
	ASMsupport_group1 = 0.0;
	ASMsupport_group2 = 0.0;

	major_alter_paths = NULL;
	major_alter_paths_num = 0;
}

GTvertex::~GTvertex()
{
// 	if (junctionInRange != NULL)
// 	{
// 		delete junctionInRange;
// 	}
// 
// 	delete representative;
// 
// 	//... delete children
	
	delete [] support;
	delete [] proportion;
	delete [] MSE_estimation;
	delete [] min_path_support;
	delete [] obs_support;
}

GenomeTree::GenomeTree()
{
	root = NULL;
}

GenomeTree::~GenomeTree()
{
	//
}

deletedSites::deletedSites()
{
	sites = 0;
	next = NULL;
}

/************************************************************************/
/* Vertex Stack for DFS                                                 */
/************************************************************************/
void stack_initial()
{
	stack_GTvertex.clear();
	stack_GTvertex.reserve(DEFAULT_JUNCTION_NUM);

	return;
}

bool stack_empty()
{
	return stack_GTvertex.empty();
}

void stack_push(GTvertex *x)
{
	if (stack_GTvertex.size() >= stack_GTvertex.capacity())
		stack_GTvertex.reserve(stack_GTvertex.size() + DEFAULT_JUNCTION_NUM);

	stack_GTvertex.push_back(x);

	return;
}

GTvertex* stack_pop()
{
	if (stack_empty() == true)
	{
		return NULL;
	} 
	else
	{
		GTvertex *vertex_pop = stack_GTvertex[stack_GTvertex.size()-1];
		stack_GTvertex.pop_back();
		return vertex_pop;
	}
}





/************************************************************************/
/* Construct GTree                                                      */
/************************************************************************/

void vertexForAlterSpliceSite(fragment *mergedJunction, GTvertex *fatherVertex)
{
	//build a representative GTvertex for the merged alternative splice site junctions
	//add the vertex to the fatherVertex
	if (mergedJunction->alter == NULL)
	{
		return;
	}

	GTvertex *newVertex, *resultVertex;
	newVertex = new GTvertex;

	alter_junction *cur_alter_junction;
	GTedge *newEdge;
	RangeJunctionList *newRangeJunctionList;
	rangeJunction *newRangeJunction;
	bool separable = false;
	long iLoop, jLoop;
	fragment **alterFragList;

	alterFragList = new fragment * [mergedJunction->alterFragCnt + 1];

	//build the junction list for the vertex
	newRangeJunctionList = new RangeJunctionList;
	newRangeJunctionList->rangeLow = mergedJunction->start;
	newRangeJunctionList->rangeHigh = mergedJunction->end;
	iLoop = 1;
	cur_alter_junction = mergedJunction->alter;
	while (cur_alter_junction != NULL)
	{
		if (cur_alter_junction->juncInfo->alter != NULL)
		{
			cout << "Warning: nested alternative splice site." << endl;
		}
		alterFragList[iLoop] = cur_alter_junction->juncInfo->clone();
		
		//reload the real start and end!
		if (alterFragList[iLoop]->start_real != 0)
		{
			alterFragList[iLoop]->start = alterFragList[iLoop]->start_real;
		}
		if (alterFragList[iLoop]->end_real != 0)
		{
			alterFragList[iLoop]->end = alterFragList[iLoop]->end_real;
		}

		iLoop++;
		cur_alter_junction = cur_alter_junction->next;
	}

	//selection sort
	fragment *curFragment, *tmpFragment;
	long minIndex;
	for (iLoop = 1; iLoop < mergedJunction->alterFragCnt; iLoop++)
	{
		minIndex = iLoop;
		for (jLoop = iLoop + 1; jLoop <= mergedJunction->alterFragCnt; jLoop++)
		{
			curFragment = alterFragList[jLoop];
			if (curFragment->start < alterFragList[minIndex]->start || curFragment->start == alterFragList[minIndex]->start && curFragment->end < alterFragList[minIndex]->end)
			{
				//curFragment is less than alterFragList[minIndex], replace
				minIndex = jLoop;
			}
		}
		//exchange selection
		if (minIndex != iLoop)
		{
			tmpFragment = alterFragList[iLoop];
			alterFragList[iLoop] = alterFragList[minIndex];
			alterFragList[minIndex] = tmpFragment;
		}
	}

	for (iLoop = mergedJunction->alterFragCnt; iLoop >= 1 ; iLoop--)
	{
		newRangeJunction = new rangeJunction;
		newRangeJunction->junc = alterFragList[iLoop];
		newRangeJunction->next = newRangeJunctionList->list;
		newRangeJunctionList->list = newRangeJunction;
	}

	delete [] alterFragList;


	//fill the vertex
	newVertex->level = fatherVertex->level; //same level as fatherVertex
	newVertex->rangeLow = mergedJunction->start;
	newVertex->rangeHigh = mergedJunction->end;
	if (fatherVertex->childType == 0)
	{
		//the merged junction is isolated
		newVertex->level = fatherVertex->level; //same level as fatherVertex
		newVertex->prevSibling = fatherVertex->prevSibling; //same siblings as the father
		newVertex->nextSibling = fatherVertex->nextSibling; //same siblings as the father
	} 
	else
	{
		//the merged junction is in a path (dependent path)
		newVertex->level = fatherVertex->level + 1; //same level as fatherVertex
		newVertex->prevSibling = NULL;
		newVertex->nextSibling = NULL;
	}
	newVertex->junctionInRange = newRangeJunctionList;
	resultVertex = newVertex;

	//construct the paths
	newRangeJunctionList = separateDepPath(newRangeJunctionList, separable);

	if (separable == true)
	{
		resultVertex->childType = 3;

		while (newRangeJunctionList != NULL)
		{
			newVertex = new GTvertex;
			newVertex->level = resultVertex->level + 1;
			newVertex->junctionInRange = newRangeJunctionList;
			newVertex->rangeLow = newRangeJunctionList->rangeLow;
			newVertex->rangeHigh = newRangeJunctionList->rangeHigh;
			newVertex->childType = 0;

			newRangeJunctionList = newRangeJunctionList->nextList;
			newVertex->junctionInRange->nextList = NULL;

			newEdge = new GTedge;
			newEdge->linkedVertex = newVertex;
			newEdge->next = resultVertex->child;
			resultVertex->child = newEdge;
		}
	} 
	else
	{
		resultVertex->junctionInRange = newRangeJunctionList;
		resultVertex->childType = 0;
	}

	//count the vertex
	countGTree(resultVertex);

	resultVertex->ASMcategory = alter_splice_site;

	newEdge = new GTedge;
	newEdge->linkedVertex = resultVertex;
	newEdge->next = fatherVertex->alterSpliceSite;
	fatherVertex->alterSpliceSite = newEdge;

	return;
}



bool calculateGeneMeanCoverage(RangeJunctionList *fragmentList, double *meanCoverageList)
{
	rangeJunction *curJunc;
	
	//count number of exons in the gene
	long exonCnt = 0, exonCntLoop;
	int sampleLoop;

	curJunc = fragmentList->list;
	while (curJunc != NULL)
	{
		if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
			++exonCnt;
		curJunc = curJunc->next;
	}

	if (exonCnt < 5)
		return false; //no need to do filtering
	
	//get expression array
	double **coverageArray = new double* [SUPPORT_VECTOR_SIZE];
	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		coverageArray[sampleLoop] = new double [exonCnt + 1];
	
		curJunc = fragmentList->list;
		exonCntLoop = 0;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
				coverageArray[sampleLoop][++exonCntLoop] = curJunc->junc->support[sampleLoop];
			curJunc = curJunc->next;
		}
	}
	long *quicksortStack = new long [2*exonCnt + 10];

	//calculate mean coverage in each sample
	long index_q25, index_q75;
	double total_coverage;
	index_q25 = long(floor(0.25 * (exonCnt - 1) + 1));
	index_q75 = long(floor(0.75 * (exonCnt - 1) + 1));

	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		quicksort(coverageArray[sampleLoop], exonCnt, quicksortStack);

		total_coverage = 0.0;
		for (exonCntLoop = index_q25 + 1; exonCntLoop <= index_q75; ++exonCntLoop)
			total_coverage += coverageArray[sampleLoop][exonCntLoop];

		meanCoverageList[sampleLoop] = total_coverage / (index_q75 - index_q25);
	}

	//clean up
	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		delete [] coverageArray[sampleLoop];
	}
	delete [] coverageArray;
	delete [] quicksortStack;

	return true;
}


void filterGeneLowCoverageExon(GTvertex *targetVertex)
{
	double *meanCoverageList = new double [SUPPORT_VECTOR_SIZE];
	rangeJunction *curJunc, *prevJunc;
	long iLoop;
	int thresh_lowCovSampleCnt = int(ceil(SUPPORT_VECTOR_SIZE * (1 - coverageThreshold_GeneExon))), lowCovSampleCnt, sampleLoop;

	//find mean coverage of 25% to 75% percentile expression in this gene
	if (calculateGeneMeanCoverage(targetVertex->junctionInRange, meanCoverageList) == false)
		return;

	for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
	{
		targetVertex->proportion[sampleLoop] = meanCoverageList[sampleLoop];
	}

	//do not run this for now
	delete [] meanCoverageList;
	return;
	
	
	//filter exons with coverage less than threshold * mean coverage
	curJunc = targetVertex->junctionInRange->list; prevJunc = NULL;
	while (curJunc != NULL)
	{
		if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
		{
			lowCovSampleCnt = 0;
			for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
			{
				if (curJunc->junc->support[sampleLoop] < coverageThreshold_GeneExon * meanCoverageList[sampleLoop])
					++lowCovSampleCnt;
			}

			if (lowCovSampleCnt >= thresh_lowCovSampleCnt)
			{
				//filter
				if (prevJunc != NULL)
				{
					prevJunc->next = curJunc->next;
					delete curJunc;
					curJunc = prevJunc->next;
				}
				else
				{
					targetVertex->junctionInRange->list = curJunc->next;
					delete curJunc;
					curJunc = targetVertex->junctionInRange->list;
				}
			}
			else
			{
				prevJunc = curJunc;
				curJunc = curJunc->next;
			}
		}
		else
		{
			prevJunc = curJunc;
			curJunc = curJunc->next;
		}
	}

	//filter junctions connected to nothing
	sortList_Junction_Num = 0;
	curJunc = targetVertex->junctionInRange->list;
	while (curJunc != NULL)
	{
		if (sortList_Junction_Num+1 >= sortList_Junction.size())
		{
			sortList_Junction.resize(sortList_Junction.size() + DEFAULT_JUNCTION_NUM, NULL);
			sortKey_Junction.resize(sortKey_Junction.size() + DEFAULT_JUNCTION_NUM, 0);
			mergeSort_Larray.resize(mergeSort_Larray.size() + DEFAULT_JUNCTION_NUM, 0);
			mergeSort_Rarray.resize(mergeSort_Rarray.size() + DEFAULT_JUNCTION_NUM, 0);
			mergeSort_LorderedList.resize(mergeSort_LorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
			mergeSort_RorderedList.resize(mergeSort_RorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
		}
		sortList_Junction[++sortList_Junction_Num] = curJunc;
		curJunc = curJunc->next;
	}
	filterFalseFragments(true);

	//reform the fragment list
	for (iLoop = 1; iLoop <= sortList_Junction_Num; ++iLoop)
	{
		sortKey_Junction[iLoop] = sortList_Junction[iLoop]->junc->end;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);

	for (iLoop = 1; iLoop <= sortList_Junction_Num; ++iLoop)
	{
		sortKey_Junction[iLoop] = sortList_Junction[iLoop]->junc->start;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);

	targetVertex->junctionInRange->list = NULL;
	for (iLoop = sortList_Junction_Num; iLoop >= 1; --iLoop)
	{
		sortList_Junction[iLoop]->next = targetVertex->junctionInRange->list;
		targetVertex->junctionInRange->list = sortList_Junction[iLoop];
	}	

	return;
}

void constructGTree(RangeJunctionList* origList)
{
	//construct GTree
	//input: sorted range fragment list
	bool separable_region, separable_path, separable_depPath;
	RangeJunctionList *curList, *curPathList, *gtfPathList;
	GTvertex *newVertex;
	GTedge *newEdge;

	stack_initial();

	GTvertex *curVertex, *lastVertex;
	curVertex = new GTvertex;
	curVertex->junctionInRange = origList;
	curVertex->level = 0;
	curVertex->rangeLow = CHROMOSOME_START;
	curVertex->rangeHigh = CHROMOSOME_END;
	curVertex->child = NULL;
	curVertex->childType = 1;

	gTree->root = curVertex;

	//separate genes
	separable_region = true;
	curList = NULL;

	curList = separateGene(curVertex->junctionInRange, separable_region);

	if (separable_region == true)
	{
		while (curList != NULL)
		{
			//(1) if two genes are overlapping, e.g., positive strand and negative strand genes on the same locus
			//(2) if there are some standalone stuff, like false region, microRNA, etc.   they are not paths, but may encode some important information
			//therefore, separate them first.
			curPathList = curList;
			curList = curList->nextList;
			curPathList->nextList = NULL;
			separable_path = true;

			curPathList = separateIndepPath(curPathList, separable_path);
			
			if (separable_path == true)
			{
				while (curPathList != NULL)
				{
					newVertex = new GTvertex;
					newVertex->level = curVertex->level; //genes, same level as the chromosome, 0
					newVertex->junctionInRange = curPathList;
					newVertex->rangeLow = curPathList->rangeLow;
					newVertex->rangeHigh = curPathList->rangeHigh;
					newVertex->child = NULL;
					newVertex->childType = 1;

					curPathList = curPathList->nextList;
					newVertex->junctionInRange->nextList = NULL;

					newEdge = new GTedge;
					newEdge->linkedVertex = newVertex;
					newEdge->next = curVertex->child;
					curVertex->child = newEdge;

					filterGeneLowCoverageExon(newVertex);
					stack_push(newVertex);
				}
			}
			else
			{
				newVertex = new GTvertex;
				newVertex->level = curVertex->level; //genes, same level as the chromosome, 0
				newVertex->junctionInRange = curPathList;
				newVertex->rangeLow = curPathList->rangeLow;
				newVertex->rangeHigh = curPathList->rangeHigh;
				newVertex->child = NULL;
				newVertex->childType = 1;
				
				newEdge = new GTedge;
				newEdge->linkedVertex = newVertex;
				newEdge->next = curVertex->child;
				curVertex->child = newEdge;

				filterGeneLowCoverageExon(newVertex);
				stack_push(newVertex);
			}			
		}
	} 
	else
	{
		curVertex->junctionInRange = curList;
		curVertex->childType = 1;

		//(1) if two genes are overlapping, e.g., positive strand and negative strand genes on the same locus
		//(2) if there are some standalone stuff, like false region, microRNA, etc.   they are not paths, but may encode some important information
		//therefore, separate them first.
		separable_path = true;
		curList = NULL;

		curList = separateIndepPath(curVertex->junctionInRange, separable_path);

		if (separable_path == true)
		{
			while (curList != NULL)
			{
				newVertex = new GTvertex;
				newVertex->level = curVertex->level; //genes, same level as the chromosome, 0
				newVertex->junctionInRange = curList;
				newVertex->rangeLow = curList->rangeLow;
				newVertex->rangeHigh = curList->rangeHigh;
				newVertex->child = NULL;
				newVertex->childType = 1;

				curList = curList->nextList;
				newVertex->junctionInRange->nextList = NULL;

				newEdge = new GTedge;
				newEdge->linkedVertex = newVertex;
				newEdge->next = curVertex->child;
				curVertex->child = newEdge;
				
				filterGeneLowCoverageExon(newVertex);
				stack_push(newVertex);
			}
		}
		else
		{
			curVertex->junctionInRange = curList;
			filterGeneLowCoverageExon(curVertex);
			stack_push(curVertex);
		}
	}


	//separate ASMs
	curVertex = stack_pop();
	while (curVertex != NULL)
	{
		separable_region = true;
		separable_path = true;
		separable_depPath = true;

		lastVertex = NULL;
		curList = NULL;

		if (curVertex->level < 1)
		{
			curVertex->ID = ++GENEcount;
		}
		
		if (curVertex->childType == 1)
		{
			//independent regions
//			cout << curVertex->rangeLow << "\t" << curVertex->rangeHigh << endl;
//			if (curVertex->rangeLow == 50569593)
//				cout << "catch --";
			curList = separateIndepRegion(curVertex->junctionInRange, separable_region);
//			cout << "done" << endl;

			if (separable_region == true)
			{
				while (curList != NULL)
				{
					newVertex = new GTvertex;
					newVertex->level = curVertex->level + 1;
					newVertex->junctionInRange = curList;
					newVertex->rangeLow = curList->rangeLow;
					newVertex->rangeHigh = curList->rangeHigh;
					newVertex->childType = 2;

					curList = curList->nextList;
					newVertex->junctionInRange->nextList = NULL;

					newEdge = new GTedge;
					newEdge->linkedVertex = newVertex;
					newEdge->next = curVertex->child;
					curVertex->child = newEdge;

					if (lastVertex != NULL)
					{
						lastVertex->prevSibling = newVertex;
						newVertex->nextSibling = lastVertex;
					}
					lastVertex = newVertex;

					stack_push(newVertex);
				}
			} 
			else
			{
				curVertex->junctionInRange = curList;
				curVertex->childType = 2;

				stack_push(curVertex);
			}
		} 
		else if (curVertex->childType == 2)
		{
			//independent paths

			/////////////////////////////////////////////////////////////////////////////
			//output gtf tracks for level 1 ASMs 
			//the 2nd condition check whether there are multiple fragments, i.e., an ASM
			if (curVertex->level == 1 && curVertex->junctionInRange->list->next != NULL)
			{		
				curVertex->ID = ++ASMcount;
				gtfPathList = separateDepPath(curVertex->junctionInRange, separable_depPath);
				if (separable_depPath == true)
				{
					output_ASMpath_gtf(curVertex, gtfPathList);
					delete gtfPathList;
				}
				else
				{
					output_ASMpath_gtf(curVertex, gtfPathList);
					curVertex->junctionInRange = gtfPathList;
				}
			}
			/////////////////////////////////////////////////////////////////////////////



			curList = separateIndepPath(curVertex->junctionInRange, separable_path);

			if (separable_path == true)
			{
				while (curList != NULL)
				{
					newVertex = new GTvertex;
					newVertex->level = curVertex->level + 1;
					newVertex->junctionInRange = curList;
					newVertex->rangeLow = curList->rangeLow;
					newVertex->rangeHigh = curList->rangeHigh;
					newVertex->childType = 1;

					curList = curList->nextList;
					newVertex->junctionInRange->nextList = NULL;

					newEdge = new GTedge;
					newEdge->linkedVertex = newVertex;
					newEdge->next = curVertex->child;
					curVertex->child = newEdge;

// 					if (lastVertex != NULL)
// 					{
// 						lastVertex->prevSibling = newVertex;
// 						newVertex->nextSibling = lastVertex;
// 					}
// 					lastVertex = newVertex;

					stack_push(newVertex);
				}
			} 
			else
			{
				//reach a leaf
				curVertex->junctionInRange = curList;
				curList = separateDepPath(curList, separable_depPath);

				if (separable_depPath == true)
				{
					curVertex->childType = 3;

					while (curList != NULL)
					{
						newVertex = new GTvertex;
						newVertex->level = curVertex->level + 1;
						newVertex->junctionInRange = curList;
						newVertex->rangeLow = curList->rangeLow;
						newVertex->rangeHigh = curList->rangeHigh;
						newVertex->childType = 0;

						curList = curList->nextList;
						newVertex->junctionInRange->nextList = NULL;

						newEdge = new GTedge;
						newEdge->linkedVertex = newVertex;
						newEdge->next = curVertex->child;
						curVertex->child = newEdge;

// 						if (lastVertex != NULL)
// 						{
// 							lastVertex->prevSibling = newVertex;
// 							newVertex->nextSibling = lastVertex;
// 						}
// 						lastVertex = newVertex;
					}
				} 
				else
				{
					curVertex->junctionInRange = curList;
					curVertex->childType = 0;
				}
			}
		}

		curVertex = stack_pop();
	}

	return;
}





/************************************************************************/
/* Genome Statistics                                                    */
/************************************************************************/

// fragment* makeRepresentative(GTvertex *rootVertex)
// {
// 	//make representative for rootVertex
// 	//the representative will consist of a junction (if needed), an exon (length = all exonic region, support = all sum), and a junction (if needed)
// 	//currently only use an exon to represent
// 	GTvertex *curVertex;
// 	GTedge *curEdge;
// 	rangeJunction *curJunc;
// 	double exonicLengthSum = 0.0;
// 	double exonicSupportSum[SUPPORT_VECTOR_SIZE];
// 	int iLoop;
// 
// 	for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 	{
// 		exonicSupportSum[iLoop] = 0.0;
// 	}
// 
// 	if (rootVertex->childType == 3)
// 	{
// 		//count exonic length and support from fragment list
// 		curJunc = rootVertex->junctionInRange->list;
// 		while (curJunc != NULL)
// 		{
// 			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
// 			{
// 				exonicLengthSum += fabs(curJunc->junc->end - curJunc->junc->start);
// 				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 				{
// 					exonicSupportSum[iLoop] += curJunc->junc->support[iLoop]; //use sum of observational count, or estimated count??
// 				}
// 			}
// 			curJunc = curJunc->next;
// 		}
// 	}
// 	else if (rootVertex->childType == 2 || rootVertex->childType == 1)
// 	{
// 		//count exonic length and support from all children
// 		curEdge = rootVertex->child;
// 		while (curEdge != NULL)
// 		{
// 			curVertex = curEdge->linkedVertex;
// 			
// 			if (curVertex->estimated == false)
// 			{
// 				//leaf node
// 				curJunc = curVertex->junctionInRange->list;
// 				while (curJunc != NULL)
// 				{
// 					if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
// 					{
// 						exonicLengthSum += fabs(curJunc->junc->end - curJunc->junc->start);
// 						for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 						{
// 							exonicSupportSum[iLoop] += curJunc->junc->support[iLoop]; //use sum of observational count, or estimated count??
// 						}
// 					}
// 					curJunc = curJunc->next;
// 				}
// 			} 
// 			else
// 			{
// 				exonicLengthSum += abs(curVertex->representative->end - curVertex->representative->start);
// 				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 				{
// 					exonicSupportSum[iLoop] += curVertex->representative->support[iLoop];
// 				}
// 			}
// 
// 			curEdge = curEdge->next;
// 		}
// 	}
// 
// 	fragment *newFragment;
// 	newFragment = new fragment;
// 	newFragment->ID = ++fragmentID_Cnt;
// 	newFragment->end = newFragment->start + exonicLengthSum;
// 	newFragment->type = frag_exon;
// 	for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 	{
// 		newFragment->support[iLoop] = exonicSupportSum[iLoop];
// 		rootVertex->support[iLoop] = exonicSupportSum[iLoop];
// 	}
// 
// 	return newFragment;
// }

fragment* makeRepresentative(GTvertex *rootVertex)
{
	//make representative for rootVertex
	//the representative will consist of a junction (if needed), an exon (length = all exonic region, support = all sum), and a junction (if needed)
	//currently only use an exon to represent
	GTvertex *curVertex;
	GTedge *curEdge;
	rangeJunction *curJunc;
	long exonicLength = 0;
	double *exonicSupportSum = new double [SUPPORT_VECTOR_SIZE];
	int exonicChildNum = 0;
	int iLoop;

	for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
	{
		exonicSupportSum[iLoop] = 0.0;
	}

	if (rootVertex->childType == 3 || rootVertex->childType == 2)
	{
		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;

			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
			{
				exonicSupportSum[iLoop] += curVertex->support[iLoop];
			}

			curEdge = curEdge->next;
		}
	}
	else if (rootVertex->childType == 1)
	{
		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;

			if (curVertex->estimated == false)
			{
				curJunc = curVertex->junctionInRange->list;
				while (curJunc != NULL)
				{
					if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
					{
						for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
						{
							exonicSupportSum[iLoop] += curJunc->junc->support[iLoop];
						}
						exonicChildNum++;
					}
					curJunc = curJunc->next;
				}
			} 
			else
			{
				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
				{
					exonicSupportSum[iLoop] += curVertex->representative->support[iLoop];
				}
				exonicChildNum++;
			}

			if (curVertex->childNum != 0 || curVertex->exonNum != 0)
				exonicLength += curVertex->rangeHigh - curVertex->rangeLow + 1;
			
			curEdge = curEdge->next;
		}

		for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
		{
			if (exonicChildNum > 0)
			{
				exonicSupportSum[iLoop] = exonicSupportSum[iLoop] / exonicChildNum;
			} 
			else
			{
				exonicSupportSum[iLoop] = 0.0;
			}
		}
	}


	fragment *newFragment;
	newFragment = new fragment;
	newFragment->ID = ++fragmentID_Cnt;
	newFragment->start = 0;
	newFragment->end = 200;
	newFragment->type = frag_exon;
	for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
	{
		newFragment->support[iLoop] = exonicSupportSum[iLoop];
		rootVertex->support[iLoop] = exonicSupportSum[iLoop];
	}

	delete [] exonicSupportSum;

	return newFragment;
}


// int alterStart_onetest(GTvertex *vertex1, GTvertex *vertex2, GTvertex *vertex3)
// {
// 	//return whether there is an alternative start/end
// 	GTedge *curEdge;
// 	rangeJunction *curJunc, *mergedJunc, *newJunc;
// 	RangeJunctionList *curJuncList;
// 	double *LL_null, *LL_alterstart, *LL_alterend, MSE;
// 	LL_null = new double [SUPPORT_VECTOR_SIZE];
// 	LL_alterstart = new double [SUPPORT_VECTOR_SIZE];
// 	LL_alterend = new double [SUPPORT_VECTOR_SIZE];
// 	int iLoop;
// 
// 	//null model -- one transcript, no alternative start
// 	GTvertex *nullVertex, *childVertex1;
// 	
// 	nullVertex = new GTvertex;
// 	nullVertex->level = vertex1->level;
// 	nullVertex->rangeLow = vertex1->rangeLow;
// 	nullVertex->rangeHigh = vertex3->rangeHigh;
// 	nullVertex->childType = 2;
// 	nullVertex->childNum = 1;
// 	nullVertex->nextSibling = NULL;
// 	nullVertex->prevSibling = NULL;
// 
// 	mergedJunc = NULL;
// 	curJunc = vertex3->junctionInRange->list;
// 	while (curJunc != NULL)
// 	{
// 		newJunc = new rangeJunction;
// 		newJunc->junc = curJunc->junc;
// 		newJunc->next = mergedJunc;
// 		mergedJunc = newJunc;
// 		curJunc = curJunc->next;
// 	}
// 	curJunc = vertex2->junctionInRange->list;
// 	while (curJunc != NULL)
// 	{
// 		newJunc = new rangeJunction;
// 		newJunc->junc = curJunc->junc;
// 		newJunc->next = mergedJunc;
// 		mergedJunc = newJunc;
// 		curJunc = curJunc->next;
// 	}
// 	curJunc = vertex1->junctionInRange->list;
// 	while (curJunc != NULL)
// 	{
// 		newJunc = new rangeJunction;
// 		newJunc->junc = curJunc->junc;
// 		newJunc->next = mergedJunc;
// 		mergedJunc = newJunc;
// 		curJunc = curJunc->next;
// 	}
// 	curJuncList = new RangeJunctionList;
// 	curJuncList->rangeLow = nullVertex->rangeLow;
// 	curJuncList->rangeHigh = nullVertex->rangeHigh;
// 	curJuncList->list = mergedJunc;
// 	nullVertex->junctionInRange = curJuncList;
// 	nullVertex->junctionNum = vertex1->junctionNum + vertex2->junctionNum + vertex3->junctionNum;
// 	nullVertex->exonNum = vertex1->exonNum + vertex2->exonNum + vertex3->exonNum;
// 
// 
// 	childVertex1 = new GTvertex;
// 	childVertex1->level = nullVertex->level + 1;
// 	childVertex1->rangeLow = nullVertex->rangeLow;
// 	childVertex1->rangeHigh = nullVertex->rangeHigh;
// 	childVertex1->childType = 0;
// 	childVertex1->childNum = 0;
// 	childVertex1->nextSibling = NULL;
// 	childVertex1->prevSibling = NULL;
// 
// 	childVertex1->junctionInRange = nullVertex->junctionInRange->clone();
// 	childVertex1->junctionNum = nullVertex->junctionNum;
// 	childVertex1->exonNum = nullVertex->exonNum;
// 
// 	curEdge = new GTedge;
// 	curEdge->linkedVertex = childVertex1;
// 	curEdge->next = NULL;
// 	nullVertex->child = curEdge;
// 
// 	//calculate loglikelihood
// 	for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 	{
// 		LL_null[iLoop] = abundance_estimation(nullVertex, iLoop, MSE);
// 	}
// 
// 
// 	//alternative model -- two transcript, alternative start
// 	GTvertex *alterVertex_start, *childVertex2;
// 
// 	alterVertex_start = new GTvertex;
// 	alterVertex_start->level = vertex1->level;
// 	alterVertex_start->rangeLow = vertex1->rangeLow;
// 	alterVertex_start->rangeHigh = vertex3->rangeHigh;
// 	alterVertex_start->childType = 3;
// 	alterVertex_start->childNum = 2;
// 	alterVertex_start->nextSibling = NULL;
// 	alterVertex_start->prevSibling = NULL;
// 
// 	alterVertex_start->junctionInRange = nullVertex->junctionInRange->clone();
// 	alterVertex_start->junctionNum = vertex1->junctionNum + vertex2->junctionNum + vertex3->junctionNum;
// 	alterVertex_start->exonNum = vertex1->exonNum + vertex2->exonNum + vertex3->exonNum;
// 
// 
// 	childVertex1 = new GTvertex;
// 	childVertex1->level = alterVertex_start->level + 1;
// 	childVertex1->rangeLow = alterVertex_start->rangeLow;
// 	childVertex1->rangeHigh = alterVertex_start->rangeHigh;
// 	childVertex1->childType = 0;
// 	childVertex1->childNum = 0;
// 	childVertex1->nextSibling = NULL;
// 	childVertex1->prevSibling = NULL;
// 
// 	childVertex1->junctionInRange = alterVertex_start->junctionInRange->clone();
// 	childVertex1->junctionNum = alterVertex_start->junctionNum;
// 	childVertex1->exonNum = alterVertex_start->exonNum;
// 
// 	childVertex2 = new GTvertex;
// 	childVertex2->level = alterVertex_start->level + 1;
// 	childVertex2->rangeLow = vertex3->rangeLow;
// 	childVertex2->rangeHigh = vertex3->rangeHigh;
// 	childVertex2->childType = 0;
// 	childVertex2->childNum = 0;
// 	childVertex2->nextSibling = NULL;
// 	childVertex2->prevSibling = NULL;
// 
// 	childVertex2->junctionInRange = vertex3->junctionInRange->clone();
// 	childVertex2->junctionNum = vertex3->junctionNum;
// 	childVertex2->exonNum = vertex3->exonNum;
// 
// 	curEdge = new GTedge;
// 	curEdge->linkedVertex = childVertex2;
// 	curEdge->next = NULL;
// 	alterVertex_start->child = curEdge;
// 	curEdge = new GTedge;
// 	curEdge->linkedVertex = childVertex1;
// 	curEdge->next = alterVertex_start->child;
// 	alterVertex_start->child = curEdge;
// 
// 	//calculate loglikelihood
// 	for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 	{
// 		LL_alterstart[iLoop] = abundance_estimation(alterVertex_start, iLoop, MSE);
// 	}
// 
// 
// 	//alternative model -- two transcript, alternative end
// 	GTvertex *alterVertex_end;
// 
// 	alterVertex_end = new GTvertex;
// 	alterVertex_end->level = vertex1->level;
// 	alterVertex_end->rangeLow = vertex1->rangeLow;
// 	alterVertex_end->rangeHigh = vertex3->rangeHigh;
// 	alterVertex_end->childType = 3;
// 	alterVertex_end->childNum = 2;
// 	alterVertex_end->nextSibling = NULL;
// 	alterVertex_end->prevSibling = NULL;
// 
// 	alterVertex_end->junctionInRange = nullVertex->junctionInRange->clone();
// 	alterVertex_end->junctionNum = vertex1->junctionNum + vertex2->junctionNum + vertex3->junctionNum;
// 	alterVertex_end->exonNum = vertex1->exonNum + vertex2->exonNum + vertex3->exonNum;
// 
// 
// 	childVertex1 = new GTvertex;
// 	childVertex1->level = alterVertex_end->level + 1;
// 	childVertex1->rangeLow = alterVertex_end->rangeLow;
// 	childVertex1->rangeHigh = alterVertex_end->rangeHigh;
// 	childVertex1->childType = 0;
// 	childVertex1->childNum = 0;
// 	childVertex1->nextSibling = NULL;
// 	childVertex1->prevSibling = NULL;
// 
// 	childVertex1->junctionInRange = alterVertex_end->junctionInRange->clone();
// 	childVertex1->junctionNum = alterVertex_end->junctionNum;
// 	childVertex1->exonNum = alterVertex_end->exonNum;
// 
// 	childVertex2 = new GTvertex;
// 	childVertex2->level = alterVertex_end->level + 1;
// 	childVertex2->rangeLow = vertex1->rangeLow;
// 	childVertex2->rangeHigh = vertex1->rangeHigh;
// 	childVertex2->childType = 0;
// 	childVertex2->childNum = 0;
// 	childVertex2->nextSibling = NULL;
// 	childVertex2->prevSibling = NULL;
// 
// 	childVertex2->junctionInRange = vertex1->junctionInRange->clone();
// 	childVertex2->junctionNum = vertex1->junctionNum;
// 	childVertex2->exonNum = vertex1->exonNum;
// 
// 	curEdge = new GTedge;
// 	curEdge->linkedVertex = childVertex2;
// 	curEdge->next = NULL;
// 	alterVertex_end->child = curEdge;
// 	curEdge = new GTedge;
// 	curEdge->linkedVertex = childVertex1;
// 	curEdge->next = alterVertex_end->child;
// 	alterVertex_end->child = curEdge;
// 
// 	//calculate loglikelihood
// 	for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 	{
// 		LL_alterend[iLoop] = abundance_estimation(alterVertex_end, iLoop, MSE);
// 	}
// 
// 
// 	//now the test...
// 	int diff_alterstart = 0;
// 	for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 	{
// 		if (LL_alterstart[iLoop] > LL_null[iLoop] && chisquare(2 * (LL_alterstart[iLoop] - LL_null[iLoop])) < 0.01)
// 		{
// 			//reject null and accept alternative
// 			diff_alterstart++;
// 		}
// 	}
// 
// 	if (diff_alterstart >= 1 && diff_alterstart >= SUPPORT_VECTOR_SIZE / 2 - 1)
// 	{
// 		//make the alternative
// 		delete vertex1->junctionInRange;
// 		
// 		//replace vertex1 with alterVertex_start
// 		vertex1->rangeLow = alterVertex_start->rangeLow;
// 		vertex1->rangeHigh = alterVertex_start->rangeHigh;
// 		vertex1->child = alterVertex_start->child;
// 		vertex1->childType = alterVertex_start->childType;
// 		vertex1->childNum = alterVertex_start->childNum;
// 		vertex1->junctionInRange = alterVertex_start->junctionInRange;
// 		vertex1->exonNum = alterVertex_start->exonNum;
// 		vertex1->junctionNum = alterVertex_start->junctionNum;
// 		vertex1->nextSibling = vertex3->nextSibling;
// 		if (vertex1->nextSibling != NULL)
// 		{
// 			vertex1->nextSibling->prevSibling = vertex1;
// 		}
// 
// 		delete vertex2->junctionInRange;
// 		delete vertex2;
// 		delete vertex3->junctionInRange;
// 		delete vertex3;
// 		delete alterVertex_start;
// 		delete nullVertex->junctionInRange;
// 		delete nullVertex->child->linkedVertex->junctionInRange;
// 		delete nullVertex->child->linkedVertex;
// 		delete nullVertex->child;
// 		delete nullVertex;
// 		delete alterVertex_end->junctionInRange;
// 		curEdge = alterVertex_end->child;
// 		while (curEdge != NULL)
// 		{
// 			delete curEdge->linkedVertex->junctionInRange;
// 			delete curEdge->linkedVertex;
// 			alterVertex_end->child = curEdge->next;
// 			delete curEdge;
// 			curEdge = alterVertex_end->child;
// 		}
// 		delete alterVertex_end;
// 
// 		return 1;
// 	}
// 
// 
// 	int diff_alterend = 0;
// 	for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 	{
// 		if (LL_alterend[iLoop] > LL_null[iLoop] && chisquare(2 * (LL_alterend[iLoop] - LL_null[iLoop])) < 0.01)
// 		{
// 			//reject null and accept alternative
// 			diff_alterend++;
// 		}
// 	}
// 
// 	if (diff_alterend >= 1 && diff_alterend >= SUPPORT_VECTOR_SIZE / 2 - 1)
// 	{
// 		//make the alternative
// 		delete vertex1->junctionInRange;
// 
// 		//replace vertex1 with alterVertex_end
// 		vertex1->rangeLow = alterVertex_end->rangeLow;
// 		vertex1->rangeHigh = alterVertex_end->rangeHigh;
// 		vertex1->child = alterVertex_end->child;
// 		vertex1->childType = alterVertex_end->childType;
// 		vertex1->childNum = alterVertex_end->childNum;
// 		vertex1->junctionInRange = alterVertex_end->junctionInRange;
// 		vertex1->exonNum = alterVertex_end->exonNum;
// 		vertex1->junctionNum = alterVertex_end->junctionNum;
// 		vertex1->nextSibling = vertex3->nextSibling;
// 		if (vertex1->nextSibling != NULL)
// 		{
// 			vertex1->nextSibling->prevSibling = vertex1;
// 		}
// 
// 		delete vertex2->junctionInRange;
// 		delete vertex2;
// 		delete vertex3->junctionInRange;
// 		delete vertex3;
// 		delete alterVertex_end;
// 		delete nullVertex->junctionInRange;
// 		delete nullVertex->child->linkedVertex->junctionInRange;
// 		delete nullVertex->child->linkedVertex;
// 		delete nullVertex->child;
// 		delete nullVertex;
// 		delete alterVertex_start->junctionInRange;
// 		curEdge = alterVertex_start->child;
// 		while (curEdge != NULL)
// 		{
// 			delete curEdge->linkedVertex->junctionInRange;
// 			delete curEdge->linkedVertex;
// 			alterVertex_start->child = curEdge->next;
// 			delete curEdge;
// 			curEdge = alterVertex_start->child;
// 		}
// 		delete alterVertex_start;
// 
// 		return 2;
// 	}
// 
// 
// 	delete nullVertex->junctionInRange;
// 	delete nullVertex->child->linkedVertex->junctionInRange;
// 	delete nullVertex->child->linkedVertex;
// 	delete nullVertex->child;
// 	delete nullVertex;
// 	delete alterVertex_start->junctionInRange;
// 	curEdge = alterVertex_start->child;
// 	while (curEdge != NULL)
// 	{
// 		delete curEdge->linkedVertex->junctionInRange;
// 		delete curEdge->linkedVertex;
// 		alterVertex_start->child = curEdge->next;
// 		delete curEdge;
// 		curEdge = alterVertex_start->child;
// 	}
// 	delete alterVertex_start;
// 	delete alterVertex_end->junctionInRange;
// 	curEdge = alterVertex_end->child;
// 	while (curEdge != NULL)
// 	{
// 		delete curEdge->linkedVertex->junctionInRange;
// 		delete curEdge->linkedVertex;
// 		alterVertex_end->child = curEdge->next;
// 		delete curEdge;
// 		curEdge = alterVertex_end->child;
// 	}
// 	delete alterVertex_end;
// 
// 	return 0;
// }
// 
// void alterStart(GTvertex *rootVertex)
// {
// 	GTedge *curEdge1, *curEdge2, *curEdge3;
// 	long alterStart_cnt = 0, alterEnd_cnt = 0; 
// 	int result;
// 
// 	if (rootVertex->child == NULL || rootVertex->child->next == NULL || rootVertex->child->next->next == NULL)
// 	{
// 		return;
// 	}
// 
// 	curEdge1 = rootVertex->child;
// 	curEdge2 = curEdge1->next;
// 	curEdge3 = curEdge2->next;
// 	while (curEdge3 != NULL)
// 	{
// 		result = 0;
// 		if (curEdge1->linkedVertex->childType == 0 && curEdge2->linkedVertex->childType == 0 && curEdge3->linkedVertex->childType == 0 && 
// 			curEdge1->linkedVertex->exonNum == 1 && curEdge1->linkedVertex->junctionNum == 0 && curEdge2->linkedVertex->exonNum == 0 && curEdge2->linkedVertex->junctionNum == 1 
// 			&& curEdge3->linkedVertex->exonNum == 1 && curEdge3->linkedVertex->junctionNum == 0
// 			&& curEdge1->linkedVertex->junctionInRange->list->junc->type == frag_exon
// 			&& curEdge2->linkedVertex->junctionInRange->list->junc->type == frag_junction && curEdge3->linkedVertex->junctionInRange->list->junc->type == frag_exon
// 			&& curEdge1->linkedVertex->rangeHigh == curEdge2->linkedVertex->rangeLow && curEdge2->linkedVertex->rangeHigh == curEdge3->linkedVertex->rangeLow)
// 		{
// 			result = alterStart_onetest(curEdge1->linkedVertex, curEdge2->linkedVertex, curEdge3->linkedVertex);
// 		}
// 
// 		if (result == 1)
// 		{
// 			curEdge1->linkedVertex->ASMcategory = alter_start;
// 			alterStartEnd << "alterStart\t" << curEdge1->linkedVertex->rangeLow << "\t" << curEdge1->linkedVertex->rangeHigh << endl;
// 			alterStart_cnt++;
// 		}
// 		if (result == 2)
// 		{
// 			curEdge1->linkedVertex->ASMcategory = alter_end;
// 			alterStartEnd << "alterEnd\t" << curEdge1->linkedVertex->rangeLow << "\t" << curEdge1->linkedVertex->rangeHigh << endl;
// 			alterEnd_cnt++;
// 		}
// 		
// 		if (result == 0)
// 		{
// 			curEdge1 = curEdge1->next;
// 			curEdge2 = curEdge2->next;
// 			curEdge3 = curEdge3->next;
// 		}
// 		else
// 		{
// 			curEdge1->next = curEdge3->next;
// 			delete curEdge2;
// 			delete curEdge3;
// 			curEdge1 = curEdge1->next;
// 			if (curEdge1 != NULL && curEdge1->next != NULL)
// 			{
// 				curEdge2 = curEdge1->next;
// 				curEdge3 = curEdge2->next;
// 			}
// 			else
// 			{
// 				curEdge2 = NULL;
// 				curEdge3 = NULL;
// 			}
// 		}
// 	}
// 
// 	if (alterStart_cnt > 0 || alterEnd_cnt > 0)
// 		cout << "*** alter start  " << alterStart_cnt << "    alter end  " << alterEnd_cnt << "***\t";
// 	
// 	return;
// }

bool preCountGTree(GTvertex *rootVertex)
{
	//collect basic counts on GTree for further statistics
	//collect child count and fragment count for every vertex
	//return true for leaves
	GTvertex *curVertex;
	GTedge *curEdge;
	rangeJunction *curJunc;
	double *fragSupportSum = new double [SUPPORT_VECTOR_SIZE], *exonSupportSum = new double [SUPPORT_VECTOR_SIZE]; 
	int juncCnt, exonCnt, fragCnt, tmp, iLoop;

	if (rootVertex->child == NULL)
	{
		//leaf
		juncCnt = 0;
		exonCnt = 0;
		fragCnt = 0;
		for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
		{
			fragSupportSum[iLoop] = 0;
			exonSupportSum[iLoop] = 0;
		}

		curJunc = rootVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			fragCnt++;
			if (curJunc->junc->type == frag_junction)
			{
				juncCnt++;
			}
			else if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
			{
				exonCnt++;
				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
				{
					exonSupportSum[iLoop] += curJunc->junc->support[iLoop];
				}
			}

			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
			{
				fragSupportSum[iLoop] += curJunc->junc->support[iLoop];
			}

			curJunc = curJunc->next;
		}
		rootVertex->junctionNum = juncCnt;
		rootVertex->exonNum = exonCnt;

		rootVertex->estimated = false;
		rootVertex->estimate_exonNum = 0; //leaf has nothing to estimate

// 		if (juncCnt <= 1)
// 		{
// 			//one junction or one exon
// 			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 			{
// 				rootVertex->support[iLoop] = rootVertex->junctionInRange->list->junc->support[iLoop];
// 			}
// 		}
// 		else
// 		{
// 			//more junctions
// 			//just calculate average for now, will change to transcript proportion later
// 			// 			for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
// 			// 			{
// 			// 				rootVertex->support[tmp] = fragSupportSum[tmp] / fragCnt;
// 			// 			}
// 		}
		for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
		{
			//count the exonic expression
			if (exonCnt > 0)
			{
				rootVertex->support[iLoop] = double(exonSupportSum[iLoop]) / exonCnt;
			} 
			else
			{
				rootVertex->support[iLoop] = 0.0;
			}
		}

		//rootVertex->anovaScore_support = Anova_test(rootVertex->support);

		delete [] fragSupportSum; delete [] exonSupportSum;
		return true;
	} 
	else
	{
		//extend children
		rootVertex->estimate_exonNum = 0;

		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			(rootVertex->childNum)++;
			curVertex = curEdge->linkedVertex;
			if (preCountGTree(curVertex) == false)
			{
			}
			if (rootVertex->childType != 3)
			{
				rootVertex->junctionNum += curVertex->junctionNum;
				rootVertex->exonNum += curVertex->exonNum;

				if (curVertex->estimated == false)
				{
					rootVertex->estimate_exonNum += curVertex->exonNum;
				} 
				else
				{
					rootVertex->estimate_exonNum += 1; //representative exon
				}
			} 

			curEdge = curEdge->next;
		}

		if (rootVertex->childType == 3)
		{
			juncCnt = 0;
			exonCnt = 0;
			fragCnt = 0;

			curJunc = rootVertex->junctionInRange->list;
			while (curJunc != NULL)
			{
				fragCnt++;
				if (curJunc->junc->type == frag_junction)
				{
					juncCnt++;
				}
				else if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
				{
					exonCnt++;
				}
				curJunc = curJunc->next;
			}
			rootVertex->junctionNum = juncCnt;
			rootVertex->exonNum = exonCnt;

			rootVertex->estimate_exonNum = exonCnt;
		}

		delete [] fragSupportSum; delete [] exonSupportSum;
		return false;
	}
}

long select_major_paths(GTvertex *rootVertex)
{
	long major_path_num = 0;
	double path_proportion_sum = 0.0, path_expression_sum = 0.0;
	GTedge *curEdge;
	GTvertex *curVertex;
	alternative_path *newAlterPath;
	int iLoop, expressedCnt;

	curEdge = rootVertex->child;
	while (curEdge != NULL)
	{
		curVertex = curEdge->linkedVertex;

		if (COUNT_MAJOR_PATHS_ONLY == true)
		{
			path_proportion_sum = 0.0;
			path_expression_sum = 0.0;
			expressedCnt = 0;
			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
			{
				path_proportion_sum += curVertex->proportion[iLoop];
				path_expression_sum += curVertex->support[iLoop];
				if (curVertex->support[iLoop] >= thresh_path_expressed)
				{
					++expressedCnt;
				}
			}

			if (path_proportion_sum / SUPPORT_VECTOR_SIZE >= 0.05 && path_expression_sum / SUPPORT_VECTOR_SIZE >= 3 && expressedCnt >= SUPPORT_VECTOR_SIZE/5)
			{
				//major path
				newAlterPath = new alternative_path;
				newAlterPath->path_start = curVertex->rangeLow;
				newAlterPath->path_end   = curVertex->rangeHigh;
				if (curVertex->prevSibling != NULL)
					newAlterPath->whole_path_start = curVertex->prevSibling->rangeLow;
				else
					newAlterPath->whole_path_start = curVertex->rangeLow;
				if (curVertex->nextSibling != NULL)
					newAlterPath->whole_path_end = curVertex->nextSibling->rangeHigh;
				else
					newAlterPath->whole_path_end = curVertex->rangeHigh;
				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
				{
					newAlterPath->support[iLoop] = curVertex->support[iLoop];
					newAlterPath->proportion[iLoop] = curVertex->proportion[iLoop];
				}
				newAlterPath->junctionNum = curVertex->junctionNum;
				newAlterPath->exonNum = curVertex->exonNum;
				newAlterPath->pathVertex = curVertex;

				++major_path_num;
				newAlterPath->next = rootVertex->major_alter_paths;
				rootVertex->major_alter_paths = newAlterPath;
			}
		} 
		else
		{
			//all paths
			newAlterPath = new alternative_path;
			newAlterPath->path_start = curVertex->rangeLow;
			newAlterPath->path_end   = curVertex->rangeHigh;
			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
			{
				newAlterPath->support[iLoop] = curVertex->support[iLoop];
				newAlterPath->proportion[iLoop] = curVertex->proportion[iLoop];
			}
			newAlterPath->junctionNum = curVertex->junctionNum;
			newAlterPath->exonNum = curVertex->exonNum;
			newAlterPath->pathVertex = curVertex;

			++major_path_num;
			newAlterPath->next = rootVertex->major_alter_paths;
			rootVertex->major_alter_paths = newAlterPath;
		}

		curEdge = curEdge->next;
	}

	rootVertex->major_alter_paths_num = major_path_num;
	return major_path_num;
}


bool countGTree(GTvertex *rootVertex)
{
	//collect basic counts on GTree for further statistics
	//collect child count and fragment count for every vertex
	//return true for leaves
	GTvertex *curVertex;
	GTedge *curEdge;
	alternative_path *curAlterPath;
	rangeJunction *curJunc;
	double *fragSupportSum = new double [SUPPORT_VECTOR_SIZE], *exonSupportSum = new double [SUPPORT_VECTOR_SIZE], totalSupport;
	int juncCnt, exonCnt, fragCnt, tmp, iLoop, jLoop;
	bool allChildrenAreLeaves;
	double *Parray_support, *Parray_proportion, *Qarray_support, *Qarray_proportion, **proportion_matrix, MSE_estimation;


	
	if (rootVertex->child == NULL)
	{
		//leaf
		juncCnt = 0;
		exonCnt = 0;
		fragCnt = 0;
		for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
		{
			fragSupportSum[iLoop] = 0;
			exonSupportSum[iLoop] = 0;
		}

		curJunc = rootVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			fragCnt++;
			if (curJunc->junc->type == frag_junction)
			{
				juncCnt++;
			}
			else if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
			{
				exonCnt++;
				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
				{
					exonSupportSum[iLoop] += curJunc->junc->support[iLoop];
				}
			}

			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
			{
				fragSupportSum[iLoop] += curJunc->junc->support[iLoop];
			}

			curJunc = curJunc->next;
		}
		rootVertex->junctionNum = juncCnt;
		rootVertex->exonNum = exonCnt;

		rootVertex->estimated = false;
		rootVertex->estimate_exonNum = 0; //leaf has nothing to estimate

// 		if (juncCnt <= 1)
// 		{
// 			//one junction or one exon
// 			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 			{
// 				rootVertex->support[iLoop] = rootVertex->junctionInRange->list->junc->support[iLoop];
// 			}
// 		}
// 		else
// 		{
// 			//more junctions
// 			//just calculate average for now, will change to transcript proportion later
// 			// 			for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
// 			// 			{
// 			// 				rootVertex->support[tmp] = fragSupportSum[tmp] / fragCnt;
// 			// 			}
// 		}
		for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
		{
			//count the exonic expression
			if (exonCnt > 0)
			{
				rootVertex->support[iLoop] = double(exonSupportSum[iLoop]) / exonCnt;
			} 
			else
			{
				rootVertex->support[iLoop] = 0.0;
			}
		}

		//rootVertex->anovaScore_support = Anova_test(rootVertex->support);


		//alternative splice sites
		if (fragCnt == 1 && juncCnt == 1 && rootVertex->junctionInRange->list->junc->alter != NULL)
		{
			vertexForAlterSpliceSite(rootVertex->junctionInRange->list->junc, rootVertex);
		}

		Parray_support = new double [SUPPORT_VECTOR_SIZE / 2 + 1];
		Qarray_support = new double [SUPPORT_VECTOR_SIZE / 2 + 1];

		for (iLoop = 1; iLoop <= SUPPORT_VECTOR_SIZE / 2; iLoop++)
		{
			Parray_support[iLoop] = rootVertex->support[iLoop - 1];
			Qarray_support[iLoop] = rootVertex->support[SUPPORT_VECTOR_SIZE / 2 + iLoop - 1];
		}

		delete [] Parray_support;
		delete [] Qarray_support;	

		calculate_ASM_group_meanExpression(rootVertex);

		delete [] fragSupportSum; delete [] exonSupportSum;

		return true;
	} 
	else
	{
		//extend children
		allChildrenAreLeaves = true;
		rootVertex->estimate_exonNum = 0;
		rootVertex->childNum = 0;
		rootVertex->junctionNum = 0;
		rootVertex->exonNum = 0;

		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			(rootVertex->childNum)++;
			curVertex = curEdge->linkedVertex;
			if (countGTree(curVertex) == false)
			{
				allChildrenAreLeaves = false;
			}
			if (rootVertex->childType != 3)
			{
				rootVertex->junctionNum += curVertex->junctionNum;
				rootVertex->exonNum += curVertex->exonNum;

				if (curVertex->estimated == false)
				{
					rootVertex->estimate_exonNum += curVertex->exonNum;
				} 
				else
				{
					rootVertex->estimate_exonNum += 1; //representative exon
				}
			} 

			curEdge = curEdge->next;
		}

		if (rootVertex->childType == 3)
		{
			juncCnt = 0;
			exonCnt = 0;
			fragCnt = 0;

			curJunc = rootVertex->junctionInRange->list;
			while (curJunc != NULL)
			{
				fragCnt++;
				if (curJunc->junc->type == frag_junction)
				{
					//alternative splice sites
					if (curJunc->junc->alter != NULL)
					{
						vertexForAlterSpliceSite(curJunc->junc, rootVertex);
					}

					juncCnt++;
				}
				else if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
				{
					exonCnt++;
				}
				curJunc = curJunc->next;
			}
			rootVertex->junctionNum = juncCnt;
			rootVertex->exonNum = exonCnt;

			rootVertex->estimate_exonNum = exonCnt;
		}

#ifdef DO_ESTIMATION
		//		if (allChildrenAreLeaves == true && (rootVertex->childType == 2 || rootVertex->childType == 3))
		if (rootVertex->childType == 2 || rootVertex->childType == 3)
		{
			for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
			{
				abundance_estimation(rootVertex, iLoop, MSE_estimation);
				rootVertex->MSE_estimation[iLoop] = MSE_estimation;
			}
			
			select_major_paths(rootVertex);

// 			if (SUPPORT_VECTOR_SIZE > 1 && rootVertex->major_alter_paths_num > 1)
// 			{
// 				Parray_support = new double [rootVertex->major_alter_paths_num + 1];
// 				Parray_proportion = new double [rootVertex->major_alter_paths_num + 1];
// 				Qarray_support = new double [rootVertex->major_alter_paths_num + 1];
// 				Qarray_proportion = new double [rootVertex->major_alter_paths_num + 1];
// 
// 				proportion_matrix = new double * [rootVertex->major_alter_paths_num + 1];
// 				for (iLoop = 0; iLoop <= rootVertex->major_alter_paths_num; iLoop++)
// 				{
// 					proportion_matrix[iLoop] = new double[SUPPORT_VECTOR_SIZE + 1];
// 				}
// 
// 				for (iLoop = 0; iLoop < rootVertex->major_alter_paths_num + 1; iLoop++)
// 				{
// 					Parray_support[iLoop] = 0.0;
// 					Parray_proportion[iLoop] = 0.0;
// 					Qarray_support[iLoop] = 0.0;
// 					Qarray_proportion[iLoop] = 0.0;
// 
// 					for (jLoop = 0; jLoop <= SUPPORT_VECTOR_SIZE; jLoop++)
// 					{
// 						proportion_matrix[iLoop][jLoop] = 0.0;
// 					}
// 				}
// 
// 				tmp = 1;
// 				curAlterPath = rootVertex->major_alter_paths;
// 				while (curAlterPath != NULL)
// 				{
// 					for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE/2; iLoop++)
// 					{
// 						Parray_support[tmp] += curAlterPath->support[iLoop];
// 						Parray_proportion[tmp] += curAlterPath->proportion[iLoop];
// 						Qarray_support[tmp] += curAlterPath->support[SUPPORT_VECTOR_SIZE / 2 + iLoop];
// 						Qarray_proportion[tmp] += curAlterPath->proportion[SUPPORT_VECTOR_SIZE / 2 + iLoop];
// 					}
// 
// 					for (iLoop = 1; iLoop <= SUPPORT_VECTOR_SIZE; iLoop++)
// 					{
// 						proportion_matrix[tmp][iLoop] = curAlterPath->proportion[iLoop - 1];
// 					}
// 
// 					tmp++;
// 					curAlterPath = curAlterPath->next;
// 				}
// 
// 				totalSupport = 0.0;
// 				for (iLoop = 1; iLoop <= rootVertex->major_alter_paths_num; iLoop++)
// 				{
// 					Parray_support[iLoop] = Parray_support[iLoop] / (SUPPORT_VECTOR_SIZE / 2);
// 					Parray_proportion[iLoop] = Parray_proportion[iLoop] / (SUPPORT_VECTOR_SIZE / 2);
// 					Qarray_support[iLoop] = Qarray_support[iLoop] / (SUPPORT_VECTOR_SIZE / 2);
// 					Qarray_proportion[iLoop] = Qarray_proportion[iLoop] / (SUPPORT_VECTOR_SIZE / 2);
// 
// 					totalSupport += Parray_support[iLoop] + Qarray_support[iLoop];
// 				}
// 
// // 				outputToPool(proportion_matrix, SUPPORT_VECTOR_SIZE, rootVertex->major_alter_paths_num);
// // 
// // 				rootVertex->JSDsqrt_mean = JSD_test(Parray_proportion, Qarray_proportion, rootVertex->major_alter_paths_num);
// // 				rootVertex->JSDsqrt_combination = JSDsqrt_matrix_combination(proportion_matrix, SUPPORT_VECTOR_SIZE, rootVertex->major_alter_paths_num);
// // 
// // 				rootVertex->JSD_Pvalue_mean = JSD_significance_level(rootVertex->JSDsqrt_mean, rootVertex->major_alter_paths_num, totalSupport);
// // 				rootVertex->JSD_Pvalue_combination = JSD_significance_level(pow(rootVertex->JSDsqrt_combination, 2), rootVertex->major_alter_paths_num, totalSupport);
// // 
// // 				rootVertex->divergence_stat = D_statistics(proportion_matrix, SUPPORT_VECTOR_SIZE, rootVertex->major_alter_paths_num);
// 
// 				if (SUPPORT_VECTOR_SIZE > 2)
// 				{
// // 					rootVertex->JSD_withinGroupRank_mean = withinMatrixRankTest_fullspace_JSD(proportion_matrix, SUPPORT_VECTOR_SIZE, rootVertex->major_alter_paths_num, rootVertex->JSDsqrt_mean, totalSupport);
// // 					rootVertex->divergence_stat_rank = withinMatrixRankTest_fullspace_Dstatistics(proportion_matrix, SUPPORT_VECTOR_SIZE, rootVertex->major_alter_paths_num, rootVertex->divergence_stat);
// 				} 
// 				else
// 				{
// 					rootVertex->JSD_withinGroupRank_mean = 1;
// 				}
// 
// 
// /*				rootVertex->JSDsqrt_mean = sqrt(rootVertex->JSDsqrt_mean);*/
// 
// 				//				rootVertex->spearmanCorr_support = Spearman_correlation(Parray_support, Qarray_support, rootVertex->childNum);
// 				//				rootVertex->spearmanCorr_proportion = Spearman_correlation(Parray_proportion, Qarray_proportion, rootVertex->childNum);
// // 				rootVertex->pearsonCorr_support = Pearson_correlation(Parray_support, Qarray_support, rootVertex->major_alter_paths_num);
// // 				rootVertex->pearsonCorr_proportion = Pearson_correlation(Parray_proportion, Qarray_proportion, rootVertex->major_alter_paths_num);
// 
// 				delete [] Parray_support;
// 				delete [] Parray_proportion;
// 				delete [] Qarray_support;
// 				delete [] Qarray_proportion;
// 
// 				for (iLoop = 0; iLoop <= rootVertex->major_alter_paths_num; iLoop++)
// 				{
// 					delete [] proportion_matrix[iLoop];
// 				}
// 				delete [] proportion_matrix;
// 			}
		}
		else if (rootVertex->childType == 1)
		{
// 			Parray_support = new double [rootVertex->childNum + 1];
// 			Qarray_support = new double [rootVertex->childNum + 1];
// 
// 			for (iLoop = 0; iLoop < rootVertex->childNum + 1; iLoop++)
// 			{
// 				Parray_support[iLoop] = 0.0;
// 				Qarray_support[iLoop] = 0.0;
// 			}
// 
// 			tmp = 1;
// 			curEdge = rootVertex->child;
// 			while (curEdge != NULL)
// 			{
// 				curVertex = curEdge->linkedVertex;
// 
// 				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE/2; iLoop++)
// 				{
// 					Parray_support[tmp] += curVertex->support[iLoop];
// 					Qarray_support[tmp] += curVertex->support[SUPPORT_VECTOR_SIZE / 2 + iLoop];
// 				}
// 
// 				tmp++;
// 				curEdge = curEdge->next;
// 			}
// 
// 			for (iLoop = 1; iLoop < rootVertex->childNum + 1; iLoop++)
// 			{
// 				Parray_support[iLoop] = Parray_support[iLoop] / (SUPPORT_VECTOR_SIZE / 2);
// 				Qarray_support[iLoop] = Qarray_support[iLoop] / (SUPPORT_VECTOR_SIZE / 2);
// 			}
// 
// //			rootVertex->spearmanCorr_support = Spearman_correlation(Parray_support, Qarray_support, rootVertex->childNum);
// 			rootVertex->pearsonCorr_support = Pearson_correlation(Parray_support, Qarray_support, rootVertex->childNum);
// 
// 			delete [] Parray_support;
// 			delete [] Qarray_support;
		}

		rootVertex->estimated = true;
		rootVertex->representative = makeRepresentative(rootVertex);


// 		Parray_support = new double [SUPPORT_VECTOR_SIZE / 2 + 1];
// 		Qarray_support = new double [SUPPORT_VECTOR_SIZE / 2 + 1];
// 
// 		for (iLoop = 1; iLoop <= SUPPORT_VECTOR_SIZE / 2; iLoop++)
// 		{
// 			Parray_support[iLoop] = rootVertex->support[iLoop - 1];
// 			Qarray_support[iLoop] = rootVertex->support[SUPPORT_VECTOR_SIZE / 2 + iLoop - 1];
// 		}
// 
// 		delete [] Parray_support;
// 		delete [] Qarray_support;	
		
#endif
		
		calculate_ASM_group_meanExpression(rootVertex);
//		rootVertex->JSD_reliable = alterSpliceReliability(rootVertex);

		delete [] fragSupportSum; delete [] exonSupportSum;

		return false;
	}
}

void output_vertex(GTvertex *rootVertex, ofstream *outputfile)
{
	//output given vertex
	GTvertex *curVertex;
	GTedge *curEdge;
	rangeJunction *curJunc;

	if (rootVertex->child == NULL)
	{
		//print fragment list
		curJunc = rootVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_exon)
			{
				(*outputfile) << "exon\t";
			} 
			else
			{
				(*outputfile) << "junc\t";
			}
			(*outputfile) << curJunc->junc->start << "\t" << curJunc->junc->end;
			curJunc = curJunc->next;
		}
		(*outputfile) << endl;
	} 
	else
	{
		//extend children
		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;
			output_vertex(curVertex, outputfile);

			curEdge = curEdge->next;
		}
	}

	return;
}


// bool alterSpliceReliability(GTvertex *curVertex)
// {
// 	int reliableCnt = 0;
// 
// 	for (int tmpCnt = 0; tmpCnt < SUPPORT_VECTOR_SIZE; tmpCnt++)
// 	{
// 		if (curVertex->support[tmpCnt] >= JSD_RELIABILITY_MIN_COVERAGE)
// 		{
// 			++reliableCnt;
// 		}
// 	}
// 	if (reliableCnt >= thresh_path_expressed)
// 	{
// 		return true;
// 	} 
// 	else
// 	{
// 		return false;
// 	}
// }


int alterSpliceCate(GTvertex *curVertex)
{
	//categorize alternative splicing events

	//if not, and the ASM has been assigned a category, then use that
	if (curVertex->ASMcategory > 0)
	{
		return curVertex->ASMcategory;
	}

	if ((curVertex->childType == 2 || curVertex->childType == 3) && curVertex->childNum > 1)
	{
		//ASM
	}
	else
	{
		return -1;
	}

	//SECOND, categorize the alternative splicing
	alternative_path *pathA, *pathB;

	if (curVertex->childType == 2 && curVertex->major_alter_paths_num > 1)
	{
		pathA = curVertex->major_alter_paths;
		while (pathA != NULL)
		{
			pathB = pathA->next;
			while (pathB != NULL)
			{
				if (abs(pathA->path_start - pathB->path_start) < 2 && abs(pathA->path_end - pathB->path_end) < 2)
				{
					if (pathA->exonNum == 0 && pathA->junctionNum == 1 && pathB->junctionNum >= 2 && pathB->exonNum >= 1 || pathA->junctionNum >= 2 && pathA->exonNum >= 1 && pathB->exonNum == 0 && pathB->junctionNum == 1)
					{
						curVertex->ASMcategory = exon_skipping;						
						return exon_skipping;
					}
					if (pathA->exonNum >= 1 && pathA->junctionNum >= 2 && pathB->exonNum >= 1 && pathB->junctionNum >= 2)
					{
						curVertex->ASMcategory = mutual_exclusive;						
					}
					if (curVertex->major_alter_paths_num == 2 && (pathA->exonNum == 1 && pathA->junctionNum == 0 && pathB->exonNum == 0 && pathB->junctionNum == 1 || pathA->exonNum == 0 && pathA->junctionNum == 1 && pathB->exonNum == 1 && pathB->junctionNum == 0))
					{
						curVertex->ASMcategory = intron_retention;
						return intron_retention;
					}
				}
				else if (abs(pathA->path_start - pathB->path_start) < 2)
				{
					curVertex->ASMcategory = alter_end;	
				}
				else if (abs(pathA->path_end - pathB->path_end) < 2)
				{
					curVertex->ASMcategory = alter_start;	
				}

				pathB = pathB->next;
			}

			pathA = pathA->next;
		}
	}

	if (curVertex->ASMcategory == mutual_exclusive || curVertex->ASMcategory == alter_end || curVertex->ASMcategory == alter_start)
	{
		//let mutual exclusive have lower priority than exon_skipping, because in the loop above we check children pair by pair
		return curVertex->ASMcategory;
	}

// 	if (curVertex->childType == 2 && curVertex->childNum > 1)
// 	{
// 		edgeA = curVertex->child;
// 		while (edgeA != NULL)
// 		{
// 			vertexA = edgeA->linkedVertex;
// 
// 			edgeB = edgeA->next;
// 			while (edgeB != NULL)
// 			{
// 				vertexB = edgeB->linkedVertex;
// 
// 				if (vertexA->rangeLow == vertexB->rangeLow && vertexA->rangeHigh == vertexB->rangeHigh)
// 				{
// 					if (vertexA->childNum == 0 && vertexA->junctionNum == 1 && vertexB->childNum == 3 && vertexB->junctionNum == 2 || vertexA->childNum == 3 && vertexA->junctionNum == 2 && vertexB->childNum == 0 && vertexB->junctionNum == 1)
// 					{
// 						output_vertex(curVertex, &exon_skipping_file);
// 						exon_skipping_file << endl;
// 						return exon_skipping;
// 					}
// 					if (vertexA->childNum == 3 && vertexA->junctionNum == 2 && vertexB->childNum == 3 && vertexB->junctionNum == 2)
// 					{
// 						//if (vertexA->junctionInRange->list->junc->end > vertexB->junctionInRange->list->next->junc->start || vertexB->junctionInRange->list->junc->end > vertexA->junctionInRange->list->next->junc->start)
// 						{
// 							output_vertex(curVertex, &mutual_exclusive_file);
// 							mutual_exclusive_file << endl;
// 							return mutual_exclusive;
// 						}
// 					}
// 				}
// 
// 				edgeB = edgeB->next;
// 			}
// 
// 			edgeA = edgeA->next;
// 		}
// 	}

	curVertex->ASMcategory = unknown;
	return unknown;
}

double logFactorial(long n)
{
	//calculate the log of the factorial of n
	double logfact = 0.0;

	for (long iLoop = 2; iLoop <= n; iLoop++)
	{
		logfact += log(double(iLoop));
	}

	return logfact;
}

long count_estimate_juncNum(GTvertex *targetVertex)
{
	long juncCnt = 0;
	rangeJunction *curJunc;
	GTedge *curEdge;
	GTvertex *curVertex;

	if (targetVertex->childType == 3)
	{
		curJunc = targetVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_junction)
			{
				++juncCnt;
			}
			curJunc = curJunc->next;
		}
	}
	else if (targetVertex->childType == 2)
	{
		curEdge = targetVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;
			if (curVertex->childType == 0 && curVertex->exonNum == 0 && curVertex->junctionNum == 1)
			{
				++juncCnt;
			} 

			curEdge = curEdge->next;
		}
	}

	return juncCnt;
}


double abundance_estimation(GTvertex *rootVertex, int index_tissue, double &MSE)
{
	//return the log likelihood
	//8/16/2011 MSE measure the mean squared error of the estimated coverage versus the observed coverage in the alternatively spliced region (i.e., excluding exons adjacent to the ASM)
	//9/12/2011 New estimation method

	if (rootVertex->child == NULL || rootVertex->childType == 1 || rootVertex->childType == 0)
	{
		return 0.0;
	}

	//transcript abundance estimation
	long exonCnt = 0, junctionCnt = 0, fragmentCnt_est = 0, transcriptCnt = 0;
	GTvertex *curVertex, *siblingVertex;
	GTedge *curEdge;
	rangeJunction *curJunc;
	long iLoop, jLoop, kLoop;
	bool hasPrevExon = false, hasNextExon = false;
	const int MaxPrevFragNum = 1, MaxNextFragNum = 1;
	int prevFragCnt = 0, nextFragCnt = 0, readlength = 50, fakeJunctionLength = readlength - 5;

	curEdge = rootVertex->child;
	while (curEdge != NULL)
	{
		if (curEdge->linkedVertex->childType == 0 && curEdge->linkedVertex->junctionNum > MAXJUNCNUMINDEPPATH)
		{
			return 0.0;
		}
		curEdge = curEdge->next;
	}

	if (rootVertex->estimate_exonNum > 0)
	{
		exonCnt = rootVertex->estimate_exonNum;
	} 
	else
	{
		exonCnt = rootVertex->exonNum;
	}
	junctionCnt = count_estimate_juncNum(rootVertex);
	fragmentCnt_est = exonCnt + junctionCnt;
	transcriptCnt = rootVertex->childNum;

	bool **A = new bool * [transcriptCnt]; //matrix A, transcript-exon matrix
	fragment **fragArray = new fragment * [fragmentCnt_est + MaxPrevFragNum + MaxNextFragNum]; //exon array. exonArray[0] is the exon prev to rootVertex, and exonArray[1+exonCnt] is the exon after the rootVertex.
	long *transFragCnt = new long [transcriptCnt]; //frag count for every transcript
	double *transLength = new double [transcriptCnt]; //transcript length
	for (iLoop = 0; iLoop < transcriptCnt; ++iLoop)
	{
		A[iLoop] = new bool [fragmentCnt_est + MaxPrevFragNum + MaxNextFragNum];
		transFragCnt[iLoop] = 0;
		transLength[iLoop] = 0.0;
	}

	double *errorRatio_transExpr = new double [transcriptCnt];
	long *errorRatio_fragCnt = new long [transcriptCnt];
	for (iLoop = 0; iLoop < transcriptCnt; ++iLoop)
	{
		errorRatio_transExpr[iLoop] = 0.0;
		errorRatio_fragCnt[iLoop] = 0;
	}


	//build exonArray
	siblingVertex = rootVertex;
	while (prevFragCnt < MaxPrevFragNum && siblingVertex->prevSibling != NULL && siblingVertex->prevSibling->childType == 0)
	{
		//there is an exon right before the siblingVertex
		if (siblingVertex->prevSibling->childNum != 0 || siblingVertex->prevSibling->junctionInRange == NULL)// || rootVertex->prevSibling->junctionNum > 0)
		{
			cout << "Abnormal case in abundance estimation: Unrecognized prev exon." << endl;
			exit(1);
		} 
		else if (siblingVertex->prevSibling->junctionInRange->list->junc->type == frag_exon || siblingVertex->prevSibling->junctionInRange->list->junc->type == frag_retained_intron
			)//|| siblingVertex->prevSibling->junctionInRange->list->junc->type == frag_junction)
		{
			fragArray[prevFragCnt++] = siblingVertex->prevSibling->junctionInRange->list->junc;
			hasPrevExon = true;
		}

		siblingVertex = siblingVertex->prevSibling;
	}

	if (rootVertex->childType == 3)
	{
		if (hasPrevExon == true)
			iLoop = prevFragCnt;
		else
			iLoop = 0;
		//build fragArray from range junction list
		curJunc = rootVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
			{
				fragArray[iLoop++] = curJunc->junc;
			}
			curJunc = curJunc->next;
		}
	}
	else if (rootVertex->childType == 2)
	{
		if (hasPrevExon == true)
			iLoop = prevFragCnt;
		else
			iLoop = 0;
		//build fragArray from all children
		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;
			if (curVertex->estimated == false)
			{
				curJunc = curVertex->junctionInRange->list;
				while (curJunc != NULL)
				{
					if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
					{
						fragArray[iLoop++] = curJunc->junc;
					}
					curJunc = curJunc->next;
				}
			} 
			else
			{
				fragArray[iLoop++] = curVertex->representative;
			}

			curEdge = curEdge->next;
		}
	}


	siblingVertex = rootVertex;
	while (nextFragCnt < MaxNextFragNum && siblingVertex->nextSibling != NULL && siblingVertex->nextSibling->childType == 0)
	{
		//there is an exon right after the siblingVertex
		if (siblingVertex->nextSibling->childNum != 0 || siblingVertex->nextSibling->junctionInRange == NULL)// || rootVertex->prevSibling->junctionNum > 0)
		{
			cout << "Abnormal case in abundance estimation: Unrecognized next exon." << endl;
			exit(1);
		} 
		else if (siblingVertex->nextSibling->junctionInRange->list->junc->type == frag_exon || siblingVertex->nextSibling->junctionInRange->list->junc->type == frag_retained_intron
			)//|| siblingVertex->nextSibling->junctionInRange->list->junc->type == frag_junction)
		{
			fragArray[prevFragCnt + fragmentCnt_est + nextFragCnt++] = siblingVertex->nextSibling->junctionInRange->list->junc;
			hasNextExon = true;
		}

		siblingVertex = siblingVertex->nextSibling;
	}


	//build transcript-exon matrix
	curEdge = rootVertex->child;
	iLoop = 0;
	while (curEdge != NULL)
	{
		curVertex = curEdge->linkedVertex;

		for (jLoop = 0; jLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; jLoop++)
		{
			A[iLoop][jLoop] = false;
		}
		if (hasPrevExon == true)
		{
			if (abs(curVertex->rangeLow - rootVertex->rangeLow) <= 1)
			{
				for (jLoop = 0; jLoop < prevFragCnt; ++jLoop)
				{
					A[iLoop][jLoop] = true;
					transFragCnt[iLoop] += 1;
					if (fragArray[jLoop]->type == frag_exon || fragArray[jLoop]->type == frag_retained_intron)
						transLength[iLoop] += abs(fragArray[jLoop]->end - fragArray[jLoop]->start + 1);
					if (fragArray[jLoop]->type == frag_junction)
						transLength[iLoop] += fakeJunctionLength;
				}
			}
		} 
		if (hasNextExon == true)
		{
			if (abs(rootVertex->rangeHigh - curVertex->rangeHigh) <= 1)
			{
				for (jLoop = 0; jLoop < nextFragCnt; ++jLoop)
				{
					A[iLoop][prevFragCnt + fragmentCnt_est + jLoop] = true;
					transFragCnt[iLoop] += 1;
					if (fragArray[prevFragCnt + fragmentCnt_est + jLoop]->type == frag_exon || fragArray[prevFragCnt + fragmentCnt_est + jLoop]->type == frag_retained_intron)
						transLength[iLoop] += abs(fragArray[prevFragCnt + fragmentCnt_est + jLoop]->end - fragArray[prevFragCnt + fragmentCnt_est + jLoop]->start + 1);
					if (fragArray[prevFragCnt + fragmentCnt_est + jLoop]->type == frag_junction)
						transLength[iLoop] += fakeJunctionLength;
				}
			}
		} 

		if (curVertex->estimated == false)
		{
			curJunc = curVertex->junctionInRange->list;
			while (curJunc != NULL)
			{
				if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
				{
					for (jLoop = prevFragCnt; jLoop < prevFragCnt + fragmentCnt_est; jLoop++)
					{
						if (fragArray[jLoop]->ID == curJunc->junc->ID)
						{
							A[iLoop][jLoop] = true;
							transFragCnt[iLoop] += 1;
							if (fragArray[jLoop]->type == frag_exon || fragArray[jLoop]->type == frag_retained_intron)
								transLength[iLoop] += abs(fragArray[jLoop]->end - fragArray[jLoop]->start + 1);
							if (fragArray[jLoop]->type == frag_junction)
								transLength[iLoop] += fakeJunctionLength;

							break;
						} 
					}

					errorRatio_transExpr[iLoop] += curJunc->junc->support[index_tissue];
					errorRatio_fragCnt[iLoop] += 1;
				}
				curJunc = curJunc->next;
			}
		} 
		else
		{
			for (jLoop = prevFragCnt; jLoop < prevFragCnt + fragmentCnt_est; jLoop++)
			{
				if (fragArray[jLoop]->ID == curVertex->representative->ID)
				{
					A[iLoop][jLoop] = true;
					transFragCnt[iLoop] += 1;
					transLength[iLoop] += abs(fragArray[jLoop]->end - fragArray[jLoop]->start + 1);

					break;
				} 
			}

			errorRatio_transExpr[iLoop] += curVertex->representative->support[index_tissue];
			errorRatio_fragCnt[iLoop] += 1;
		}

		iLoop++;
		curEdge = curEdge->next;
	}



	//estimation
	double *P = new double [fragmentCnt_est + MaxPrevFragNum + MaxNextFragNum]; // probability that reads fall into each exon
	double *Y = new double [fragmentCnt_est + MaxPrevFragNum + MaxNextFragNum]; // observed exon expression
	double *Q = new double [transcriptCnt]; //proportions of alternative transcripts
	double *Qold = new double [transcriptCnt]; //new Q vector
	double *Qtmp = new double [transcriptCnt];
	double *Qdelta = new double [transcriptCnt]; //change on Q
	double *Z = new double [transcriptCnt]; //expectation of read count, will be derived from C
	double *C = new double [transcriptCnt]; //expectation of read coverage
	double lambda; 


	//fill in Y
	for (iLoop = 0; iLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; iLoop++)
	{
		if ((fragArray[iLoop]->support)[index_tissue] < 0.0)
			(fragArray[iLoop]->support)[index_tissue] = 0.0;
		Y[iLoop] = (fragArray[iLoop]->support)[index_tissue];
	}

	//initiate Q, Z, transExonCnt
	for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
	{
		Q[iLoop] = 1.0 / transcriptCnt;
		Qtmp[iLoop] = 0.0;
		Qold[iLoop] = 0.0;
		Qdelta[iLoop] = 0.0;

		Z[iLoop] = 0.0;
		C[iLoop] = 0.0;
	}


	//calculate P
	long lengthSum = 0;
	for (iLoop = 0; iLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; iLoop++)
	{
		if (fragArray[iLoop]->type == frag_exon || fragArray[iLoop]->type == frag_retained_intron)
		{
			lengthSum += abs(fragArray[iLoop]->end - fragArray[iLoop]->start + 1);
		} 
	}
	for (iLoop = 0; iLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; iLoop++)
	{
		if (fragArray[iLoop]->type == frag_exon || fragArray[iLoop]->type == frag_retained_intron)
		{
			P[iLoop] = double(abs(fragArray[iLoop]->end - fragArray[iLoop]->start + 1)) / lengthSum;
		}
		else if (fragArray[iLoop]->type == frag_junction)
		{
			P[iLoop] = double(fakeJunctionLength) / lengthSum;
		}
	}



	//EM
	double sum_count, sum_denominator, sum_numerator, tmpvalue;
	bool stoppable = false, tmpFlag;
	long MAX_ITER = 100, iterCnt = 0;

	//12/18/2012 limit the max # of transcripts in an estimation
	if (transcriptCnt > 30)
		iterCnt = MAX_ITER;

	while (iterCnt < MAX_ITER && stoppable == false)
	{
		for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
		{
			Qold[iLoop] = Q[iLoop];
		}

		//E-step
		for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
		{
			sum_count = 0.0;

			for (jLoop = 0; jLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; jLoop++)
			{
				if (A[iLoop][jLoop] == true)
				{
					sum_denominator = 0.0;
					for (kLoop = 0; kLoop < transcriptCnt; kLoop++)
					{
						if (A[kLoop][jLoop] == true)
						{
							sum_denominator += P[jLoop] * Q[kLoop];
						}
					}

					if (sum_denominator > 0.0)
					{
						tmpvalue = P[jLoop] * Q[iLoop] * Y[jLoop] / sum_denominator;
						sum_count += tmpvalue;
					}
				}
			}

			C[iLoop] = sum_count / transFragCnt[iLoop];
			Z[iLoop] = C[iLoop] * transLength[iLoop] / readlength;
		}


		//M-step
		//calculate lambda
		lambda = 0.0;
		for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
		{
			lambda += Z[iLoop];
		}

		//calculate Q
		//get Q iteratively until stable
		if (transcriptCnt > 1)
		{
			tmpFlag = false;
			while (tmpFlag == false)
			{
				for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
				{
					Qtmp[iLoop] = Q[iLoop];
				}
				for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
				{
					sum_denominator = 0.0;
					for (jLoop = 0; jLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; jLoop++)
					{
						if (A[iLoop][jLoop] == true)
						{
							sum_denominator += P[jLoop];
						}
					}
					sum_denominator *= (lambda - Z[iLoop]);

					sum_numerator = 0.0;
					for (jLoop = 0; jLoop < transcriptCnt; jLoop++)
					{
						if (jLoop == iLoop)
						{
							continue;
						}

						sum_count = 0.0;
						for (kLoop = 0; kLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; kLoop++)
						{
							if (A[jLoop][kLoop] == true)
							{
								sum_count += P[kLoop];
							}
						}
						sum_count *= Qtmp[jLoop];

						sum_numerator += sum_count;
					}
					sum_numerator *= Z[iLoop];

					if (sum_denominator > 0.0)
					{
						Qtmp[iLoop] = sum_numerator / sum_denominator;
					} 
					else
					{
						if (lambda == Z[iLoop])
						{
							Qtmp[iLoop] = 1.0;
						}
						else
						{
							Qtmp[iLoop] = 0.0;
						}
					}
				}

				sum_count = 0.0;
				for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
				{
					sum_count += Qtmp[iLoop];
				}
				for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
				{
					//calculate change on Q
					if (sum_count > 0.0)
					{
						Qtmp[iLoop] = Qtmp[iLoop] / sum_count;
					} 
					else
					{
						Qtmp[iLoop] = 1.0 / transcriptCnt;
					}
					Qdelta[iLoop] = fabs(Qtmp[iLoop] - Q[iLoop]);
					Q[iLoop] = Qtmp[iLoop];
				}

				tmpFlag = true;
				for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
				{
					if (Q[iLoop] > 0.0 && fabs(Qdelta[iLoop] / Q[iLoop]) > 0.001)
					{
						tmpFlag = false;
						break;
					}
				}
			}

			for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
			{
				//calculate change on Q
				Qdelta[iLoop] = fabs(Qold[iLoop] - Q[iLoop]);
			}

			iterCnt++;
			tmpFlag = true;
			for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
			{
				if (Q[iLoop] > 0.0 && fabs(Qdelta[iLoop] / Q[iLoop]) > 0.001)
				{
					tmpFlag = false;
					break;
				}
			}
			if (iterCnt >= 3 && tmpFlag == true)
			{
				stoppable = true;
			}	
		} 
		else
		{
			stoppable = true;
		}

	}

	//one more E-step
	for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
	{
		sum_count = 0.0;

		for (jLoop = 0; jLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; jLoop++)
		{
			if (A[iLoop][jLoop] == true)
			{
				sum_denominator = 0.0;
				for (kLoop = 0; kLoop < transcriptCnt; kLoop++)
				{
					if (A[kLoop][jLoop] == true)
					{
						sum_denominator += P[jLoop] * Q[kLoop];
					}
				}

				if (sum_denominator > 0.0)
				{
					sum_count += P[jLoop] * Q[iLoop] * Y[jLoop] / sum_denominator;
				}
			}
		}

		C[iLoop] = sum_count / transFragCnt[iLoop];
		Z[iLoop] = C[iLoop] * transLength[iLoop] / readlength;
	}


	for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
	{
		if (Z[iLoop] < 0)
		{
			cout << "neg_Z\t";
		}
	}


	//copy results & calculate MSE
	MSE = 0.0;
	double min_path_support = MAX_CHR_LENGTH;
	curEdge = rootVertex->child;
	iLoop = 0;
	while (curEdge != NULL)
	{
		curVertex = curEdge->linkedVertex;

		//!!!Do not use support for the whole transcript fragment here! Use only the portion belonging to the vertex!
		//		curVertex->support[index_tissue] = Z[iLoop];
		if (transFragCnt[iLoop] > 0)
		{
			curVertex->support[index_tissue] = C[iLoop]; //Z[iLoop] / transFragCnt[iLoop]; //evenly divided, for coverage only

			if (C[iLoop] < min_path_support)
				min_path_support = C[iLoop];
		} 
		else
		{
			curVertex->support[index_tissue] = 0.0;
		}

		if (errorRatio_fragCnt[iLoop] > 0 && curVertex->support[index_tissue] > 10 && abs(errorRatio_transExpr[iLoop]/errorRatio_fragCnt[iLoop] - curVertex->support[index_tissue]) / curVertex->support[index_tissue] > MSE)
		{
			curVertex->obs_support[index_tissue] = errorRatio_transExpr[iLoop]/errorRatio_fragCnt[iLoop];
			MSE = abs(errorRatio_transExpr[iLoop]/errorRatio_fragCnt[iLoop] - curVertex->support[index_tissue]) / curVertex->support[index_tissue];
		}

		curVertex->proportion[index_tissue] = Q[iLoop];

		if (Q[iLoop] < 0)
		{
			cout << "neg_Q\t";
		}

		iLoop++;
		curEdge = curEdge->next;
	}

	rootVertex->min_path_support[index_tissue] = min_path_support;

	//calculate the log likelihood
	double loglikelihood = 0.0, lambda_i;

	// 	lambda = 0.0;
	// 	for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
	// 	{
	// 		lambda += Z[iLoop];
	// 	}
	// 
	// 	sum_numerator = 0.0;
	// 	for (jLoop = 0; jLoop < transcriptCnt; jLoop++)
	// 	{
	// 		sum_count = 0.0;
	// 		for (kLoop = 0; kLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; kLoop++)
	// 		{
	// 			if (A[jLoop][kLoop] == true)
	// 			{
	// 				sum_count += P[kLoop];
	// 			}
	// 		}
	// 		sum_count *= Q[jLoop];
	// 
	// 		sum_numerator += sum_count;
	// 	}
	// 
	// 	for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
	// 	{
	// 		sum_count = 0.0;
	// 		for (kLoop = 0; kLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; kLoop++)
	// 		{
	// 			if (A[iLoop][kLoop] == true)
	// 			{
	// 				sum_count += P[kLoop];
	// 			}
	// 		}
	// 		sum_count *= Q[iLoop];
	// 		if (sum_numerator > 0.0)
	// 		{
	// 			lambda_i = sum_count / sum_numerator;
	// 		} 
	// 		else
	// 		{
	// 			lambda_i = 0.0;
	// 		}
	// 
	// 		sum_count = 0.0;
	// 		for (kLoop = 0; kLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; kLoop++)
	// 		{
	// 			if (A[iLoop][kLoop] == true)
	// 			{
	// 				sum_count += P[kLoop];
	// 			}
	// 		}
	// 
	// 		loglikelihood += -lambda_i + Z[iLoop] * log(lambda_i);
	// 
	// 		for (kLoop = 0; kLoop < prevFragCnt + fragmentCnt_est + nextFragCnt; kLoop++)
	// 		{
	// 			if (A[iLoop][kLoop] == true)
	// 			{
	// 				if (long(Z[iLoop]/transFragCnt[iLoop] >= 0))
	// 				{
	// 					// 					if (ceil(Z[iLoop]/transExonCnt[iLoop]) > 1000)
	// 					// 					{
	// 					// 						cout << "Warning: calculating factorial of a number larger than 10" << endl;
	// 					// 					}
	// 					//loglikelihood += - log(double(factorial(long(ceil(Z[iLoop]/transExonCnt[iLoop])))));
	// 					loglikelihood += - logFactorial(long(ceil(Z[iLoop]/transFragCnt[iLoop])));
	// 				}
	// 				else
	// 				{
	// 					cout << "ERROR: estimation, factorial" << endl;
	// 				}
	// 				loglikelihood += Z[iLoop]/transFragCnt[iLoop] * log(P[kLoop] / sum_count);
	// 			}
	// 		}
	// 	}


	//clean up
	delete [] fragArray;
	delete [] transFragCnt;
	delete [] errorRatio_fragCnt;
	delete [] errorRatio_transExpr;
	for (iLoop = 0; iLoop < transcriptCnt; iLoop++)
	{
		delete [] A[iLoop];
	}
	delete [] A;
	delete [] P;
	delete [] Y;
	delete [] Z;
	delete [] Q;
	delete [] Qold;
	delete [] Qtmp;
	delete [] Qdelta;
	delete [] C;
	delete [] transLength;

	return loglikelihood;
}




/************************************************************************/
/* MAIN                                                                 */
/************************************************************************/



void input_junction(string junctionFilename)
{
	ifstream junctionFile;
	junctionFile.open(junctionFilename.c_str());

	fragment *newJunction;
	alter_junction *newAlterJunc;
	rangeJunction *newRangeJunc;
	string name;
	bool junctionExist, thresholdValid, strand;
	long temp, support1, support2;
	string info;

	while (junctionFile >> name)
	{
		newJunction = new fragment;
		newJunction->type = frag_junction;

		newJunction->chromosome_start = name;
		newJunction->chromosome_end = name;
		junctionFile >> newJunction->start;
		junctionFile >> newJunction->end;
		junctionFile >> strand; newJunction->transDirection = strand == true ? sense : antisense;
		getline(junctionFile, info);

#ifndef TRANSCRIPTION_DIRECTION
		newJunction->transDirection = sense;
#endif
		
		thresholdValid = true;

		//remove small indel
		if (newJunction->end - newJunction->start - 1 < min_junction_intron_length)
			thresholdValid = false;

		//if more ...


		if (thresholdValid == true)
		{
			junctionExist = false;
			for (temp = sortList_Junction_Num>0?sortList_Junction_Num:1; temp <= sortList_Junction_Num; temp++)
			{
				if (sortList_Junction[temp]->junc->chromosome_start.compare(newJunction->chromosome_start) == 0 && sortList_Junction[temp]->junc->chromosome_end.compare(newJunction->chromosome_end) == 0 && sortList_Junction[temp]->junc->start == newJunction->start && sortList_Junction[temp]->junc->end == newJunction->end)
				{
					if (sortList_Junction[temp]->junc->transDirection == newJunction->transDirection)
					{
						junctionExist = true;
						delete newJunction;
						break;
					}
					else
					{
						CLEAN_SPLICE_DIRECTION = true;
						//cout << "Warning: two junctions at the same location but of different transcription direction: " << newJunction->start << " - " << newJunction->end << "\n";
					}
				}

				//search for alternative splicing junctions
				// 				if (strcmp(sortList_Junction[temp]->junc->chromosome_start, newJunction->chromosome_start) == 0 && strcmp(sortList_Junction[temp]->junc->chromosome_end, newJunction->chromosome_end) == 0 && abs(sortList_Junction[temp]->junc->start - newJunction->start) <= 100 && sortList_Junction[temp]->junc->end == newJunction->end)
				// 				{
				// 					junctionExist = true;
				// 					//alter 5'
				// 					newAlterJunc = new alter_junction;
				// 					newAlterJunc->category = true;
				// 					newAlterJunc->juncInfo = newJunction;
				// 					newAlterJunc->next = sortList_Junction[temp]->junc->alter;
				// 					sortList_Junction[temp]->junc->alter = newAlterJunc;
				// 				}
				// 				else if (strcmp(sortList_Junction[temp]->junc->chromosome_start, newJunction->chromosome_start) == 0 && strcmp(sortList_Junction[temp]->junc->chromosome_end, newJunction->chromosome_end) == 0 && sortList_Junction[temp]->junc->start == newJunction->start && abs(sortList_Junction[temp]->junc->end - newJunction->end) <= 100)
				// 				{
				// 					junctionExist = true;
				// 					//alter 3'
				// 					newAlterJunc = new alter_junction;
				// 					newAlterJunc->category = false;
				// 					newAlterJunc->juncInfo = newJunction;
				// 					newAlterJunc->next = sortList_Junction[temp]->junc->alter;
				// 					sortList_Junction[temp]->junc->alter = newAlterJunc;
				// 				}
			}

			if (junctionExist == false)
			{
				newJunction->ID = ++fragmentID_Cnt;
				newRangeJunc = new rangeJunction;
				newRangeJunc->junc = newJunction;

				sortList_Junction_Num++;
				if (sortList_Junction_Num >= sortList_Junction.size())
				{
					sortList_Junction.resize(sortList_Junction.size() + DEFAULT_JUNCTION_NUM, NULL);
					sortKey_Junction.resize(sortKey_Junction.size() + DEFAULT_JUNCTION_NUM, 0);
					mergeSort_Larray.resize(mergeSort_Larray.size() + DEFAULT_JUNCTION_NUM, 0);
					mergeSort_Rarray.resize(mergeSort_Rarray.size() + DEFAULT_JUNCTION_NUM, 0);
					mergeSort_LorderedList.resize(mergeSort_LorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
					mergeSort_RorderedList.resize(mergeSort_RorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
				}
				sortList_Junction[sortList_Junction_Num] = newRangeJunc;
			}
		}	
	}

	junctionFile.close();

	return;
}


void getJunctionSupport(int index)
{
	//get junction support from fragent file
	ifstream junctionSupportFile;
	string inputfilename;
#ifdef UNIX
	inputfilename = inputPath + "junction_" + itostr(index+1) + ".txt";
#else
	inputfilename = "tmp\\junction_" + itostr(index+1) + ".txt";
#endif
	junctionSupportFile.open(inputfilename.c_str());

	long i;
	for (i = 1; i <= sortList_Junction_Num; i++)
	{
		sortKey_Junction[i] = sortList_Junction[i]->junc->end;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);

	for (i = 1; i <= sortList_Junction_Num; i++)
	{
		sortKey_Junction[i] = sortList_Junction[i]->junc->start;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);


	long startPosition, endPosition, juncIndex, tmpCount;
	string chromosome;
	string info;
	bool strand;
	startPosition = 0;
	juncIndex = 1;

	while (junctionSupportFile >> chromosome)
	{
		junctionSupportFile >> startPosition;
		junctionSupportFile >> endPosition;
		junctionSupportFile >> strand;

#ifndef TRANSCRIPTION_DIRECTION
		strand = true;
#endif

		getline(junctionSupportFile, info);

		while (juncIndex <= sortList_Junction_Num && startPosition > sortList_Junction[juncIndex]->junc->start)
		{
			juncIndex++;
		}

		while (juncIndex <= sortList_Junction_Num && startPosition == sortList_Junction[juncIndex]->junc->start && endPosition > sortList_Junction[juncIndex]->junc->end)
		{
			juncIndex++;
		}

		if (juncIndex > sortList_Junction_Num)
		{
			break;
		}

		for (tmpCount = 0; juncIndex + tmpCount <= sortList_Junction_Num && startPosition == sortList_Junction[juncIndex+tmpCount]->junc->start && endPosition == sortList_Junction[juncIndex+tmpCount]->junc->end; ++tmpCount)
		{
			if ((strand == true ? sense : antisense) == sortList_Junction[juncIndex+tmpCount]->junc->transDirection)
			{
				((sortList_Junction[juncIndex+tmpCount]->junc->support)[index]) += 1;
			}
		}
	}

	junctionSupportFile.close();

	return;
}

//sort splice sites
void merge_SpliceSiteSort(long p, long q, long r)
{
	long n1, n2, i, j, k;

	n1 = q - p + 1;
	n2 = r - q;

	for (i = 1; i <= n1; i++)
	{
		mergeSort_Larray[i] = sortKey_spliceSite[p + i - 1];
		mergeSort_LorderedList_spliceSite[i] = sortList_spliceSite[p + i - 1];
	}
	for (j = 1; j <= n2; j++)
	{
		mergeSort_Rarray[j] = sortKey_spliceSite[q + j];
		mergeSort_RorderedList_spliceSite[j] = sortList_spliceSite[q + j];
	}

	mergeSort_Larray[n1 + 1] = MAX_CHR_LENGTH * 2;
	mergeSort_Rarray[n2 + 1] = MAX_CHR_LENGTH * 2;

	i = 1;
	j = 1;

	for (k = p; k <= r; k++)
	{
		if (mergeSort_Larray[i] <= mergeSort_Rarray[j])
		{
			sortKey_spliceSite[k] = mergeSort_Larray[i];
			sortList_spliceSite[k] = mergeSort_LorderedList_spliceSite[i];

			i++;
		} 
		else
		{
			sortKey_spliceSite[k] = mergeSort_Rarray[j];
			sortList_spliceSite[k] = mergeSort_RorderedList_spliceSite[j];

			j++;
		}
	}

	return;
}


void mergeSort_SpliceSiteSort(long sortList_size)
{
	//non-recursive merge sort for sorting junctions
	long m, n, i, r;
	m = 1;
	n = sortList_size;

	while (m <= n)
	{
		i = 1;
		while (i <= n - m)
		{
			r = (i + 2 * m - 1) < n ? (i + 2 * m - 1) : n;
			merge_SpliceSiteSort(i, i + m - 1, r);
			i = i + 2 * m;
		}

		m = m * 2;
	}

	return;
}

void getExons()
{
	spliceSite *curSite;
	fragment *newFragment;
	rangeJunction *newRangeJunc;
	long i;

	if (sortList_Junction_Num == 0)
	{
		newFragment = new fragment;
		newFragment->type = frag_exon;
		newFragment->altersite = 0;
		newFragment->frag_name = "Exon_1";
		newFragment->ID = ++fragmentID_Cnt;

		//curSite = sortList_spliceSite[1];
		newFragment->chromosome_start = resultNamePrefix;
		newFragment->chromosome_end = resultNamePrefix;
		newFragment->start = CHROMOSOME_START;
		newFragment->end = CHROMOSOME_END;

		//insert exonic fragments to list
		newRangeJunc = new rangeJunction;
		newRangeJunc->junc = newFragment;

		sortList_Junction_Num++;
		if (sortList_Junction_Num >= sortList_Junction.size())
		{
			sortList_Junction.resize(sortList_Junction.size() + DEFAULT_JUNCTION_NUM, NULL);
			sortKey_Junction.resize(sortKey_Junction.size() + DEFAULT_JUNCTION_NUM, 0);
			mergeSort_Larray.resize(mergeSort_Larray.size() + DEFAULT_JUNCTION_NUM, 0);
			mergeSort_Rarray.resize(mergeSort_Rarray.size() + DEFAULT_JUNCTION_NUM, 0);
			mergeSort_LorderedList.resize(mergeSort_LorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
			mergeSort_RorderedList.resize(mergeSort_RorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
		}
		sortList_Junction[sortList_Junction_Num] = newRangeJunc;

		return;
	}

	sortList_spliceSite.resize(sortList_Junction_Num*2+10);
	sortKey_spliceSite.assign(sortList_Junction_Num*2+10, 0);
	mergeSort_LorderedList_spliceSite.resize(sortList_Junction_Num*2+10);
	mergeSort_RorderedList_spliceSite.resize(sortList_Junction_Num*2+10);

	//input splice sites
	sortList_spliceSite_Num = 0;
	for (i = 1; i <= sortList_Junction_Num; i++)
	{
		curSite = new spliceSite;		
		curSite->chromosome = sortList_Junction[i]->junc->chromosome_start;
		curSite->position = sortList_Junction[i]->junc->start;
		curSite->directionOut = true;
		sortList_spliceSite_Num++;
		sortList_spliceSite[sortList_spliceSite_Num] = curSite;

		curSite = new spliceSite;		
		curSite->chromosome = sortList_Junction[i]->junc->chromosome_end;
		curSite->position = sortList_Junction[i]->junc->end;
		curSite->directionOut = false;
		sortList_spliceSite_Num++;
		sortList_spliceSite[sortList_spliceSite_Num] = curSite;
	}


	//sort splice sites
	for (i = 1; i <= sortList_spliceSite_Num; i++)
	{
		sortKey_spliceSite[i] = sortList_spliceSite[i]->position;
	}
	mergeSort_SpliceSiteSort(sortList_spliceSite_Num);


	//build exonic fragments from adjacent splice sites 
	//exonic if not out-in

	//first exon
	newFragment = new fragment;
	newFragment->type = frag_exon;
	newFragment->altersite = 0;
	newFragment->frag_name = "Exon_1";
	newFragment->ID = ++fragmentID_Cnt;

	curSite = sortList_spliceSite[1];
	newFragment->chromosome_start = curSite->chromosome;
	newFragment->chromosome_end = curSite->chromosome;
	newFragment->start = CHROMOSOME_START;
	newFragment->end = curSite->position;

	//insert exonic fragments to list
	newRangeJunc = new rangeJunction;
	newRangeJunc->junc = newFragment;

	sortList_Junction_Num++;
	if (sortList_Junction_Num >= sortList_Junction.size())
	{
		sortList_Junction.resize(sortList_Junction.size() + DEFAULT_JUNCTION_NUM, NULL);
		sortKey_Junction.resize(sortKey_Junction.size() + DEFAULT_JUNCTION_NUM, 0);
		mergeSort_Larray.resize(mergeSort_Larray.size() + DEFAULT_JUNCTION_NUM, 0);
		mergeSort_Rarray.resize(mergeSort_Rarray.size() + DEFAULT_JUNCTION_NUM, 0);
		mergeSort_LorderedList.resize(mergeSort_LorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
		mergeSort_RorderedList.resize(mergeSort_RorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
	}
	sortList_Junction[sortList_Junction_Num] = newRangeJunc;

	for (i = 1; i < sortList_spliceSite_Num; i++)
	{
		if (sortList_spliceSite[i]->position == sortList_spliceSite[i+1]->position)
		{
			curSite = sortList_spliceSite[i];
			delete curSite;
		} 
		else
		{
			if (sortList_spliceSite[i]->directionOut == true && sortList_spliceSite[i+1]->directionOut == false)
			{
				curSite = sortList_spliceSite[i];
				delete curSite;
			} 
			else
			{
				newFragment = new fragment;
				newFragment->type = frag_exon;

				if (sortList_spliceSite[i]->directionOut == false && sortList_spliceSite[i+1]->directionOut == true)
					newFragment->altersite = 0;
				else if (sortList_spliceSite[i]->directionOut == true && sortList_spliceSite[i+1]->directionOut == true)
					newFragment->altersite = 1;
				else if (sortList_spliceSite[i]->directionOut == false && sortList_spliceSite[i+1]->directionOut == false)
					newFragment->altersite = 2;

				newFragment->frag_name = "Exon_" + itostr(i/2 + 2);
				newFragment->ID = ++fragmentID_Cnt;

				curSite = sortList_spliceSite[i];
				newFragment->chromosome_start = curSite->chromosome;
				if (curSite->directionOut == true)
				{
					newFragment->start = curSite->position + 1;
				} 
				else
				{
					newFragment->start = curSite->position;
				}
				delete curSite;

				curSite = sortList_spliceSite[i+1];
				newFragment->chromosome_end = curSite->chromosome;
				if (curSite->directionOut == true)
				{
					newFragment->end = curSite->position;
				} 
				else
				{
					newFragment->end = curSite->position - 1;
				}

				//insert exonic fragments to list
				newRangeJunc = new rangeJunction;
				newRangeJunc->junc = newFragment;

				sortList_Junction_Num++;
				if (sortList_Junction_Num >= sortList_Junction.size())
				{
					sortList_Junction.resize(sortList_Junction.size() + DEFAULT_JUNCTION_NUM, NULL);
					sortKey_Junction.resize(sortKey_Junction.size() + DEFAULT_JUNCTION_NUM, 0);
					mergeSort_Larray.resize(mergeSort_Larray.size() + DEFAULT_JUNCTION_NUM, 0);
					mergeSort_Rarray.resize(mergeSort_Rarray.size() + DEFAULT_JUNCTION_NUM, 0);
					mergeSort_LorderedList.resize(mergeSort_LorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
					mergeSort_RorderedList.resize(mergeSort_RorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
				}
				sortList_Junction[sortList_Junction_Num] = newRangeJunc;
			}
		}
	}

	//last exon
	newFragment = new fragment;
	newFragment->type = frag_exon;
	newFragment->altersite = 0;
	newFragment->frag_name = "Exon_" + itostr(sortList_spliceSite_Num/2 + 2);
	newFragment->ID = ++fragmentID_Cnt;

	curSite = sortList_spliceSite[sortList_spliceSite_Num];
	newFragment->chromosome_start = curSite->chromosome;
	newFragment->chromosome_end = curSite->chromosome;
	newFragment->start = curSite->position;
	newFragment->end = CHROMOSOME_END;
	delete curSite;

	//insert exonic fragments to list
	newRangeJunc = new rangeJunction;
	newRangeJunc->junc = newFragment;

	sortList_Junction_Num++;
	if (sortList_Junction_Num >= sortList_Junction.size())
	{
		sortList_Junction.resize(sortList_Junction.size() + DEFAULT_JUNCTION_NUM, NULL);
		sortKey_Junction.resize(sortKey_Junction.size() + DEFAULT_JUNCTION_NUM, 0);
		mergeSort_Larray.resize(mergeSort_Larray.size() + DEFAULT_JUNCTION_NUM, 0);
		mergeSort_Rarray.resize(mergeSort_Rarray.size() + DEFAULT_JUNCTION_NUM, 0);
		mergeSort_LorderedList.resize(mergeSort_LorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
		mergeSort_RorderedList.resize(mergeSort_RorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
	}
	sortList_Junction[sortList_Junction_Num] = newRangeJunc;

	sortList_spliceSite.clear(); free_vector(sortList_spliceSite);
	sortKey_spliceSite.clear(); free_vector(sortKey_spliceSite);
	mergeSort_LorderedList_spliceSite.clear(); free_vector(mergeSort_LorderedList_spliceSite);
	mergeSort_RorderedList_spliceSite.clear(); free_vector(mergeSort_RorderedList_spliceSite);

	return;
}

double windowChangeRatio(long *coverageVector, long vectorLength, long curPos)
{
	double ratio = 1.0;

	if (coverageVector == NULL)
		return 1.0;
	else if (curPos <= 0)
		return 2 * COVERAGE_CHANGE_THRESH;
	else if (curPos >= vectorLength - 1)
		return 0.0;

	double curSupp_left = 0, curSupp_right = 0, maxSupp_left, minSupp_left, maxSupp_right, minSupp_right; //support of the left window and the right window
	int windowSize_left = 0, windowSize_right = 0; //size of the left window and the right window 

	for (maxSupp_left = minSupp_left = coverageVector[curPos - 1]; windowSize_left < COVERAGE_CHANGE_WINDOW && curPos - windowSize_left > 0; ++windowSize_left)
	{
		curSupp_left += coverageVector[curPos - 1 - windowSize_left];
		if (coverageVector[curPos - 1 - windowSize_left] > maxSupp_left)
			maxSupp_left = coverageVector[curPos - 1 - windowSize_left];
		if (coverageVector[curPos - 1 - windowSize_left] < minSupp_left)
			minSupp_left = coverageVector[curPos - 1 - windowSize_left];
	}

	for (maxSupp_right = minSupp_right = coverageVector[curPos]; windowSize_right < COVERAGE_CHANGE_WINDOW && curPos + windowSize_right < vectorLength; ++windowSize_right)
	{
		curSupp_right += coverageVector[curPos + windowSize_right];
		if (coverageVector[curPos + windowSize_right] > maxSupp_right)
			maxSupp_right = coverageVector[curPos + windowSize_right];
		if (coverageVector[curPos + windowSize_right] < minSupp_right)
			minSupp_right = coverageVector[curPos + windowSize_right];
	}

	if (curSupp_left == 0.0)
		ratio = 2 * COVERAGE_CHANGE_THRESH;
	else
		ratio = (curSupp_right / windowSize_right) / (curSupp_left / windowSize_left);

	if (minSupp_left > maxSupp_right && ratio < 1)
		return ratio;
	else if (maxSupp_left < minSupp_right && ratio > 1)
		return ratio;
	else
		return 1.0;
}


double cutAlterSite(long *coverageVector, long vectorLength, long startIndex, int altersitetype, long &totalSupport, long &cutFragSupport, long &cutFragStart, long &cutFragEnd)
{
	double maxRatio = 0, supp_left, supp_right, curRatio;
	long size_left, size_right, curPos;
	totalSupport = 0; cutFragSupport = 0; cutFragStart = startIndex; cutFragEnd = startIndex + vectorLength - 1;
	
	if (coverageVector == NULL || altersitetype == 0)
		return 0;
	else if (altersitetype == 1)
	{
		for (size_left = 0, supp_left = 0.0; size_left < vectorLength - 1; ++size_left)
			supp_left += coverageVector[size_left];

		size_right = 1;
		supp_right = coverageVector[vectorLength - 1];

		totalSupport = long(supp_left + supp_right);

		for (curPos = vectorLength - 1; curPos > 0; --curPos)
		{
			if (supp_left == 0.0)
			{
				if (curPos < MAX_NOCOVERAGE_LENGTH)
				{
					//tolerate this gap
					return maxRatio>1? maxRatio : 0;
				} 
				else
				{
					//return this position
					cutFragStart = startIndex + curPos;
					cutFragSupport = supp_right;
					return COVERAGE_CHANGE_THRESH + 1; //definitely accept
				}				
			}
			else
				curRatio = (supp_right / size_right) / (supp_left / size_left);

			if (curRatio > maxRatio)
			{
				maxRatio = curRatio;
				cutFragStart = startIndex + curPos;				
				cutFragSupport = long(supp_right);
			}

			--size_left;
			supp_left -= coverageVector[curPos - 1];
			++size_right;
			supp_right += coverageVector[curPos - 1];
		}

		return maxRatio>1? maxRatio : 0;
	}
	else if (altersitetype == 2)
	{
		size_left = 0;
		supp_left = 0.0;

		for (size_right = 0, supp_right = 0.0; size_right < vectorLength; ++size_right)
			supp_right += coverageVector[size_right];

		totalSupport = long(supp_left + supp_right);

		for (curPos = 0; curPos < vectorLength; ++curPos)
		{
			if (supp_right == 0.0)
			{
				if (vectorLength - curPos <= MAX_NOCOVERAGE_LENGTH)
				{
					//tolerate this gap
					return maxRatio>1? maxRatio : 0;
				} 
				else
				{
					//return this position
					cutFragEnd = startIndex + curPos - 1;
					cutFragSupport = supp_left;
					return COVERAGE_CHANGE_THRESH + 1; //definitely accept
				}				
			}
			else if (size_left == 0)
				curRatio = 0;
			else
				curRatio = (supp_left / size_left) / (supp_right / size_right);
			
			if (curRatio > maxRatio)
			{
				maxRatio = curRatio;
				cutFragEnd = startIndex + curPos - 1;
				cutFragSupport = long(supp_left);
			}

			++size_left;
			supp_left += coverageVector[curPos];
			--size_right;
			supp_right -= coverageVector[curPos];
		}

		return maxRatio>1? maxRatio : 0;
	}
	else
	{
		cout << "Abnormal alter site type" << endl;
		exit(1);
	}
}


void getExonSupport(int index, long junctionListStart)
{
	//get junction support from fragent file
	ifstream exonSupportFile;
	string inputfilename;
	if (index < 0)
	{
#ifdef UNIX
		inputfilename = inputPath + "exon_all.txt";
#else
		inputfilename = "tmp\\exon_all.txt";	
#endif
		index = SUPPORT_VECTOR_SIZE;
	} 
	else
	{
#ifdef UNIX
#ifdef COVERAGE
		inputfilename = inputPath + "exon_" + itostr(index+1) + ".txt";
#else
		inputfilename = inputPath + "read_" + itostr(index+1) + ".txt";
#endif 
#else
#ifdef COVERAGE
		inputfilename = "tmp\\exon_" + itostr(index+1) + ".txt";
#else
		inputfilename = "tmp\\read_" + itostr(index+1) + ".txt";	
#endif
#endif
	}
	exonSupportFile.open(inputfilename.c_str());

	long startPosition, endPosition, juncIndex, tmp, tmpVectorLength, iLoop, jLoop, curFragStart, curSupport, read_retainCheck_start, read_retainCheck_end, fragAddCnt, tmpExonAddLength, altersiteFragStart, altersiteFragEnd, altersiteFragSupport;
	double changeRatio;
	fragment *newFragment;
	vector <fragment*> fragBuffer;
	fragBuffer.reserve(1000000);
	rangeJunction *newRangeJunc;
	deletedSites *newDelSite;
	string chromosome;
	int spliced;
	bool makingNewFrag, *read_inExon, doCheck, largeGap;
	startPosition = 0;
	juncIndex = junctionListStart;

	while (exonSupportFile >> chromosome)
	{
		exonSupportFile >> startPosition;
		exonSupportFile >> endPosition;
		endPosition--;
		exonSupportFile >> spliced;

		if (DORETAIN == true)
		{
			if (endPosition - startPosition + 1 + 1 > 10000)
			{
				cout << "error_" << index << "_" << startPosition << "_" << endPosition << " ";
			}
			read_inExon = new bool [endPosition - startPosition + 1 + 1]; //another 1 for sentinel
			for (iLoop = 0; iLoop <= endPosition - startPosition; iLoop++)
			{
				read_inExon[iLoop] = false;
			}
			read_inExon[endPosition - startPosition + 1] = true;
		}

		while (juncIndex <= sortList_Junction_Num && startPosition > sortList_Junction[juncIndex]->junc->end)
		{
			if (DOTRIMMING == true)// && juncIndex != junctionListStart && juncIndex != sortList_Junction_Num)
			{
				//check & separate
				fragAddCnt = 0;
				fragBuffer.clear();

				if (sortList_Junction[juncIndex]->junc->coverage != NULL)
				{
					curSupport = 0;
					makingNewFrag = false;
					curFragStart = sortList_Junction[juncIndex]->junc->start;
					tmpVectorLength = sortList_Junction[juncIndex]->junc->end - sortList_Junction[juncIndex]->junc->start + 1;

					if (sortList_Junction[juncIndex]->junc->type == frag_exon && sortList_Junction[juncIndex]->junc->altersite > 0)
					{
						changeRatio = cutAlterSite(sortList_Junction[juncIndex]->junc->coverage, tmpVectorLength, sortList_Junction[juncIndex]->junc->start, sortList_Junction[juncIndex]->junc->altersite, curSupport, altersiteFragSupport, altersiteFragStart, altersiteFragEnd);

						if (changeRatio < COVERAGE_CHANGE_THRESH)
						{
							//alternative splice site
							tmpExonAddLength = tmpVectorLength;
							if (tmpExonAddLength >= MIN_ALTER_SPLICE_SITE_LENGTH && curSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
							{
								//get an end, make a new fragment
								newFragment = new fragment;
								newFragment->type = sortList_Junction[juncIndex]->junc->type;
								newFragment->frag_name = "ExonAdd_" + itostr(fragmentID_Cnt + 1);
								newFragment->ID = ++fragmentID_Cnt;

								newFragment->chromosome_start = sortList_Junction[juncIndex]->junc->chromosome_start;
								newFragment->chromosome_end = sortList_Junction[juncIndex]->junc->chromosome_end;
								newFragment->start = sortList_Junction[juncIndex]->junc->start;
								newFragment->end = sortList_Junction[juncIndex]->junc->end;
								newFragment->support[index] += curSupport;

								++fragAddCnt;
								if (fragAddCnt >= fragBuffer.capacity())
									fragBuffer.reserve(fragBuffer.capacity() + 1000000);
								fragBuffer.push_back(newFragment); 
							}
						}
						else
						{
							//alternative start/end
							tmpExonAddLength = altersiteFragEnd - altersiteFragStart + 1;
							if (tmpExonAddLength >= MIN_EXON_LENGTH && altersiteFragSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
							{
								//get an end, make a new fragment
								newFragment = new fragment;
								newFragment->type = sortList_Junction[juncIndex]->junc->type;
								newFragment->frag_name = "ExonAdd_" + itostr(fragmentID_Cnt + 1);
								newFragment->ID = ++fragmentID_Cnt;

								newFragment->chromosome_start = sortList_Junction[juncIndex]->junc->chromosome_start;
								newFragment->chromosome_end = sortList_Junction[juncIndex]->junc->chromosome_end;
								newFragment->start = altersiteFragStart;
								newFragment->end = altersiteFragEnd;
								newFragment->support[index] += altersiteFragSupport;

								++fragAddCnt;
								if (fragAddCnt >= fragBuffer.capacity())
									fragBuffer.reserve(fragBuffer.capacity() + 1000000);
								fragBuffer.push_back(newFragment); 
							}
						}
					}
					else
					{

						for (iLoop = 0; iLoop <= tmpVectorLength; iLoop++)
						{
							if (iLoop + MAX_NOCOVERAGE_LENGTH - 1 > tmpVectorLength)
							{
								largeGap = true;
							} 
							else
							{
								largeGap = true;
								for (tmp = 0; tmp < MAX_NOCOVERAGE_LENGTH; tmp++)
								{
									if (sortList_Junction[juncIndex]->junc->coverage[iLoop + tmp] > coverageThreshold_exon)
									{
										largeGap = false;
										break;
									}
								}
							}

							changeRatio = windowChangeRatio(sortList_Junction[juncIndex]->junc->coverage, tmpVectorLength, iLoop);

							if (sortList_Junction[juncIndex]->junc->type == frag_exon)
							{
								if (sortList_Junction[juncIndex]->junc->coverage[iLoop] <= coverageThreshold_exon * SUPPORT_VECTOR_SIZE)
								{
									if (makingNewFrag == true && largeGap == true)
									{
										tmpExonAddLength = sortList_Junction[juncIndex]->junc->start + iLoop - curFragStart;
										if (tmpExonAddLength >= MIN_EXON_LENGTH && curSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
										{
											//get an end, make a new fragment
											newFragment = new fragment;
											newFragment->type = sortList_Junction[juncIndex]->junc->type;
											newFragment->frag_name = "ExonAdd_" + itostr(fragmentID_Cnt + 1);
											newFragment->ID = ++fragmentID_Cnt;

											newFragment->chromosome_start = sortList_Junction[juncIndex]->junc->chromosome_start;
											newFragment->chromosome_end = sortList_Junction[juncIndex]->junc->chromosome_end;
											newFragment->start = curFragStart;
											if (iLoop < tmpVectorLength)
												newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop;
											else
												newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop - 1; //sentinel case
											newFragment->support[index] += curSupport;

											++fragAddCnt;
											if (fragAddCnt >= fragBuffer.capacity())
												fragBuffer.reserve(fragBuffer.capacity() + 1000000);
											fragBuffer.push_back(newFragment); 
										}

										curSupport = 0;
										makingNewFrag = false;								
									}
									else
										curSupport += sortList_Junction[juncIndex]->junc->coverage[iLoop];
								}
								else
								{
									if (makingNewFrag == false)// && changeRatio > COVERAGE_CHANGE_THRESH)
										//if (makingNewFrag == false)
									{
										curFragStart = sortList_Junction[juncIndex]->junc->start + iLoop;
										makingNewFrag = true;
									}
									curSupport += sortList_Junction[juncIndex]->junc->coverage[iLoop];
								}
							}
							else if (sortList_Junction[juncIndex]->junc->type == frag_retained_intron)
							{
								if ((changeRatio < 1/COVERAGE_CHANGE_THRESH || sortList_Junction[juncIndex]->junc->coverage[iLoop] <= coverageThreshold_intron * SUPPORT_VECTOR_SIZE) && largeGap == true)
								{
									if (makingNewFrag == true)
									{
										tmpExonAddLength = sortList_Junction[juncIndex]->junc->start + iLoop - curFragStart;
										if (tmpExonAddLength >= MIN_EXON_LENGTH && curSupport >= coverageThreshold_intron * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
											//if (tmpExonAddLength >= MIN_EXON_LENGTH && tmpExonAddLength < tmpVectorLength - 1 && curSupport >= coverageThreshold_intron * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
										{
											//get an end, make a new fragment
											newFragment = new fragment;
											newFragment->type = sortList_Junction[juncIndex]->junc->type;
											newFragment->frag_name = "ExonAdd_" + itostr(fragmentID_Cnt + 1);
											newFragment->ID = ++fragmentID_Cnt;

											newFragment->chromosome_start = sortList_Junction[juncIndex]->junc->chromosome_start;
											newFragment->chromosome_end = sortList_Junction[juncIndex]->junc->chromosome_end;
											newFragment->start = curFragStart;
											if (iLoop < tmpVectorLength)
												newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop;
											else
												newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop - 1; //sentinel case
											newFragment->support[index] += curSupport;

											++fragAddCnt;
											if (fragAddCnt >= fragBuffer.capacity())
												fragBuffer.reserve(fragBuffer.capacity() + 1000000);
											fragBuffer.push_back(newFragment); 
										}

										curSupport = 0;
										makingNewFrag = false;								
									}
								}
								else
								{
									if (makingNewFrag == false && changeRatio > COVERAGE_CHANGE_THRESH) //for real data
									//if (makingNewFrag == false) //for simulated data
									{
										curFragStart = sortList_Junction[juncIndex]->junc->start + iLoop;
										makingNewFrag = true;
									}
									curSupport += sortList_Junction[juncIndex]->junc->coverage[iLoop];
								}
							}
						}
					}
				} 

				if (fragAddCnt > 1)
				{
					//shift
					sortList_Junction_Num += fragAddCnt - 1;
					if (sortList_Junction_Num >= sortList_Junction.size())
					{
						sortList_Junction.resize(sortList_Junction.size() + DEFAULT_JUNCTION_NUM, NULL);
						sortKey_Junction.resize(sortKey_Junction.size() + DEFAULT_JUNCTION_NUM, 0);
						mergeSort_Larray.resize(mergeSort_Larray.size() + DEFAULT_JUNCTION_NUM, 0);
						mergeSort_Rarray.resize(mergeSort_Rarray.size() + DEFAULT_JUNCTION_NUM, 0);
						mergeSort_LorderedList.resize(mergeSort_LorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
						mergeSort_RorderedList.resize(mergeSort_RorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
					}
					for (iLoop = sortList_Junction_Num; iLoop > juncIndex + fragAddCnt - 1; iLoop--)
					{
						sortList_Junction[iLoop] = sortList_Junction[iLoop - fragAddCnt + 1];
					}
				}

// 				if (juncIndex + fragAddCnt >= 9297)
// 				{
// 					cout << "error";
// 				}

				delete sortList_Junction[juncIndex]->junc;
				delete sortList_Junction[juncIndex];

				for (iLoop = 0; iLoop < fragAddCnt; iLoop++)
				{
					newRangeJunc = new rangeJunction;
					newRangeJunc->junc = fragBuffer[iLoop];

					sortList_Junction[juncIndex + iLoop] = newRangeJunc;
				}

				if (fragAddCnt == 0)
				{
					//no expression, the exon is deleted
					//shift back
					sortList_Junction_Num--;
					for (iLoop = juncIndex; iLoop <= sortList_Junction_Num; iLoop++)
					{
						sortList_Junction[iLoop] = sortList_Junction[iLoop + 1];
					}
				}

				juncIndex += fragAddCnt;
			}
			else
			{
				juncIndex++;
			}

		}

		if (juncIndex > sortList_Junction_Num)
		{
			break;
		}

		tmp = 0;
		while (juncIndex + tmp <= sortList_Junction_Num && startPosition <= sortList_Junction[juncIndex + tmp]->junc->end && endPosition >= sortList_Junction[juncIndex + tmp]->junc->start)
		{
			if (DOTRIMMING == true)
			{
				if (sortList_Junction[juncIndex + tmp]->junc->coverage == NULL)
				{
					tmpVectorLength = sortList_Junction[juncIndex + tmp]->junc->end - sortList_Junction[juncIndex + tmp]->junc->start + 1;
// 					if (tmpVectorLength > 1000000)
// 					{
// 						cout << "error";
// 					}
					sortList_Junction[juncIndex + tmp]->junc->coverage = new long [tmpVectorLength + 1]; //add one sentinel
					for (iLoop = 0; iLoop <= tmpVectorLength; iLoop++)
						sortList_Junction[juncIndex + tmp]->junc->coverage[iLoop] = 0;
				}
			}

			if (startPosition >= sortList_Junction[juncIndex + tmp]->junc->start)
			{
				if (endPosition <= sortList_Junction[juncIndex + tmp]->junc->end)
				{
					((sortList_Junction[juncIndex + tmp]->junc->support)[index]) += 1 * (endPosition - startPosition + 1);
					if (DOTRIMMING == true)
					{
						for (iLoop = startPosition - sortList_Junction[juncIndex + tmp]->junc->start; iLoop <= endPosition - sortList_Junction[juncIndex + tmp]->junc->start; iLoop++)
							sortList_Junction[juncIndex + tmp]->junc->coverage[iLoop] += 1;
					}
					if (DORETAIN == true)
					{
						for (iLoop = 0; iLoop <= endPosition - startPosition; iLoop++)
							read_inExon[iLoop] = true;
					}
				}
				else if (endPosition > sortList_Junction[juncIndex + tmp]->junc->end)
				{
					((sortList_Junction[juncIndex + tmp]->junc->support)[index]) += 1 * (sortList_Junction[juncIndex + tmp]->junc->end - startPosition + 1);
					if (DOTRIMMING == true)
					{
						for (iLoop = startPosition - sortList_Junction[juncIndex + tmp]->junc->start; iLoop <= sortList_Junction[juncIndex + tmp]->junc->end - sortList_Junction[juncIndex + tmp]->junc->start; iLoop++)
							sortList_Junction[juncIndex + tmp]->junc->coverage[iLoop] += 1;	
					}
					if (DORETAIN == true)
					{
						for (iLoop = 0; iLoop <= sortList_Junction[juncIndex + tmp]->junc->end - startPosition; iLoop++)
							read_inExon[iLoop] = true;
					}
				}
			} 
			else
			{
				if (endPosition <= sortList_Junction[juncIndex + tmp]->junc->end)
				{
					((sortList_Junction[juncIndex + tmp]->junc->support)[index]) += 1 * (endPosition - sortList_Junction[juncIndex + tmp]->junc->start + 1);
					if (DOTRIMMING == true)
					{
						for (iLoop = 0; iLoop <= endPosition - sortList_Junction[juncIndex + tmp]->junc->start; iLoop++)
							sortList_Junction[juncIndex + tmp]->junc->coverage[iLoop] += 1;
					}
					if (DORETAIN == true)
					{
						for (iLoop = sortList_Junction[juncIndex + tmp]->junc->start - startPosition; iLoop <= endPosition - startPosition; iLoop++)
							read_inExon[iLoop] = true;
					}
				}
				else if (endPosition > sortList_Junction[juncIndex + tmp]->junc->end)
				{
					((sortList_Junction[juncIndex + tmp]->junc->support)[index]) += 1 * (sortList_Junction[juncIndex + tmp]->junc->end - sortList_Junction[juncIndex + tmp]->junc->start + 1);
					if (DOTRIMMING == true)
					{
						for (iLoop = 0; iLoop <= sortList_Junction[juncIndex + tmp]->junc->end - sortList_Junction[juncIndex + tmp]->junc->start; iLoop++)
							sortList_Junction[juncIndex + tmp]->junc->coverage[iLoop] += 1;	
					}
					if (DORETAIN == true)
					{
						for (iLoop = sortList_Junction[juncIndex + tmp]->junc->start - startPosition; iLoop <= sortList_Junction[juncIndex + tmp]->junc->end - startPosition; iLoop++)
							read_inExon[iLoop] = true;
					}
				}
			}
			tmp++;
		}



		// 		if (tmp == 0)
		// 		{
		// 			//retained intron
		// 			retainedIntron << startPosition << '\t' << endPosition << endl;
		// 
		// 			sortList_Junction_Num++;
		// 			for (iLoop = sortList_Junction_Num; iLoop > juncIndex; iLoop--)
		// 			{
		// 				sortList_Junction[iLoop] = sortList_Junction[iLoop - 1];
		// 			}
		// 
		// 			//make an exon
		// 			newFragment = new fragment;
		// 			newFragment->type = frag_retained_intron;
		// 			sprintf(newFragment->name, "RetainedIntron_&ld", juncIndex);
		// 			newFragment->ID = ++fragmentID_Cnt;
		// 
		// 			strcpy(newFragment->chromosome_start, sortList_Junction[juncIndex-1]->junc->chromosome_start);
		// 			strcpy(newFragment->chromosome_end, sortList_Junction[juncIndex-1]->junc->chromosome_end);
		// 			newFragment->start = sortList_Junction[juncIndex-1]->junc->end + 1;
		// 			newFragment->end = sortList_Junction[juncIndex+1]->junc->start - 1;
		// 			newFragment->support[index] += 1 * (endPosition - startPosition);
		// 
		// 			//insert exonic fragments to list
		// 			newRangeJunc = new rangeJunction;
		// 			newRangeJunc->junc = newFragment;
		// 
		// 			sortList_Junction[juncIndex] = newRangeJunc;
		// 		}

// 		if (startPosition == 19117218 && endPosition == 19117317)
// 		{
// 			cout << "a";
// 		}

		if (DORETAIN == true)
		{
			doCheck = false;
			tmp = 0;
			for (iLoop = 0; iLoop <= endPosition - startPosition + 1; iLoop++)
			{
				if (read_inExon[iLoop] == true && doCheck == true)
				{
					//do check
					read_retainCheck_end = iLoop - 1;

					while (juncIndex + tmp <= sortList_Junction_Num)
					{
						if ((juncIndex + tmp - 1 < junctionListStart || juncIndex + tmp - 1 >= junctionListStart && startPosition + read_retainCheck_start > sortList_Junction[juncIndex + tmp - 1]->junc->end) && startPosition + read_retainCheck_end < sortList_Junction[juncIndex + tmp]->junc->start)
						{
							break;
						}
						tmp++;
					}
					

					//retained intron
//					retainedIntron << startPosition + read_retainCheck_start << '\t' << startPosition + read_retainCheck_end << endl;

					sortList_Junction_Num++;
					if (sortList_Junction_Num >= sortList_Junction.size())
					{
						sortList_Junction.resize(sortList_Junction.size() + DEFAULT_JUNCTION_NUM, NULL);
						sortKey_Junction.resize(sortKey_Junction.size() + DEFAULT_JUNCTION_NUM, 0);
						mergeSort_Larray.resize(mergeSort_Larray.size() + DEFAULT_JUNCTION_NUM, 0);
						mergeSort_Rarray.resize(mergeSort_Rarray.size() + DEFAULT_JUNCTION_NUM, 0);
						mergeSort_LorderedList.resize(mergeSort_LorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
						mergeSort_RorderedList.resize(mergeSort_RorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
					}
					for (jLoop = sortList_Junction_Num; jLoop > juncIndex + tmp; jLoop--)
					{
						sortList_Junction[jLoop] = sortList_Junction[jLoop - 1];
					}

					//make an exon
					newFragment = new fragment;
					newFragment->type = frag_retained_intron;
					newFragment->frag_name = "RetainedIntron_" + itostr(juncIndex);
					newFragment->ID = ++fragmentID_Cnt;

					newFragment->chromosome_start = sortList_Junction[juncIndex + tmp -1]->junc->chromosome_start;
					newFragment->chromosome_end = sortList_Junction[juncIndex + tmp -1]->junc->chromosome_end;
					if (juncIndex + tmp - 1 >= junctionListStart)
					{
						newFragment->start = sortList_Junction[juncIndex + tmp -1]->junc->end + 1;
						//newFragment->start = sortList_Junction[juncIndex + tmp -1]->junc->end;// + 1;
					} 
					else
					{
						newFragment->start = startPosition;
					}
					if (juncIndex + tmp + 1 <= sortList_Junction_Num)
					{
						newFragment->end = sortList_Junction[juncIndex + tmp +1]->junc->start - 1;
						//newFragment->end = sortList_Junction[juncIndex + tmp +1]->junc->start;// - 1;
					} 
					else
					{
						newFragment->end = endPosition;
					}
					newFragment->support[index] += 1 * (read_retainCheck_end - read_retainCheck_start + 1);

					tmpVectorLength = newFragment->end - newFragment->start + 1;
					newFragment->coverage = new long [tmpVectorLength + 1]; //add one sentinel
// 					if (tmpVectorLength > 1000000)
// 					{
// 						cout << "error";
// 					}
					for (jLoop = 0; jLoop <= tmpVectorLength; jLoop++)
						newFragment->coverage[jLoop] = 0;
					for (jLoop = startPosition + read_retainCheck_start - newFragment->start; jLoop <= startPosition + read_retainCheck_end - newFragment->start; jLoop++)
						newFragment->coverage[jLoop] += 1;


					//insert exonic fragments to list
					newRangeJunc = new rangeJunction;
					newRangeJunc->junc = newFragment;

					sortList_Junction[juncIndex + tmp] = newRangeJunc;


					doCheck = false;
				}
				else if (read_inExon[iLoop] == false && doCheck == false)
				{
					read_retainCheck_start = iLoop;
					doCheck = true;
				}
			}
		}


		if (DORETAIN == true)
		{
			delete [] read_inExon;
		}
	}

	exonSupportFile.close();

	while (juncIndex <= sortList_Junction_Num)
	{
		if (DOTRIMMING == true)// && juncIndex != junctionListStart && juncIndex != sortList_Junction_Num)
		{
////////////////////////////////////begin modification at 8/3/2012 to fix a boundary case
////////////////////////////////////the following the original code
// 			//check & separate
// 			fragAddCnt = 0;
// 
// 			if (sortList_Junction[juncIndex]->junc->coverage != NULL)
// 			{
// 				curSupport = 0;
// 				makingNewFrag = false;
// 				curFragStart = sortList_Junction[juncIndex]->junc->start;
// 				tmpVectorLength = sortList_Junction[juncIndex]->junc->end - sortList_Junction[juncIndex]->junc->start + 1;
// 				for (iLoop = 0; iLoop <= tmpVectorLength; iLoop++)
// 				{
// 					if (sortList_Junction[juncIndex]->junc->coverage[iLoop] <= coverageThreshold_exon)
// 					{
// 						// 							if (iLoop == 0)
// 						// 							{
// 						// 								//deleted 3' sites (Note: 5' for the exon, but 3' for the junction)
// 						// 								newDelSite = new deletedSites;
// 						// 								newDelSite->sites = sortList_Junction[juncIndex]->junc->start;
// 						// 							}
// 
// 						if (makingNewFrag == true)
// 						{
// 							tmpExonAddLength = sortList_Junction[juncIndex]->junc->start + iLoop - curFragStart;
// 							if (tmpExonAddLength >= MIN_EXON_LENGTH && curSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
// 							{
// 								//get an end, make a new fragment
// 								newFragment = new fragment;
// 								newFragment->type = sortList_Junction[juncIndex]->junc->type;
// 								sprintf(newFragment->name, "ExonAdd_%ld", fragmentID_Cnt + 1);
// 								newFragment->ID = ++fragmentID_Cnt;
// 
// 								strcpy(newFragment->chromosome_start, sortList_Junction[juncIndex]->junc->chromosome_start);
// 								strcpy(newFragment->chromosome_end, sortList_Junction[juncIndex]->junc->chromosome_end);
// 								newFragment->start = curFragStart;
// 								newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop - 1;
// 								newFragment->support[index] += curSupport;
// 
// 								fragBuffer[fragAddCnt++] = newFragment; 
// 							}
// 
// 							curSupport = 0;
// 							makingNewFrag = false;								
// 						}
// 					}
// 					else
// 					{
// 						if (makingNewFrag == false)
// 						{
// 							curFragStart = sortList_Junction[juncIndex]->junc->start + iLoop;
// 							makingNewFrag = true;
// 						}
// 						curSupport += sortList_Junction[juncIndex]->junc->coverage[iLoop];
// 					}
// 				}
// 			} 
////////////////////////////////////now the new code
			//check & separate
			fragAddCnt = 0;
			fragBuffer.clear();

			if (sortList_Junction[juncIndex]->junc->coverage != NULL)
			{
				curSupport = 0;
				makingNewFrag = false;
				curFragStart = sortList_Junction[juncIndex]->junc->start;
				tmpVectorLength = sortList_Junction[juncIndex]->junc->end - sortList_Junction[juncIndex]->junc->start + 1;

				if (sortList_Junction[juncIndex]->junc->type == frag_exon && sortList_Junction[juncIndex]->junc->altersite > 0)
				{
					changeRatio = cutAlterSite(sortList_Junction[juncIndex]->junc->coverage, tmpVectorLength, sortList_Junction[juncIndex]->junc->start, sortList_Junction[juncIndex]->junc->altersite, curSupport, altersiteFragSupport, altersiteFragStart, altersiteFragEnd);

					if (changeRatio < COVERAGE_CHANGE_THRESH)
					{
						//alternative splice site
						tmpExonAddLength = tmpVectorLength;
						if (tmpExonAddLength >= MIN_ALTER_SPLICE_SITE_LENGTH && curSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
						{
							//get an end, make a new fragment
							newFragment = new fragment;
							newFragment->type = sortList_Junction[juncIndex]->junc->type;
							newFragment->frag_name = "ExonAdd_" + itostr(fragmentID_Cnt + 1);
							newFragment->ID = ++fragmentID_Cnt;

							newFragment->chromosome_start = sortList_Junction[juncIndex]->junc->chromosome_start;
							newFragment->chromosome_end = sortList_Junction[juncIndex]->junc->chromosome_end;
							newFragment->start = sortList_Junction[juncIndex]->junc->start;
							newFragment->end = sortList_Junction[juncIndex]->junc->end;
							newFragment->support[index] += curSupport;

							++fragAddCnt;
							if (fragAddCnt >= fragBuffer.capacity())
								fragBuffer.reserve(fragBuffer.capacity() + 1000000);
							fragBuffer.push_back(newFragment); 
						}
					}
					else
					{
						//alternative start/end
						tmpExonAddLength = altersiteFragEnd - altersiteFragStart + 1;
						if (tmpExonAddLength >= MIN_EXON_LENGTH && altersiteFragSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
						{
							//get an end, make a new fragment
							newFragment = new fragment;
							newFragment->type = sortList_Junction[juncIndex]->junc->type;
							newFragment->frag_name = "ExonAdd_" + itostr(fragmentID_Cnt + 1);
							newFragment->ID = ++fragmentID_Cnt;

							newFragment->chromosome_start = sortList_Junction[juncIndex]->junc->chromosome_start;
							newFragment->chromosome_end = sortList_Junction[juncIndex]->junc->chromosome_end;
							newFragment->start = altersiteFragStart;
							newFragment->end = altersiteFragEnd;
							newFragment->support[index] += altersiteFragSupport;

							++fragAddCnt;
							if (fragAddCnt >= fragBuffer.capacity())
								fragBuffer.reserve(fragBuffer.capacity() + 1000000);
							fragBuffer.push_back(newFragment); 
						}
					}
				}
				else
				{

					for (iLoop = 0; iLoop <= tmpVectorLength; iLoop++)
					{
						if (iLoop + MAX_NOCOVERAGE_LENGTH - 1 > tmpVectorLength)
						{
							largeGap = true;
						} 
						else
						{
							largeGap = true;
							for (tmp = 0; tmp < MAX_NOCOVERAGE_LENGTH; tmp++)
							{
								if (sortList_Junction[juncIndex]->junc->coverage[iLoop + tmp] > coverageThreshold_exon)
								{
									largeGap = false;
									break;
								}
							}
						}

						changeRatio = windowChangeRatio(sortList_Junction[juncIndex]->junc->coverage, tmpVectorLength, iLoop);

						if (sortList_Junction[juncIndex]->junc->type == frag_exon)
						{
							if (sortList_Junction[juncIndex]->junc->coverage[iLoop] <= coverageThreshold_exon * SUPPORT_VECTOR_SIZE)
							{
								if (makingNewFrag == true && largeGap == true)
								{
									tmpExonAddLength = sortList_Junction[juncIndex]->junc->start + iLoop - curFragStart;
									if (tmpExonAddLength >= MIN_EXON_LENGTH && curSupport >= coverageThreshold_exon * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
									{
										//get an end, make a new fragment
										newFragment = new fragment;
										newFragment->type = sortList_Junction[juncIndex]->junc->type;
										newFragment->frag_name = "ExonAdd_" + itostr(fragmentID_Cnt + 1);
										newFragment->ID = ++fragmentID_Cnt;

										newFragment->chromosome_start = sortList_Junction[juncIndex]->junc->chromosome_start;
										newFragment->chromosome_end = sortList_Junction[juncIndex]->junc->chromosome_end;
										newFragment->start = curFragStart;
										if (iLoop < tmpVectorLength)
											newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop;
										else
											newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop - 1; //sentinel case
										newFragment->support[index] += curSupport;

										++fragAddCnt;
										if (fragAddCnt >= fragBuffer.capacity())
											fragBuffer.reserve(fragBuffer.capacity() + 1000000);
										fragBuffer.push_back(newFragment); 
									}

									curSupport = 0;
									makingNewFrag = false;								
								}
								else
									curSupport += sortList_Junction[juncIndex]->junc->coverage[iLoop];
							}
							else
							{
								if (makingNewFrag == false)// && changeRatio > COVERAGE_CHANGE_THRESH)
									//if (makingNewFrag == false)
								{
									curFragStart = sortList_Junction[juncIndex]->junc->start + iLoop;
									makingNewFrag = true;
								}
								curSupport += sortList_Junction[juncIndex]->junc->coverage[iLoop];
							}
						}
						else if (sortList_Junction[juncIndex]->junc->type == frag_retained_intron)
						{
							if ((changeRatio < 1/COVERAGE_CHANGE_THRESH || sortList_Junction[juncIndex]->junc->coverage[iLoop] <= coverageThreshold_intron * SUPPORT_VECTOR_SIZE) && largeGap == true)
							{
								if (makingNewFrag == true)
								{
									tmpExonAddLength = sortList_Junction[juncIndex]->junc->start + iLoop - curFragStart;
									if (tmpExonAddLength >= MIN_EXON_LENGTH && curSupport >= coverageThreshold_intron * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
										//if (tmpExonAddLength >= MIN_EXON_LENGTH && tmpExonAddLength < tmpVectorLength - 1 && curSupport >= coverageThreshold_intron * tmpExonAddLength * SUPPORT_VECTOR_SIZE)
									{
										//get an end, make a new fragment
										newFragment = new fragment;
										newFragment->type = sortList_Junction[juncIndex]->junc->type;
										newFragment->frag_name = "ExonAdd_" + itostr(fragmentID_Cnt + 1);
										newFragment->ID = ++fragmentID_Cnt;

										newFragment->chromosome_start = sortList_Junction[juncIndex]->junc->chromosome_start;
										newFragment->chromosome_end = sortList_Junction[juncIndex]->junc->chromosome_end;
										newFragment->start = curFragStart;
										if (iLoop < tmpVectorLength)
											newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop;
										else
											newFragment->end = sortList_Junction[juncIndex]->junc->start + iLoop - 1; //sentinel case
										newFragment->support[index] += curSupport;

										++fragAddCnt;
										if (fragAddCnt >= fragBuffer.capacity())
											fragBuffer.reserve(fragBuffer.capacity() + 1000000);
										fragBuffer.push_back(newFragment); 
									}

									curSupport = 0;
									makingNewFrag = false;								
								}
							}
							else
							{
								if (makingNewFrag == false && changeRatio > COVERAGE_CHANGE_THRESH)
								//if (makingNewFrag == false)
								{
									curFragStart = sortList_Junction[juncIndex]->junc->start + iLoop;
									makingNewFrag = true;
								}
								curSupport += sortList_Junction[juncIndex]->junc->coverage[iLoop];
							}
						}
					}
				}
			} 


			if (fragAddCnt > 1)
			{
				//shift
				sortList_Junction_Num += fragAddCnt - 1;
				if (sortList_Junction_Num >= sortList_Junction.size())
				{
					sortList_Junction.resize(sortList_Junction.size() + DEFAULT_JUNCTION_NUM, NULL);
					sortKey_Junction.resize(sortKey_Junction.size() + DEFAULT_JUNCTION_NUM, 0);
					mergeSort_Larray.resize(mergeSort_Larray.size() + DEFAULT_JUNCTION_NUM, 0);
					mergeSort_Rarray.resize(mergeSort_Rarray.size() + DEFAULT_JUNCTION_NUM, 0);
					mergeSort_LorderedList.resize(mergeSort_LorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
					mergeSort_RorderedList.resize(mergeSort_RorderedList.size() + DEFAULT_JUNCTION_NUM, NULL);
				}
				for (iLoop = sortList_Junction_Num; iLoop > juncIndex + fragAddCnt - 1; iLoop--)
				{
					sortList_Junction[iLoop] = sortList_Junction[iLoop - fragAddCnt + 1];
				}
			}

			// 				if (juncIndex + fragAddCnt >= 9297)
			// 				{
			// 					cout << "error";
			// 				}

			delete sortList_Junction[juncIndex]->junc;
			delete sortList_Junction[juncIndex];

			for (iLoop = 0; iLoop < fragAddCnt; iLoop++)
			{
				newRangeJunc = new rangeJunction;
				newRangeJunc->junc = fragBuffer[iLoop];

				sortList_Junction[juncIndex + iLoop] = newRangeJunc;
			}

			if (fragAddCnt == 0)
			{
				//no expression, the exon is deleted
				//shift back
				sortList_Junction_Num--;
				for (iLoop = juncIndex; iLoop <= sortList_Junction_Num; iLoop++)
				{
					sortList_Junction[iLoop] = sortList_Junction[iLoop + 1];
				}
			}

			juncIndex += fragAddCnt;
		}
		else
		{
			juncIndex++;
		}

	}

	fragBuffer.clear();
	free_vector(fragBuffer);

	return;
}


// void makeMergedJunc(fragment *targetJunc)
// {
// 	if (targetJunc->alter == NULL)
// 	{
// 		return;
// 	}
// 	
// 	alter_junction *curAlter;
// 	fragment *newFrag;
// 	
// 	//make a new junction for the targetJunc
// 	newFrag = new fragment;
// 	newFrag->
// 	
// 	
// 	return;
// }

void mergeAlterSpliceJunc()
{
	//merge all alternative splice site junctions
	//return a merged junction list
	long iLoop;
	alter_junction *newAlterJunc, *tmpAlterJunc, *alterList, *alterListTail;
	rangeJunction *resultList, *curJunc, *prevJunc, *curJunc_target, *curJunc_check, *prevJunc_check, *tmpJunc, *nextExon, *prevExon;
	fragment *mergedFrag;
	bool isAlter;
	double *alterSupport = new double [SUPPORT_VECTOR_SIZE];


	//////////////////////////////////////////////////////////////////////////
	//3' sites
	//////////////////////////////////////////////////////////////////////////
	for (iLoop = 1; iLoop <= sortList_Junction_Num; iLoop++)
	{
		sortKey_Junction[iLoop] = sortList_Junction[iLoop]->junc->start;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);

	for (iLoop = 1; iLoop <= sortList_Junction_Num; iLoop++)
	{
		sortKey_Junction[iLoop] = sortList_Junction[iLoop]->junc->end;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);

	//convert array to linked list;
	resultList = NULL;
	for (iLoop = 1; iLoop <= sortList_Junction_Num; iLoop++)
	{
		sortList_Junction[iLoop]->next = resultList;
		resultList = sortList_Junction[iLoop];
	}	

	prevJunc = NULL;
	curJunc = resultList;
	while (curJunc != NULL)
	{
		isAlter = false;

		if (curJunc->junc->type == frag_exon && curJunc->next != NULL && curJunc->junc->start != curJunc->junc->end)
		{
			//use exons as the trigger. exons will always come up first in such a sorted order

// 			if (curJunc->junc->start == 77938150 && curJunc->junc->end == 77938150)
// 			{
// 				cout << "a";
// 			}
			curJunc_target = curJunc->next;
			while (curJunc->junc->end == curJunc_target->junc->end)
			{
				alterList = NULL;
				alterListTail = NULL;
				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
				{
					alterSupport[iLoop] = 0.0;
				}

				prevJunc_check = curJunc;
				curJunc_check = curJunc->next;
				while (curJunc_check != NULL && curJunc_check->junc->end != curJunc->junc->start)
				{
					prevJunc_check = curJunc_check;
					curJunc_check = curJunc_check->next;
				}

				while (curJunc_check != NULL && curJunc_check->junc->end == curJunc->junc->start)
				{
					if (curJunc_check->junc->start == curJunc_target->junc->start)
					{
						if (alterList != NULL)
						{
							cout << "abnormal alter junction";
						}

						isAlter = true;

						newAlterJunc = new alter_junction;
						newAlterJunc->juncInfo = curJunc_check->junc;
//						newAlterJunc->category = 1;
						newAlterJunc->next = alterList;
						alterList = newAlterJunc;
						if (alterListTail == NULL)
						{
							alterListTail = newAlterJunc;
						}

						if (curJunc_check->junc->alter != NULL)
						{
							tmpAlterJunc = curJunc_check->junc->alter;
							while (tmpAlterJunc->next != NULL)
							{
								tmpAlterJunc = tmpAlterJunc->next;
							}

							tmpAlterJunc->next = alterList;
							alterList = curJunc_check->junc->alter;
							curJunc_check->junc->alter = NULL;
						}

						for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
						{
							alterSupport[iLoop] += curJunc_check->junc->support[iLoop];
						}

						//delete curJunc_check
						prevJunc_check->next = curJunc_check->next;
						curJunc_check->junc = NULL;
						delete curJunc_check;
						curJunc_check = prevJunc_check->next;
					}
					else
					{
						prevJunc_check = curJunc_check;
						curJunc_check = curJunc_check->next;
					}
				}


				//merge
				if (alterList != NULL)
				{
					if (curJunc_target->junc->alter == NULL)
					{
						//make
						mergedFrag = new fragment;
						mergedFrag->start = curJunc_target->junc->start;
						mergedFrag->end = curJunc_target->junc->end;
						mergedFrag->ID = ++fragmentID_Cnt;
						mergedFrag->chromosome_start = curJunc_target->junc->chromosome_start;
						mergedFrag->chromosome_end = curJunc_target->junc->chromosome_end;
						mergedFrag->type = frag_junction;

						newAlterJunc = new alter_junction;
						newAlterJunc->juncInfo = curJunc_target->junc;
//						newAlterJunc->category = 1;
						alterListTail->next = newAlterJunc;
						alterListTail = newAlterJunc;

						newAlterJunc = new alter_junction;
						newAlterJunc->juncInfo = curJunc->junc->clone();
//						newAlterJunc->category = 2;
						alterListTail->next = newAlterJunc;
						alterListTail = newAlterJunc;

						mergedFrag->alter = alterList;

						for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
						{
							alterSupport[iLoop] += curJunc_target->junc->support[iLoop];
							mergedFrag->support[iLoop] = alterSupport[iLoop];
						}

						//insert
						curJunc_target->junc = mergedFrag;
					}
					else
					{
						//use curJunc_target
						newAlterJunc = new alter_junction;
						newAlterJunc->juncInfo = curJunc->junc->clone();
//						newAlterJunc->category = 2;
						alterListTail->next = newAlterJunc;
						alterListTail = newAlterJunc;

						alterListTail->next = curJunc_target->junc->alter;
						curJunc_target->junc->alter = alterList;

						for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
						{
							alterSupport[iLoop] += curJunc_target->junc->support[iLoop];
							mergedFrag->support[iLoop] = alterSupport[iLoop];
						}
					}

					//.... 1. if curjunc_target->alter == null, then it's a real junc, merge to the alter list and make the mergedFrag
					//        otherwise it's a mergedFrag, either not making a new mergedFrag or delete the it
					//     2. if junctions in alter list have alter != null, then remember to merge their alters into the alter list!

				}

				curJunc_target = curJunc_target->next;
			}

		}

		if (isAlter == false)
		{
			prevJunc = curJunc;
			curJunc = curJunc->next;
		}
		else
		{
			//shift
			nextExon = NULL;
			prevExon = curJunc;
			tmpJunc = curJunc->next;
			while (tmpJunc != NULL && tmpJunc->junc->end != curJunc->junc->start)
			{
				prevExon = tmpJunc;
				tmpJunc = tmpJunc->next;
			}

			while (tmpJunc != NULL && tmpJunc->junc->end == curJunc->junc->start)
			{
				tmpJunc->junc->end_real = tmpJunc->junc->end;
				tmpJunc->junc->end = curJunc->junc->end; //shift the ends 

				if (tmpJunc->junc->type == frag_exon)
				{
					//save this exon and will move this exon to the next position
					nextExon = tmpJunc;
					prevExon->next = nextExon->next;

					tmpJunc = prevExon->next;
				}
				else
				{
					prevExon = tmpJunc;
					tmpJunc = tmpJunc->next;
				}
			}

			//delete
			if (nextExon == NULL)
			{
				prevJunc->next = curJunc->next;
			} 
			else
			{
				nextExon->next = curJunc->next;
				prevJunc->next = nextExon;
			}
			delete curJunc->junc;
			delete curJunc;
			curJunc = prevJunc->next;
		}
	}

	// 	for (iLoop = sortList_Junction_Num; iLoop > 0; iLoop--)
	// 	{
	// 		if (sortList_Junction[iLoop]->junc->type == frag_junction)
	// 		{
	// 			if (iLoop > 2 && sortList_Junction[iLoop - 1]->junc->type == frag_exon && sortList_Junction[iLoop]->junc->end == sortList_Junction[iLoop - 1]->junc->end
	// 				&& sortList_Junction[iLoop - 2]->junc->type == frag_junction && sortList_Junction[iLoop - 2]->junc->end == sortList_Junction[iLoop - 1]->junc->start
	// 				&& sortList_Junction[iLoop]->junc->start == sortList_Junction[iLoop - 2]->junc->start)
	// 			{
	// 				//got a 3' alternative, merge to iLoop-th junction and shift
	// 				
	// 				//merge
	// 				newAlterJunc = new alter_junction;
	// 				newAlterJunc->category = 1;
	// 				newAlterJunc->juncInfo = sortList_Junction[iLoop - 1]->junc;
	// 
	// 				//head-insert for 3' alternative
	// 				newAlterJunc->next = sortList_Junction[iLoop]->junc->alter;
	// 				sortList_Junction[iLoop]->junc->alter = newAlterJunc;
	// 
	// 				delete sortList_Junction[iLoop - 1];
	// 				
	// 				newAlterJunc = new alter_junction;
	// 				newAlterJunc->category = 1;
	// 				newAlterJunc->juncInfo = sortList_Junction[iLoop - 2]->junc;
	// 
	// 				//head-insert for 3' alternative
	// 				newAlterJunc->next = sortList_Junction[iLoop]->junc->alter;
	// 				sortList_Junction[iLoop]->junc->alter = newAlterJunc;
	// 
	// 				delete sortList_Junction[iLoop - 2];
	// 				
	// 				//shift
	// 				for (jLoop = iLoop; jLoop <= sortList_Junction_Num; jLoop++)
	// 				{
	// 					sortList_Junction[jLoop - 2] = sortList_Junction[jLoop];
	// 				}
	// 
	// 				sortList_Junction_Num -= 2;
	// 				iLoop -= 1; //examine iLoop-th one more time in case there are multiple alternative sites
	// 			}
	// 		}
	// 	}


	//convert back from junction list to array
	sortList_Junction_Num = 0;
	tmpJunc = resultList;
	while (tmpJunc != NULL)
	{
		sortList_Junction[++sortList_Junction_Num] = tmpJunc;
		tmpJunc = tmpJunc->next;
	}



	//5' sites
	for (iLoop = 1; iLoop <= sortList_Junction_Num; iLoop++)
	{
		sortKey_Junction[iLoop] = sortList_Junction[iLoop]->junc->end;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);

	for (iLoop = 1; iLoop <= sortList_Junction_Num; iLoop++)
	{
		sortKey_Junction[iLoop] = sortList_Junction[iLoop]->junc->start;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);


	//	return;

	//convert array to linked list;
	resultList = NULL;
	for (iLoop = sortList_Junction_Num; iLoop >= 1; iLoop--)
	{
		sortList_Junction[iLoop]->next = resultList;
		resultList = sortList_Junction[iLoop];
	}	

	prevJunc = NULL;
	curJunc = resultList;
	while (curJunc != NULL)
	{
		isAlter = false;

		if (curJunc->junc->type == frag_exon && curJunc->next != NULL && curJunc->junc->start != curJunc->junc->end)
		{
			//use exons as the trigger. exons will always come up first in such a sorted order

			curJunc_target = curJunc->next;
			while (curJunc->junc->start == curJunc_target->junc->start)
			{
				alterList = NULL;
				alterListTail = NULL;
				for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
				{
					alterSupport[iLoop] = 0.0;
				}

				prevJunc_check = curJunc;
				curJunc_check = curJunc->next;
				while (curJunc_check != NULL && curJunc_check->junc->start != curJunc->junc->end)
				{
					prevJunc_check = curJunc_check;
					curJunc_check = curJunc_check->next;
				}

				while (curJunc_check != NULL && curJunc_check->junc->start == curJunc->junc->end)
				{
					if (curJunc_check->junc->end == curJunc_target->junc->end)
					{
						if (alterList != NULL)
						{
							cout << "abnormal alter junction";
						}

						isAlter = true;

						//9/8/2011 if curJunc_check->junc->alter != NULL, do not import curJunc_check, just the alter list
						if (curJunc_check->junc->alter != NULL)
						{
							if (alterList == NULL)
							{
								alterList = curJunc_check->junc->alter;
							}

							tmpAlterJunc = curJunc_check->junc->alter;
							while (tmpAlterJunc->next != NULL)
							{
								tmpAlterJunc = tmpAlterJunc->next;
							}

							alterListTail = tmpAlterJunc;
							curJunc_check->junc->alter = NULL;
						}
						else
						{
							newAlterJunc = new alter_junction;
							newAlterJunc->juncInfo = curJunc_check->junc;
							//						newAlterJunc->category = -1;
							newAlterJunc->next = alterList;
							alterList = newAlterJunc;
							if (alterListTail == NULL)
							{
								alterListTail = newAlterJunc;
							}
						}

						for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
						{
							alterSupport[iLoop] += curJunc_check->junc->support[iLoop];
						}

						//delete curJunc_check
						prevJunc_check->next = curJunc_check->next;
						curJunc_check->junc = NULL;
						delete curJunc_check;
						curJunc_check = prevJunc_check->next;
					}
					else
					{
						prevJunc_check = curJunc_check;
						curJunc_check = curJunc_check->next;
					}
				}


				//merge
				if (alterList != NULL)
				{
					if (curJunc_target->junc->alter == NULL)
					{
						//make
						mergedFrag = new fragment;
						mergedFrag->start = curJunc_target->junc->start;
						mergedFrag->end = curJunc_target->junc->end;
						mergedFrag->ID = ++fragmentID_Cnt;
						mergedFrag->chromosome_start = curJunc_target->junc->chromosome_start;
						mergedFrag->chromosome_end = curJunc_target->junc->chromosome_end;
						mergedFrag->type = frag_junction;

						newAlterJunc = new alter_junction;
						newAlterJunc->juncInfo = curJunc_target->junc;
//						newAlterJunc->category = -1;
						newAlterJunc->next = alterList;
						alterList = newAlterJunc;

						newAlterJunc = new alter_junction;
						newAlterJunc->juncInfo = curJunc->junc->clone();
//						newAlterJunc->category = 2;
						newAlterJunc->next = alterList;
						alterList = newAlterJunc;

						mergedFrag->alter = alterList;

						for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
						{
							alterSupport[iLoop] += curJunc_target->junc->support[iLoop];
							mergedFrag->support[iLoop] = alterSupport[iLoop];
						}

						//insert
						curJunc_target->junc = mergedFrag;
					}
					else
					{
						//use curJunc_target
						newAlterJunc = new alter_junction;
						newAlterJunc->juncInfo = curJunc->junc->clone();
//						newAlterJunc->category = 2;
						newAlterJunc->next = alterList;
						alterList = newAlterJunc;

						alterListTail->next = curJunc_target->junc->alter;
						curJunc_target->junc->alter = alterList;

						for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
						{
							alterSupport[iLoop] += curJunc_target->junc->support[iLoop];
							mergedFrag->support[iLoop] = alterSupport[iLoop];
						}
					}

					//.... 1. if curjunc_target->alter == null, then it's a real junc, merge to the alter list and make the mergedFrag
					//        otherwise it's a mergedFrag, either not making a new mergedFrag or delete the it
					//     2. if junctions in alter list have alter != null, then remember to merge their alters into the alter list!

				}

				curJunc_target = curJunc_target->next;
			}

		}

		if (isAlter == false)
		{
			prevJunc = curJunc;
			curJunc = curJunc->next;
		}
		else
		{
			//shift
			nextExon = NULL;
			prevExon = curJunc;
			tmpJunc = curJunc->next;
			while (tmpJunc != NULL && tmpJunc->junc->start != curJunc->junc->end)
			{
				prevExon = tmpJunc;
				tmpJunc = tmpJunc->next;
			}

			while (tmpJunc != NULL && tmpJunc->junc->start == curJunc->junc->end)
			{
				tmpJunc->junc->start_real = tmpJunc->junc->start;
				tmpJunc->junc->start = curJunc->junc->start; //shift the ends 

				if (tmpJunc->junc->type == frag_exon)
				{
					//save this exon and will move this exon to the next position
					nextExon = tmpJunc;
					prevExon->next = nextExon->next;

					tmpJunc = prevExon->next;
				}
				else
				{
					prevExon = tmpJunc;
					tmpJunc = tmpJunc->next;
				}
			}

			//delete
			if (nextExon == NULL)
			{
				prevJunc->next = curJunc->next;
			} 
			else
			{
				nextExon->next = curJunc->next;
				prevJunc->next = nextExon;
			}
			delete curJunc->junc;
			delete curJunc;
			curJunc = prevJunc->next;
		}
	}

	// 	for (iLoop = 1; iLoop <= sortList_Junction_Num; iLoop++)
	// 	{
	// 		if (sortList_Junction[iLoop]->junc->type == frag_junction)
	// 		{
	// 			if (iLoop <= sortList_Junction_Num - 2 && sortList_Junction[iLoop + 1]->junc->type == frag_exon && sortList_Junction[iLoop]->junc->start == sortList_Junction[iLoop + 1]->junc->start
	// 				&& sortList_Junction[iLoop + 2]->junc->type == frag_junction && sortList_Junction[iLoop + 2]->junc->start == sortList_Junction[iLoop + 1]->junc->end
	// 				&& sortList_Junction[iLoop]->junc->end == sortList_Junction[iLoop + 2]->junc->end)
	// 			{
	// 				//got a 5' alternative, merge to iLoop-th junction and shift
	// 
	// 				//merge
	// 				newAlterJunc = new alter_junction;
	// 				newAlterJunc->category = -1;
	// 				newAlterJunc->juncInfo = sortList_Junction[iLoop + 1]->junc;
	// 
	// 				//tail-insert for 5' alternative
	// 				if (sortList_Junction[iLoop]->junc->fivePalterTail == NULL)
	// 				{
	// 					newAlterJunc->next = sortList_Junction[iLoop]->junc->alter;
	// 					sortList_Junction[iLoop]->junc->alter = newAlterJunc;
	// 					sortList_Junction[iLoop]->junc->fivePalterTail = newAlterJunc;
	// 				} 
	// 				else
	// 				{
	// 					newAlterJunc->next = sortList_Junction[iLoop]->junc->fivePalterTail->next;
	// 					sortList_Junction[iLoop]->junc->fivePalterTail->next = newAlterJunc;
	// 					sortList_Junction[iLoop]->junc->fivePalterTail = newAlterJunc;
	// 				}
	// 
	// 				delete sortList_Junction[iLoop + 1];
	// 
	// 				newAlterJunc = new alter_junction;
	// 				newAlterJunc->category = -1;
	// 				newAlterJunc->juncInfo = sortList_Junction[iLoop + 2]->junc;
	// 
	// 				//tail-insert for 5' alternative
	// 				newAlterJunc->next = sortList_Junction[iLoop]->junc->fivePalterTail->next;
	// 				sortList_Junction[iLoop]->junc->fivePalterTail->next = newAlterJunc;
	// 				sortList_Junction[iLoop]->junc->fivePalterTail = newAlterJunc;				
	// 
	// 				delete sortList_Junction[iLoop + 2];
	// 
	// 				//shift
	// 				for (jLoop = iLoop + 3; jLoop <= sortList_Junction_Num; jLoop++)
	// 				{
	// 					sortList_Junction[jLoop - 2] = sortList_Junction[jLoop];
	// 				}
	// 
	// 				sortList_Junction_Num -= 2;
	// 				iLoop--; //examine iLoop-th one more time in case there are multiple alternative sites
	// 			}
	// 		}
	// 	}

	//convert back from junction list to array
	sortList_Junction_Num = 0;
	tmpJunc = resultList;
	while (tmpJunc != NULL)
	{
		sortList_Junction[++sortList_Junction_Num] = tmpJunc;
		tmpJunc = tmpJunc->next;
	}

	delete [] alterSupport;

	return;
}


void cleanSpliceDirection()
{
	//merge splice junctions with same location but different strands
	if (CLEAN_SPLICE_DIRECTION == false || sortList_Junction_Num < 1)
		return;
	
	long iLoop;
	int indexLoop;
	double curSupp, nextSupp;

	rangeJunction *junctionList = NULL, *curJunc, *delJunc;
	for (iLoop = sortList_Junction_Num; iLoop >= 1; --iLoop)
	{
		sortList_Junction[iLoop]->next = junctionList;
		junctionList = sortList_Junction[iLoop];
	}

	curJunc = junctionList;
	while (curJunc->next != NULL)
	{
		if (curJunc->junc->start == curJunc->next->junc->start && curJunc->junc->end == curJunc->next->junc->end)
		{
			curSupp = 0.0; nextSupp = 0.0;
			for (indexLoop = 0; indexLoop < SUPPORT_VECTOR_SIZE; ++indexLoop)
			{
				curSupp += curJunc->junc->support[indexLoop];
				nextSupp += curJunc->next->junc->support[indexLoop];
			}

			curJunc->junc->transDirection = curSupp > nextSupp ? curJunc->junc->transDirection : curJunc->next->junc->transDirection;
			for (indexLoop = 0; indexLoop < SUPPORT_VECTOR_SIZE; ++indexLoop)
				curJunc->junc->support[indexLoop] += curJunc->next->junc->support[indexLoop];

			delJunc = curJunc->next;
			curJunc->next = curJunc->next->next;
			delete delJunc;
			--sortList_Junction_Num;
		}
		else
			curJunc = curJunc->next;
	}
		
	for (curJunc = junctionList, iLoop = 0; curJunc != NULL; curJunc = curJunc->next)
		sortList_Junction[++iLoop] = curJunc;

	if (iLoop != sortList_Junction_Num)
	{
		cout << "Number of splice junction does not match in cleanSpliceDirection" << endl;
		exit(1);
	}

	return;
}


void filterJunction(double thresh_maxValue, int thresh_zeroCnt, double thresh_meanValue)
{
	long iLoop, jLoop;
	int indexLoop, zeroCnt;
	bool keep;
	double total_support;

	for (iLoop = 1; iLoop <= sortList_Junction_Num; iLoop++)
	{
		keep = false;
		zeroCnt = 0;
		total_support = 0.0;

		for (indexLoop = 0; indexLoop < SUPPORT_VECTOR_SIZE; indexLoop++)
		{
			if (sortList_Junction[iLoop]->junc->support[indexLoop] > thresh_maxValue)
			{
				keep = true;
			}
			if (sortList_Junction[iLoop]->junc->support[indexLoop] == 0)
			{
				zeroCnt++;
			}
			total_support += sortList_Junction[iLoop]->junc->support[indexLoop];
		}
		if (keep == false || zeroCnt >= thresh_zeroCnt || total_support/SUPPORT_VECTOR_SIZE < thresh_meanValue)
		{
// 			gTreePlotTree << "0\t" << "1" << "\t" << sortList_Junction[iLoop]->junc->start << "\t" << sortList_Junction[iLoop]->junc->end << "\t" << "1\t0\t0\t0\t0";
// 			for (indexLoop = 0; indexLoop < SUPPORT_VECTOR_SIZE; indexLoop++)
// 				gTreePlotTree << "\t" << (sortList_Junction[iLoop]->junc->support)[indexLoop];
// 			gTreePlotTree << endl;

			//throw away
			delete sortList_Junction[iLoop]->junc;
			delete sortList_Junction[iLoop];

			sortList_Junction_Num--;
			for (jLoop = iLoop; jLoop <= sortList_Junction_Num; jLoop++)
			{
				sortList_Junction[jLoop] = sortList_Junction[jLoop + 1];
			}
			iLoop--;
		}
	}

	return;
}

void filterIntron(long start_sortList_Junction_index)
{
	long iLoop, jLoop;
	bool keep;
	int sampleLoop;
	double totalCoverage;

	for (iLoop = start_sortList_Junction_index; iLoop <= sortList_Junction_Num; ++iLoop)
	{
		if (sortList_Junction[iLoop]->junc->type == frag_retained_intron)
		{
			keep = false;
			totalCoverage = 0.0;
			for (sampleLoop = 0; sampleLoop < SUPPORT_VECTOR_SIZE; ++sampleLoop)
				totalCoverage += sortList_Junction[iLoop]->junc->support[sampleLoop];
			if (totalCoverage/SUPPORT_VECTOR_SIZE >= coverageThreshold_intron)
				keep = true;
		}
		else
			keep = true;
		
		if (keep == false)
		{
			//throw away
			delete sortList_Junction[iLoop]->junc;
			delete sortList_Junction[iLoop];

			sortList_Junction_Num--;
			for (jLoop = iLoop; jLoop <= sortList_Junction_Num; jLoop++)
			{
				sortList_Junction[jLoop] = sortList_Junction[jLoop + 1];
			}
			iLoop--;
		}
	}

	return;
}


void filterFalseFragments(bool filterStandaloneFragment)
{
	//filter standalone fragments (exons and junctions)
	if (sortList_Junction_Num <= 1)
	{
		//cout << "Warning: less than 2 fragments." << endl;
		return;
	}

	RangeJunctionList *resultList;
	rangeJunction *curJunc, *prevJunc, *delJunc, *tmpJunc;
	double maxSupport;
	long jLoop = 1;
	bool tobedeleted;

	//first, both ends of a junction must appear. filter false ones
	//filter junctions with no start exon
	for (jLoop = 1; jLoop <= sortList_Junction_Num; jLoop++)
	{
		sortKey_Junction[jLoop] = sortList_Junction[jLoop]->junc->end;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);

	resultList = new RangeJunctionList;	
	for (jLoop = 1; jLoop <= sortList_Junction_Num; jLoop++)
	{
		sortList_Junction[jLoop]->next = resultList->list;
		resultList->list = sortList_Junction[jLoop];
	}

	curJunc = resultList->list;
	prevJunc = NULL;
	while (curJunc != NULL)
	{
// 		if (curJunc->junc->start == 20889708)
// 		{
// 			cout << "a";
// 		}
		if (curJunc->junc->type == frag_junction)
		{
			tobedeleted = true;
			tmpJunc = curJunc->next;
			while (tmpJunc != NULL && tmpJunc->junc->end >= curJunc->junc->start)
			{
				if ((tmpJunc->junc->type == frag_exon || tmpJunc->junc->type == frag_retained_intron) && tmpJunc->junc->end == curJunc->junc->start)
				{
					tobedeleted = false;
					break;
				}
				tmpJunc = tmpJunc->next;
			}
			if (tobedeleted == true)
			{
				if (prevJunc != NULL)
				{
// 					gTreePlotTree << "0\t" << "1" << "\t" << curJunc->junc->start << "\t" << curJunc->junc->end << "\t" << "1\t0\t0\t0\t0";
// 					for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 						gTreePlotTree << "\t" << (curJunc->junc->support)[iLoop];
// 					gTreePlotTree << endl;

					prevJunc->next = curJunc->next;
					delete curJunc;
					sortList_Junction_Num--;
					curJunc = prevJunc->next;
				}
				else
				{
// 					gTreePlotTree << "0\t" << "1" << "\t" << curJunc->junc->start << "\t" << curJunc->junc->end << "\t" << "1\t0\t0\t0\t0";
// 					for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 						gTreePlotTree << "\t" << (curJunc->junc->support)[iLoop];
// 					gTreePlotTree << endl;

					resultList->list = curJunc->next;
					delete curJunc;
					sortList_Junction_Num--;
					curJunc = resultList->list;
				}
			}
			else
			{
				prevJunc = curJunc;
				curJunc = curJunc->next;
			}
		}
		else
		{
			prevJunc = curJunc;
			curJunc = curJunc->next;
		}
	}

	for (jLoop = 1, curJunc = resultList->list; curJunc != NULL; jLoop++, curJunc = curJunc->next)
	{
		sortList_Junction[jLoop] = curJunc;
	}
	if (jLoop != sortList_Junction_Num + 1)
	{
		cout << "Error: filter false junctions" << endl;
		exit(1);
	}
	resultList->list = NULL;

	//filter junctions with no end exon
	for (jLoop = 1; jLoop <= sortList_Junction_Num; jLoop++)
	{
		sortKey_Junction[jLoop] = sortList_Junction[jLoop]->junc->start;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);

	for (jLoop = sortList_Junction_Num; jLoop >= 1; jLoop--)
	{
		sortList_Junction[jLoop]->next = resultList->list;
		resultList->list = sortList_Junction[jLoop];
	}

	curJunc = resultList->list;
	prevJunc = NULL;
	while (curJunc != NULL)
	{
		if (curJunc->junc->type == frag_junction)
		{
			tobedeleted = true;
			tmpJunc = curJunc->next;
			while (tmpJunc != NULL && tmpJunc->junc->start <= curJunc->junc->end)
			{
				if ((tmpJunc->junc->type == frag_exon || tmpJunc->junc->type == frag_retained_intron) && tmpJunc->junc->start == curJunc->junc->end)
				{
					tobedeleted = false;
					break;
				}
				tmpJunc = tmpJunc->next;
			}
			if (tobedeleted == true)
			{
				if (prevJunc != NULL)
				{
// 					gTreePlotTree << "0\t" << "1" << "\t" << curJunc->junc->start << "\t" << curJunc->junc->end << "\t" << "1\t0\t0\t0\t0";
// 					for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 						gTreePlotTree << "\t" << (curJunc->junc->support)[iLoop];
// 					gTreePlotTree << endl;

					prevJunc->next = curJunc->next;
					delete curJunc;
					sortList_Junction_Num--;
					curJunc = prevJunc->next;
				}
				else
				{
// 					gTreePlotTree << "0\t" << "1" << "\t" << curJunc->junc->start << "\t" << curJunc->junc->end << "\t" << "1\t0\t0\t0\t0";
// 					for (iLoop = 0; iLoop < SUPPORT_VECTOR_SIZE; iLoop++)
// 						gTreePlotTree << "\t" << (curJunc->junc->support)[iLoop];
// 					gTreePlotTree << endl;

					resultList->list = curJunc->next;
					delete curJunc;
					sortList_Junction_Num--;
					curJunc = resultList->list;
				}
			}
			else
			{
				prevJunc = curJunc;
				curJunc = curJunc->next;
			}
		}
		else
		{
			prevJunc = curJunc;
			curJunc = curJunc->next;
		}
	}

	for (jLoop = 1, curJunc = resultList->list; curJunc != NULL; jLoop++, curJunc = curJunc->next)
	{
		sortList_Junction[jLoop] = curJunc;
	}
	if (jLoop != sortList_Junction_Num + 1)
	{
		cout << "Error: filter false junctions" << endl;
		exit(1);
	}

	resultList->list = NULL;


	//second, filter standalone & low-expression exons
	if (filterStandaloneFragment == true)
	{
		for (jLoop = 1; jLoop <= sortList_Junction_Num; jLoop++)
		{
			sortKey_Junction[jLoop] = sortList_Junction[jLoop]->junc->end;
		}
		mergeSort_JunctionSort(sortList_Junction_Num);

		for (jLoop = 1; jLoop <= sortList_Junction_Num; jLoop++)
		{
			sortKey_Junction[jLoop] = sortList_Junction[jLoop]->junc->start;
		}
		mergeSort_JunctionSort(sortList_Junction_Num);

		for (jLoop = sortList_Junction_Num; jLoop >= 1; jLoop--)
		{
			sortList_Junction[jLoop]->next = resultList->list;
			resultList->list = sortList_Junction[jLoop];
		}	

		jLoop = 1;
		curJunc = resultList->list;
		prevJunc = NULL;
		while (curJunc != NULL)
		{
			tobedeleted = false;
			delJunc = NULL;

			if (prevJunc != NULL && curJunc->next != NULL)
			{
				if (abs(prevJunc->junc->end	- curJunc->junc->start) <= 1 || abs(prevJunc->junc->start - curJunc->junc->start) <= 1 || abs(prevJunc->junc->end - curJunc->junc->end) <= 1 
					|| abs(curJunc->junc->end - curJunc->next->junc->start) <= 1 || abs(curJunc->junc->end - curJunc->next->junc->end) <= 1 || abs(curJunc->junc->start - curJunc->next->junc->start) <= 1)
				{//keep
					prevJunc = curJunc; 
					sortList_Junction[jLoop++] = curJunc;
					curJunc = curJunc->next;
				}
				else
				{
					prevJunc->next = curJunc->next;
					delJunc = curJunc;
					tobedeleted = true;
					sortList_Junction_Num--;
					curJunc = prevJunc->next;
				}
			}
			else if (prevJunc == NULL && curJunc->next != NULL)
			{
				if (abs(curJunc->junc->end - curJunc->next->junc->start) <= 1 || abs(curJunc->junc->end - curJunc->next->junc->end) <= 1 || abs(curJunc->junc->start - curJunc->next->junc->start) <= 1)
				{//keep
					prevJunc = curJunc; 
					sortList_Junction[jLoop++] = curJunc;
					curJunc = curJunc->next;
				}
				else
				{
					resultList->list = curJunc->next;
					delJunc = curJunc;
					tobedeleted = true;
					sortList_Junction_Num--;
					curJunc = resultList->list;
				}
			}
			else if (prevJunc != NULL && curJunc->next == NULL)
			{
				if (abs(prevJunc->junc->end	- curJunc->junc->start) <= 1 || abs(prevJunc->junc->start - curJunc->junc->start) <= 1 || abs(prevJunc->junc->end - curJunc->junc->end) <= 1)
				{//keep
					prevJunc = curJunc; 
					sortList_Junction[jLoop++] = curJunc;
					curJunc = curJunc->next;
				}
				else
				{
					prevJunc->next = curJunc->next;
					delJunc = curJunc;
					tobedeleted = true;
					sortList_Junction_Num--;
					curJunc = prevJunc->next;
				}
			}
			else
			{
				//only one fragment, filter by coverage
// 				if (maxSupport > MIN_COVERAGE_FILTER_FRAGMENT)
// 				{//keep
// 					prevJunc = curJunc; 
// 					sortList_Junction[jLoop++] = curJunc;
// 					curJunc = curJunc->next;
// 				}
//				else
				{
					resultList->list = curJunc->next;
					delete curJunc;
					sortList_Junction_Num--;
					curJunc = resultList->list;
				}
			}

			if (tobedeleted == true && delJunc != NULL)
			{
				delete delJunc;
			}
		}

		for (jLoop = 1, curJunc = resultList->list; curJunc != NULL; jLoop++, curJunc = curJunc->next)
		{
			sortList_Junction[jLoop] = curJunc;
		}
		if (jLoop != sortList_Junction_Num + 1)
		{
			cout << "Error: filter false exons" << endl;
			exit(1);
		}
	}

	resultList->list = NULL;
	delete resultList;

	return;
}


RangeJunctionList* buildOrigList()
{
	RangeJunctionList *resultList;
	alter_junction *curAlterJunc;
	long i;
	int alterJuncCnt;

//#ifdef FILTER_FRAGMENTS
	filterFalseFragments(false);
//#endif

#ifdef DO_MERGEALTERSITES
	mergeAlterSpliceJunc();

	//output alternative splice sites to leaf file
	for (i = 1; i <= sortList_Junction_Num; i++)
	{
		if (sortList_Junction[i]->junc->alter != NULL)
		{
			alterJuncCnt = 0;
			sortList_Junction[i]->junc->alterFragCnt = 0;
			curAlterJunc = sortList_Junction[i]->junc->alter;
			while (curAlterJunc != NULL)
			{
				(sortList_Junction[i]->junc->alterFragCnt)++;
				if (curAlterJunc->juncInfo->type == frag_junction)
				{
					alterJuncCnt++;
				}
				curAlterJunc = curAlterJunc->next;
			}
		}
	}
#else
	for (i = 1; i <= sortList_Junction_Num; i++)
	{
		sortKey_Junction[i] = sortList_Junction[i]->junc->end;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);
 
	for (i = 1; i <= sortList_Junction_Num; i++)
	{
		sortKey_Junction[i] = sortList_Junction[i]->junc->start;
	}
	mergeSort_JunctionSort(sortList_Junction_Num);
#endif
	
	resultList = new RangeJunctionList;
	resultList->rangeLow = CHROMOSOME_START;
	resultList->rangeHigh = CHROMOSOME_END;

	for (i = sortList_Junction_Num; i >= 1; i--)
	{
		sortList_Junction[i]->next = resultList->list;
		resultList->list = sortList_Junction[i];
	}	

	return resultList;
}


void output_alterJunc(GTvertex *alterSite, ofstream *outputfile)
{
	//output alternative splice sites
	GTedge *childedge;
	rangeJunction *curJunc;
	int tmp;

	childedge = alterSite->child;
	while (childedge != NULL)
	{
		(*outputfile) << endl;
		for (tmp = 1; tmp <= alterSite->level + 7; tmp++)
		{
			(*outputfile) << "  ";
		}
		(*outputfile) << "| ";

		for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
		{
			(*outputfile) << childedge->linkedVertex->proportion[tmp] << ";";
		}
		(*outputfile) << "=  ";		

		curJunc = childedge->linkedVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			(*outputfile) << "(";
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
			{
				(*outputfile) << "exon: ";
			} 
			else
			{
				(*outputfile) << "junc: ";
			}
			(*outputfile) << curJunc->junc->start << ", " << curJunc->junc->end << "; ";
			for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
				(*outputfile) << (curJunc->junc->support)[tmp] << "/ ";
			(*outputfile) << ") ";

			curJunc = curJunc->next;
		}

		childedge = childedge->next;
	}

	return;
}

void calculate_ASM_group_meanExpression(GTvertex *rootVertex)
{
	double ASMsupport1, ASMsupport2;
	int tmp;

	ASMsupport1 = 0.0;
	ASMsupport2 = 0.0;
	for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE/2; ++tmp)
	{
		ASMsupport1 += rootVertex->support[tmp];
		ASMsupport2 += rootVertex->support[tmp + SUPPORT_VECTOR_SIZE/2];
	}
	ASMsupport1 = ASMsupport1 / (SUPPORT_VECTOR_SIZE/2.);
	ASMsupport2 = ASMsupport2 / (SUPPORT_VECTOR_SIZE/2.);
	
	rootVertex->ASMsupport_group1 = ASMsupport1;
	rootVertex->ASMsupport_group2 = ASMsupport2;

	return;		
}



bool GTree_output(GTvertex *rootVertex, int sign, ofstream *outputfile)
{
	//output GTree in pre-order
	//return true if still in a minimum ASM, i.e., still no type 2 or 3 nodes
	GTvertex *curVertex;
	GTedge *curEdge;
	rangeJunction *curJunc;
	long i;
	int vertexCategory, tmp;
	bool isMinASM, tmpFlag;
	

	for (i = 1; i <= rootVertex->level; i++)
	{
		(*outputfile) << "  ";
	}

	if (rootVertex->level > 0)
	{
		if (sign == 1)
		{
			//children are independent regions
			(*outputfile) << "-";
		} 
		else if (sign == 2)
		{
			//children are independent paths
			(*outputfile) << "+";
		}
		else if (sign == 3)
		{
			//children are dependent paths
			(*outputfile) << "x";
		}
	}

	(*outputfile) << "[" << rootVertex->rangeLow << ", " << rootVertex->rangeHigh << "] " << rootVertex->childNum << '/' << rootVertex->junctionNum << " ";


	vertexCategory = alterSpliceCate(rootVertex);

	if (vertexCategory == exon_skipping)
	{
		exon_skipping_cnt++;
		(*outputfile) << "exon_skipping ";
	}
	else if (vertexCategory == mutual_exclusive)
	{
		mutual_exclusive_cnt++;
		(*outputfile) << "mutual_exclusive ";
	}
	else if (vertexCategory == intron_retention)
	{
		intron_retention_cnt++;
		(*outputfile) << "retained_intron ";
	}
	else if (vertexCategory == diff_expression)
	{
		(*outputfile) << "diff_expression ";
	}


	if (rootVertex->child == NULL)
	{
		(*outputfile) << "{";
		for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
		{
			(*outputfile) << rootVertex->support[tmp]  << "/" << rootVertex->proportion[tmp] << "; ";
		}
		(*outputfile) << "}";

		
		//print fragment list
		(*outputfile) << ": ";
		curJunc = rootVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			(*outputfile) << "(";
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
			{
				(*outputfile) << "exon: ";
			} 
			else
			{
				(*outputfile) << "junc: ";
			}
			(*outputfile) << curJunc->junc->start << ", " << curJunc->junc->end << "; ";
			for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
				(*outputfile) << (curJunc->junc->support)[tmp] << "/ ";
			(*outputfile) << ") ";

			curJunc = curJunc->next;
		}



		// 		if (rootVertex->anovaScore_support > 3.46)
		// 		{
		// 			(*outputfile) << "DIFF";
		// 		}
		(*outputfile) << endl;

		if (rootVertex->alterSpliceSite == NULL)
		{
			
		}
		else
		{
			//output alternative splice sites

			curEdge = rootVertex->alterSpliceSite;
			while (curEdge != NULL)
			{
				curVertex = curEdge->linkedVertex;
				for (i = 1; i <= curVertex->level; i++)
				{
					(*outputfile) << "  ";
				}
				(*outputfile) << "  alter: " << curVertex->childNum << ". ";
				output_alterJunc(curVertex, outputfile);
				(*outputfile) << endl;

				if (curVertex->ASMcategory > 0 && rootVertex->ASMcategory != diff_expression)
				{
					alterSpliceCate(curVertex);
				}

				curEdge = curEdge->next;
			}
		}


// 		if (rootVertex->junctionNum > 1 && sign != 3)
// 		{
// 			curJunc = rootVertex->junctionInRange->list;
// 			while (curJunc != NULL)
// 			{
// 				if (curJunc->junc->type == frag_exon)
// 				{
// 					abnormalfile << "1\t";
// 				} 
// 				else
// 				{
// 					abnormalfile << "0\t";
// 				}
// 				abnormalfile << curJunc->junc->start << "\t" << curJunc->junc->end << endl;
// 				curJunc = curJunc->next;
// 			}
// 			cout << "error: rootVertex->junctionNum > 1 && sign != 3" << endl;
// 			exit(1);
// 			//abnormalfile << endl << endl << endl;
// 		}
		
		return true;
	} 
	else
	{
		//extend children
		// 		(*outputfile) << "[";
		// 		for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
		// 		{
		// 			(*outputfile) << rootVertex->representative->support[tmp] << "/ ";
		// 		}
		// 		(*outputfile) << "]  ";

		(*outputfile) << "{";
		for (tmp = 0; tmp < SUPPORT_VECTOR_SIZE; tmp++)
		{
			(*outputfile) << rootVertex->support[tmp] << "/" << rootVertex->proportion[tmp] << "; ";
		}
		(*outputfile) << "}";

		if (rootVertex->ASMcategory != diff_expression)
		{
			(*outputfile) << endl;
		}
		

		isMinASM = true;

		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;
			tmpFlag = GTree_output(curVertex, rootVertex->childType, outputfile);

			if (tmpFlag == false)
			{
				isMinASM = false;
			}

			curEdge = curEdge->next;
		}


//		if (isMinASM == true)
//		{
			if (rootVertex->childType == 2 || rootVertex->childType == 3)
			{
				return false;
			} 
			else
			{
				return true;
			}
//		} 
//		else
//		{
//			return false;
//		}

		//output alternative splice sites
		if (rootVertex->alterSpliceSite != NULL)
		{
			curEdge = rootVertex->alterSpliceSite;
			while (curEdge != NULL)
			{
				curVertex = curEdge->linkedVertex;
				for (i = 1; i <= curVertex->level; i++)
				{
					(*outputfile) << "  ";
				}
				(*outputfile) << "  alter: " << curVertex->childNum << ". ";
				output_alterJunc(curVertex, outputfile);
				(*outputfile) << endl;

				if (curVertex->ASMcategory > 0 && rootVertex->ASMcategory != diff_expression)
				{
					alterSpliceCate(curVertex);
				}
								
				curEdge = curEdge->next;
			}
		}
	}
}




/************************************************************************/
/* DIFFERENTIAL ANALYSIS                                             */
/************************************************************************/

void get_gene_expression_analysis_list(GTvertex *rootVertex)
{
	GTedge *curEdge;
	int groupLoopCnt, individualLoopCnt, techRepLoopCnt, sampleSizeCnt, supportIndex;


	if (rootVertex->level < 1 && rootVertex != gTree->root && rootVertex->exonNum >= 1)
		vertexListForStatistics.push_back(rootVertex);
	
	if (rootVertex->child == NULL)
	{
		//do nothing
	}
	else
	{
		//extend children
		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			get_gene_expression_analysis_list(curEdge->linkedVertex);
			curEdge = curEdge->next;
		}
	}
}


void differential_analysis_expression_level()
{
	unsigned long geneLoopCnt;
	GTvertex *curVertex;
	int sampleLoopCnt;

	//get expression matrix
	vertexListForStatistics.clear();
	get_gene_expression_analysis_list(gTree->root);

	for (geneLoopCnt = 0; geneLoopCnt < vertexListForStatistics.size(); ++geneLoopCnt)
	{
		curVertex = vertexListForStatistics[geneLoopCnt];
		outfile_stats_expr << resultNamePrefix << "\t" << curVertex->rangeLow << "\t" << curVertex->rangeHigh;
		for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; ++sampleLoopCnt)
		{
			outfile_stats_expr << "\t" << curVertex->support[sampleLoopCnt];
		}
		for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; ++sampleLoopCnt)
		{
			outfile_stats_expr << "\t" << curVertex->proportion[sampleLoopCnt];
		}
		outfile_stats_expr << endl;
	}

	return;
}



int asm_category(GTvertex *curVertex)
{
	//categorize alternative splicing events, structure only, no expression

	//if the ASM has been assigned a category, then use that
	if (curVertex->ASMcategory > 0)
	{
		return curVertex->ASMcategory;
	}

	if ((curVertex->childType == 2 || curVertex->childType == 3) && curVertex->childNum > 1 && curVertex->major_alter_paths_num > 1)
	{
		//ASM
	}
	else
	{
		return -1;
	}

	//categorize the alternative splicing
	alternative_path *pathA, *pathB;

	if (curVertex->childType == 2)
	{
		pathA = curVertex->major_alter_paths;
		while (pathA != NULL)
		{
			pathB = pathA->next;
			while (pathB != NULL)
			{
				if (abs(pathA->path_start - pathB->path_start) < 2 && abs(pathA->path_end - pathB->path_end) < 2)
				{
					if (pathA->exonNum == 0 && pathA->junctionNum == 1 && pathB->junctionNum >= 2 && pathB->exonNum >= 1 || pathA->junctionNum >= 2 && pathA->exonNum >= 1 && pathB->exonNum == 0 && pathB->junctionNum == 1)
					{
						curVertex->ASMcategory = exon_skipping;
						return exon_skipping;
					}
					if (pathA->exonNum >= 1 && pathA->junctionNum >= 2 && pathB->exonNum >= 1 && pathB->junctionNum >= 2)
					{
						curVertex->ASMcategory = mutual_exclusive;						
					}
					if (curVertex->major_alter_paths_num == 2 && (pathA->exonNum == 1 && pathA->junctionNum == 0 && pathB->exonNum == 0 && pathB->junctionNum == 1 || pathA->exonNum == 0 && pathA->junctionNum == 1 && pathB->exonNum == 1 && pathB->junctionNum == 0))
					{
						curVertex->ASMcategory = intron_retention;
						return intron_retention;
					}
				}
				else if (abs(pathA->path_start - pathB->path_start) < 2)
				{
					curVertex->ASMcategory = alter_end;	
				}
				else if (abs(pathA->path_end - pathB->path_end) < 2)
				{
					curVertex->ASMcategory = alter_start;	
				}

				pathB = pathB->next;
			}

			pathA = pathA->next;
		}
	}

	if (curVertex->ASMcategory == mutual_exclusive || curVertex->ASMcategory == alter_end || curVertex->ASMcategory == alter_start)
	{
		return curVertex->ASMcategory;
	}

	curVertex->ASMcategory = unknown;
	return unknown;
}


void output_ASMpath_gtf(GTvertex *targetVertex, RangeJunctionList *pathList)
{
	long path_cnt = 0, exon_cnt;

	RangeJunctionList *curPath;
	rangeJunction *curJunc;
	char strand;

	curPath = pathList;
	while (curPath != NULL)
	{
		++path_cnt;
		exon_cnt = 0;

		if (curPath->transDirection == antisense)
			strand = '-';
		else
			strand = '+';

		if (targetVertex->prevSibling != NULL && abs(curPath->rangeLow - targetVertex->prevSibling->rangeHigh) < 2)
			outfile_gtf_asm_path << resultNamePrefix << "\t" << "ASM" << "\t" << "exon" << "\t" << targetVertex->prevSibling->rangeLow << "\t" << targetVertex->prevSibling->rangeHigh << "\t" << "." << "\t" << strand << "\t" << "." << "\t"
			<< "gene_id \"" << resultNamePrefix << "." << targetVertex->ID << "\"; " << "transcript_id \"" << resultNamePrefix << "." << targetVertex->ID << ".p" << path_cnt << "\";" << endl;

		curJunc = curPath->list;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron)
				outfile_gtf_asm_path << resultNamePrefix << "\t" << "ASM" << "\t" << "exon" << "\t" << curJunc->junc->start << "\t" << curJunc->junc->end << "\t" << "." << "\t" << strand << "\t" << "." << "\t"
				<< "gene_id \"" << resultNamePrefix << "." << targetVertex->ID << "\"; " << "transcript_id \"" << resultNamePrefix << "." << targetVertex->ID << ".p" << path_cnt << "\";" << endl;

			curJunc = curJunc->next;
		}

		if (targetVertex->nextSibling != NULL && abs(curPath->rangeHigh - targetVertex->nextSibling->rangeLow) < 2)
			outfile_gtf_asm_path << resultNamePrefix << "\t" << "ASM" << "\t" << "exon" << "\t" << targetVertex->nextSibling->rangeLow << "\t" << targetVertex->nextSibling->rangeHigh << "\t" << "." << "\t" << strand << "\t" << "." << "\t"
			<< "gene_id \"" << resultNamePrefix << "." << targetVertex->ID << "\"; " << "transcript_id \"" << resultNamePrefix << "." << targetVertex->ID << ".p" << path_cnt << "\";" << endl;


		curPath = curPath->nextList;
	}
	
	return;
}

// void output_asm_path(GTvertex *targetVertex, ofstream *targetFile, char *ASM_id)
// {
// 	alternative_path *asm_path;
// 	GTedge *path_component;
// 	GTvertex *cur_component;
// 	long path_cnt = 0, exon_cnt;
// 
// 	asm_path = targetVertex->major_alter_paths;
// 	while (asm_path != NULL)
// 	{
// 		++path_cnt;
// 		exon_cnt = 0;
// 
// 		if (targetVertex->prevSibling != NULL && targetVertex->prevSibling->exonNum > 0 && abs(asm_path->path_start - targetVertex->prevSibling->rangeHigh) < 2)
// 			*targetFile << resultNamePrefix << "\t" << "ASM" << "\t" << "exon" << "\t" << targetVertex->prevSibling->rangeLow << "\t" << targetVertex->prevSibling->rangeHigh << "\t" << "." << "\t" << "+\t" << "." << "\t"
// 				<< "gene_id \"" << ASM_id << "\"; " << "transcript_id \"" << ASM_id << "_p" << path_cnt << "\";" << endl;
// 
// 		path_component = asm_path->pathVertex->child;
// 		while (path_component != NULL)
// 		{
// 			cur_component = path_component->linkedVertex;
// 			if (cur_component->exonNum > 0)
// 				*targetFile << resultNamePrefix << "\t" << "ASM" << "\t" << "exon" << "\t" << cur_component->rangeLow << "\t" << cur_component->rangeHigh << "\t" << "." << "\t" << "+\t" << "." << "\t"
// 				<< "gene_id \"" << ASM_id << "\"; " << "transcript_id \"" << ASM_id << "_p" << path_cnt << "\";" << endl;
// 
// 			path_component = path_component->next;
// 		}
// 
// 		if (targetVertex->nextSibling != NULL && targetVertex->nextSibling->exonNum > 0 && abs(asm_path->path_end - targetVertex->nextSibling->rangeLow) < 2)
// 			*targetFile << resultNamePrefix << "\t" << "ASM" << "\t" << "exon" << "\t" << targetVertex->nextSibling->rangeLow << "\t" << targetVertex->nextSibling->rangeHigh << "\t" << "." << "\t" << "+\t" << "." << "\t"
// 				<< "gene_id \"" << ASM_id << "\"; " << "transcript_id \"" << ASM_id << "_p" << path_cnt << "\";" << endl;
// 
// 		asm_path = asm_path->next;
// 	}
// 
// 	return;
// }

void output_one_asm(GTvertex *targetVertex)
{
	int sampleLoopCnt;
	alternative_path *curAlterPath;
	string ASM_id;

	ASM_id = resultNamePrefix + "." + itostr(targetVertex->ID);

	//output vertex information
	outfile_stat_asm << ASM_id << "\t" << resultNamePrefix << "\t" << targetVertex->rangeLow << "\t" << targetVertex->rangeHigh << "\t" << targetVertex->major_alter_paths_num << "\t" << targetVertex->ASMcategory;
	for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; ++sampleLoopCnt)
	{
		outfile_stat_asm << "\t" << targetVertex->support[sampleLoopCnt] << "\t" << targetVertex->MSE_estimation[sampleLoopCnt];
	}
	outfile_stat_asm << endl;

	//output alternative paths
	curAlterPath = targetVertex->major_alter_paths;
	while (curAlterPath != NULL)
	{
		outfile_stat_asm << curAlterPath->whole_path_start << "\t" << curAlterPath->whole_path_end;
		for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; ++sampleLoopCnt)
		{
			outfile_stat_asm << "\t" << curAlterPath->support[sampleLoopCnt] << "\t" << curAlterPath->proportion[sampleLoopCnt];
		}
		outfile_stat_asm << endl;
		
		curAlterPath = curAlterPath->next;
	}

//	output_asm_path(targetVertex, &ASMpath_gtf_file, ASM_id);

	return;
}


void output_ASM_composition(GTvertex *targetVertex)
{
	//output the composition of an ASM, including its children and extending all its type-1 children
	double totalExpr = 0.0;
	long composCnt = 0;
	rangeJunction *curJunc;
	GTvertex *curVertex, *outputVertex;
	GTedge *curEdge, *outputEdge;

	//output basic information
	outfile_asm_composition << resultNamePrefix << "\t" << targetVertex->rangeLow << "\t" << targetVertex->rangeHigh << "\t" << targetVertex->ASMcategory << "\t";

	for (int sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; ++sampleLoopCnt)
		totalExpr += targetVertex->support[sampleLoopCnt];
	outfile_asm_composition << totalExpr / SUPPORT_VECTOR_SIZE << "\t";

	//count number of compositions & output compositions
	if (targetVertex->childType == 3)
	{
		curJunc = targetVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
				++composCnt;
			curJunc = curJunc->next;
		}

		outfile_asm_composition << composCnt << "\t";

		curJunc = targetVertex->junctionInRange->list;
		while (curJunc != NULL)
		{
			if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
			{
				if (curJunc->junc->type == frag_junction)
					outfile_asm_composition << "1\t";
				else
					outfile_asm_composition << "0\t";
				if (curJunc->junc->start_real > 0)
					outfile_asm_composition << curJunc->junc->start_real << "\t";
				else 
					outfile_asm_composition << curJunc->junc->start << "\t";
				if (curJunc->junc->end_real > 0)
					outfile_asm_composition << curJunc->junc->end_real << "\t";
				else 
					outfile_asm_composition << curJunc->junc->end << "\t";
			}
			curJunc = curJunc->next;
		}
	}
	else if (targetVertex->childType == 2)
	{
		curEdge = targetVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;
			if (curVertex->estimated == false)
			{
				curJunc = curVertex->junctionInRange->list;
				while (curJunc != NULL)
				{
					if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
						++composCnt;
					curJunc = curJunc->next;
				}
			} 
			else
			{
				if (curVertex->childType == 1)
				{
					composCnt += curVertex->childNum;
				}
				else
					++composCnt;
			}

			curEdge = curEdge->next;
		}

		outfile_asm_composition << composCnt << "\t";

		curEdge = targetVertex->child;
		while (curEdge != NULL)
		{
			curVertex = curEdge->linkedVertex;
			if (curVertex->estimated == false)
			{
				curJunc = curVertex->junctionInRange->list;
				while (curJunc != NULL)
				{
					if (curJunc->junc->type == frag_exon || curJunc->junc->type == frag_retained_intron || curJunc->junc->type == frag_junction)
					{
						if (curJunc->junc->type == frag_junction)
							outfile_asm_composition << "1\t";
						else
							outfile_asm_composition << "0\t";
						if (curJunc->junc->start_real > 0)
							outfile_asm_composition << curJunc->junc->start_real << "\t";
						else 
							outfile_asm_composition << curJunc->junc->start << "\t";
						if (curJunc->junc->end_real > 0)
							outfile_asm_composition << curJunc->junc->end_real << "\t";
						else 
							outfile_asm_composition << curJunc->junc->end << "\t";
					}
					curJunc = curJunc->next;
				}
			} 
			else
			{
				if (curVertex->childType == 1)
				{
					outputEdge = curVertex->child;
					while (outputEdge != NULL)
					{
						outputVertex = outputEdge->linkedVertex;
						if (outputVertex->childType == 0)
						{
							if (outputVertex->junctionInRange->list->junc->type == frag_junction)
								outfile_asm_composition << "1\t";
							else
								outfile_asm_composition << "0\t";
							if (outputVertex->junctionInRange->list->junc->start_real > 0)
								outfile_asm_composition << outputVertex->junctionInRange->list->junc->start_real << "\t";
							else
								outfile_asm_composition << outputVertex->junctionInRange->list->junc->start << "\t";
							if (outputVertex->junctionInRange->list->junc->end_real > 0)
								outfile_asm_composition << outputVertex->junctionInRange->list->junc->end_real << "\t";
							else
								outfile_asm_composition << outputVertex->junctionInRange->list->junc->end << "\t";
						}
						else
						{
							outfile_asm_composition << "2\t" << outputVertex->rangeLow << "\t" << outputVertex->rangeHigh << "\t";
						}
						
						outputEdge = outputEdge->next;
					}
				}
				else
					outfile_asm_composition << "2\t" << curVertex->rangeLow << "\t" << curVertex->rangeHigh << "\t";
			}

			curEdge = curEdge->next;
		}
	}

	outfile_asm_composition << endl;

	return;
}

void output_asm_analysis(GTvertex *rootVertex)
{
	//output all asms for differential transcription analysis
	GTvertex *curVertex;
	GTedge *curEdge;
	
	asm_category(rootVertex);

	if (rootVertex->level < 1)
	{
		++GENEcount;
	}

	if (rootVertex->level < 1 && rootVertex != gTree->root)
	{
		long segmentCnt = 1;
		
		if (rootVertex->child == NULL && rootVertex->exonNum > 0)
		{
			//rootVertex is a exon
			outfile_gtf_splice_graph << resultNamePrefix << "\t" << "ASM" << "\t" << "exon" << "\t" << rootVertex->rangeLow << "\t" << rootVertex->rangeHigh << "\t" << "." << "\t" << "+\t" << "." << "\t"
				<< "gene_id \"" << resultNamePrefix << "_" << GENEcount << "\"; " << "transcript_id \"" << resultNamePrefix << "_" << GENEcount << "_" << segmentCnt << "\";" << endl;
		}
		else
		{
			curEdge = rootVertex->child;
			while (curEdge != NULL)
			{
				if (curEdge->linkedVertex->child == NULL)
				{
					if (curEdge->linkedVertex->exonNum > 0)
						outfile_gtf_splice_graph << resultNamePrefix << "\t" << "ASM" << "\t" << "exon" << "\t" << curEdge->linkedVertex->rangeLow << "\t" << curEdge->linkedVertex->rangeHigh << "\t" << "." << "\t" << "+\t" << "." << "\t"
						<< "gene_id \"" << resultNamePrefix << "_" << GENEcount << "\"; " << "transcript_id \"" << resultNamePrefix << "_" << GENEcount << "_" << segmentCnt << "\";" << endl;
				} 
				else
				{
					++segmentCnt;
				}
				curEdge = curEdge->next;
			}
		}
	}

	if (rootVertex->child == NULL)
	{
		//alternative splice site
		if (rootVertex->alterSpliceSite != NULL)
		{
			curEdge = rootVertex->alterSpliceSite;
			while (curEdge != NULL)
			{
				curVertex = curEdge->linkedVertex;
				if (curVertex->ASMcategory > 0 && curVertex->major_alter_paths_num > 1)
				{
					output_one_asm(curVertex);
					output_ASM_composition(curVertex);
				}

				curEdge = curEdge->next;
			}
		}

		return;
	} 
	else
	{
		//extend children
		curEdge = rootVertex->child;
		while (curEdge != NULL)
		{
			output_asm_analysis(curEdge->linkedVertex);
			curEdge = curEdge->next;
		}

		if ((rootVertex->childType == 2 || rootVertex->childType == 3) && (rootVertex->major_alter_paths_num > 1))
		{
			output_one_asm(rootVertex);
			output_ASM_composition(rootVertex);
		}

		if (rootVertex->alterSpliceSite != NULL)
		{
			curEdge = rootVertex->alterSpliceSite;
			while (curEdge != NULL)
			{
				curVertex = curEdge->linkedVertex;
				if (curVertex->ASMcategory > 0 && curVertex->major_alter_paths_num > 1)
				{
					output_one_asm(curVertex);
					output_ASM_composition(curVertex);
				}
				
				curEdge = curEdge->next;
			}
		}

		return;
	}
}

void differential_analysis_asm()
{
	output_asm_analysis(gTree->root);

	return;
}


/************************************************************************/
/* MAIN                                                                 */
/************************************************************************/

void initialization()
{
	string filename;

#ifdef UNIX
	filename = resultPath + "detail/" + resultNamePrefix + "_GTree.txt";
	outfile_gtree.open(filename.c_str());
	filename = resultPath + "detail/" + resultNamePrefix + "_not_enumerated.txt";
	outfile_not_enumerated.open(filename.c_str());
	
	filename = resultPath + "stat/expression.txt";
	outfile_stats_expr.open(filename.c_str(), fstream::app);
	outfile_stats_expr.precision(10);
	filename = resultPath + "stat/asm.txt";
	outfile_stat_asm.open(filename.c_str(), fstream::app);
	outfile_stat_asm.precision(10);

	filename = resultPath + "asm_path.gtf";
	outfile_gtf_asm_path.open(filename.c_str(), fstream::app);
	filename = resultPath + "splice_graph.gtf";
	outfile_gtf_splice_graph.open(filename.c_str(), fstream::app);

	filename = resultPath + "stat/splice_all.bed";
	outfile_junction_all.open(filename.c_str(), fstream::app);
	filename = resultPath + "stat/splice_filtered.bed";
	outfile_junction_filtered.open(filename.c_str(), fstream::app);

	filename = resultPath + "detail/" + resultNamePrefix + "_asm.txt";
	outfile_asm_composition.open(filename.c_str());
#else
	filename = resultPath + "GTree.txt";
	outfile_gtree.open(filename.c_str());
	filename = resultPath + "unresolvable.txt";
	outfile_not_enumerated.open(filename.c_str());
	
	filename = resultPath + "stat\\expression.txt";
	outfile_stats_expr.open(filename.c_str());
	outfile_stats_expr.precision(10);
	filename = resultPath + "stat\\asm.txt";
	outfile_stat_asm.open(filename.c_str());
	outfile_stat_asm.precision(10);

	filename = resultPath + "asm_path.gtf";
	outfile_gtf_asm_path.open(filename.c_str());
	outfile_gtf_asm_path << "browser full ASM\ntrack name=\"ASM\" description=\"ASM\" visibility=2 useScore=1\n";
	filename = resultPath + "non_asm.gtf";
	outfile_gtf_splice_graph.open(filename.c_str());
	outfile_gtf_splice_graph << "browser dense Gene\ntrack name=\"Gene structure\" description=\"Gene structure\" visibility=1 useScore=1\n";
#endif
	
	//read statFile
	ifstream statfile;
#ifdef UNIX
	filename = inputPath + "stat.txt";
	statfile.open(filename.c_str());
#else
	statfile.open("tmp\\stat.txt");
#endif

	long tmp = -1, tmp_small = MAX_CHR_LENGTH, tmp_large = -1;
	while (statfile >> tmp)
	{
		if (tmp < tmp_small)
			tmp_small = tmp;
		
		statfile >> tmp;
		if (tmp > tmp_large)
			tmp_large = tmp;
	}
	CHROMOSOME_START = tmp_small;
	CHROMOSOME_END = tmp_large;
	statfile.close();

	srand ( time(NULL) );

	
	//read dataset stat file
	DatasetReadCount_total.assign(SUPPORT_VECTOR_SIZE, 0);
	DatasetNormalizationRatio.assign(SUPPORT_VECTOR_SIZE, 1.0);

	ifstream datasetStatFile;
	filename = inputPath + "DatasetStat.txt";
	datasetStatFile.open(filename.c_str());
	if (datasetStatFile.is_open() == true)
	{
		for (int tmpCnt = 0; tmpCnt < SUPPORT_VECTOR_SIZE; tmpCnt++)
		{
			datasetStatFile >> DatasetReadCount_total[tmpCnt];
			if (DatasetReadCount_total[tmpCnt] > DatasetReadCount_normalization)
				DatasetReadCount_normalization = DatasetReadCount_total[tmpCnt];
		}
		datasetStatFile.close();
	}

	sortList_Junction.resize(DEFAULT_JUNCTION_NUM*3+1);
	sortKey_Junction.assign(DEFAULT_JUNCTION_NUM*3+1, 0);
	mergeSort_Larray.assign(DEFAULT_JUNCTION_NUM*3+1, 0);
	mergeSort_Rarray.assign(DEFAULT_JUNCTION_NUM*3+1, 0);
	mergeSort_LorderedList.resize(DEFAULT_JUNCTION_NUM*3+1);
	mergeSort_RorderedList.resize(DEFAULT_JUNCTION_NUM*3+1);
	

	return;
}



void cleanAll()
{
	outfile_gtree.close();
	outfile_not_enumerated.close();
	outfile_stats_expr.close();
	outfile_stat_asm.close();
	outfile_gtf_asm_path.close();
	outfile_gtf_splice_graph.close();
	outfile_junction_all.close();
	outfile_junction_filtered.close();
	outfile_asm_composition.close();

	while (juncPathQueue.empty() == false)
	{
		juncPathQueue.pop();
	}

	vertexListForStatistics.clear(); free_vector(vertexListForStatistics);
	DatasetReadCount_total.clear(); free_vector(DatasetReadCount_total);
	DatasetNormalizationRatio.clear(); free_vector(DatasetNormalizationRatio);

	sortList_Junction.clear(); free_vector(sortList_Junction);
	sortKey_Junction.clear(); free_vector(sortKey_Junction);
	mergeSort_Larray.clear(); free_vector(mergeSort_Larray);
	mergeSort_Rarray.clear(); free_vector(mergeSort_Rarray);
	mergeSort_LorderedList.clear(); free_vector(mergeSort_LorderedList);
	mergeSort_RorderedList.clear(); free_vector(mergeSort_RorderedList);

	return;
}


void inputData()
{
	string filename;
	int sampleLoopCnt;
	long prevJuncNum;
	fragment *curFrag;

	/************************************************************************/
	/* INPUT ALL JUNCTIONS                                                  */
	/************************************************************************/
#ifdef UNIX
	filename = inputPath + "junction_all.txt";
#else
	filename = "tmp\\junction_all.txt";
#endif
	input_junction(filename);

#ifdef COUTSTEPS
	cout << "input junction done." << endl;
#endif


	/************************************************************************/
	/* get junction support                                                 */
	/************************************************************************/

#ifndef JUNCTIONONLY
	for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; sampleLoopCnt++)
	{
		getJunctionSupport(sampleLoopCnt);
	}
	cleanSpliceDirection();
#endif

#ifdef COUTSTEPS
	cout << "junction support done." << endl;
#endif
	cout << "total junction = " << sortList_Junction_Num << "... " << flush;

	/************************************************************************/
	/* junction filtering                                                   */
	/************************************************************************/

	for (long i = 1; i <= sortList_Junction_Num; i++)
		if (sortList_Junction[i]->junc->type == frag_junction)
		{
			curFrag = sortList_Junction[i]->junc;
			outfile_junction_all << curFrag->chromosome_start << "\t" << curFrag->start - 19 << "\t" << curFrag->end + 20 << "\talljunction\t1000\t";
			if (curFrag->transDirection == undetermined || curFrag->transDirection == sense)
				outfile_junction_all << "+\t";
			else
				outfile_junction_all << "-\t";
			outfile_junction_all << curFrag->start - 19 << "\t" << curFrag->end + 20 << "\t0\t2\t" << "20,20,\t0," << curFrag->end - curFrag->start + 19 << endl;
		}
		
	
#ifdef FILTER_JUNCTION
	filterJunction(thresh_junction_filter_max_read_support, SUPPORT_VECTOR_SIZE - thresh_junction_filter_num_samples_presence, thresh_junction_filter_mean_read_support); //throw a junction away if zeroCnt >= 2nd parameter or max support <= 1st parameter or mean support <= 3rd parameter
	cout << "after filtering = " << sortList_Junction_Num << "... " << flush;
#endif

	for (long i = 1; i <= sortList_Junction_Num; i++)
		if (sortList_Junction[i]->junc->type == frag_junction)
		{
			curFrag = sortList_Junction[i]->junc;
			outfile_junction_filtered << curFrag->chromosome_start << "\t" << curFrag->start - 19 << "\t" << curFrag->end + 20 << "\tfilteredjunction\t1000\t";
			if (curFrag->transDirection == undetermined || curFrag->transDirection == sense)
				outfile_junction_filtered << "+\t";
			else
				outfile_junction_filtered << "-\t";
			outfile_junction_filtered << curFrag->start - 19 << "\t" << curFrag->end + 20 << "\t0\t2\t" << "20,20,\t0," << curFrag->end - curFrag->start + 19 << endl;
		}

	prevJuncNum = sortList_Junction_Num;

	getExons();
	//cout << "get exons done." << endl;
#ifndef JUNCTIONONLY
	getExonSupport(-1, prevJuncNum + 1);
	DOTRIMMING = false;
	DORETAIN = false;
	for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; sampleLoopCnt++)
	{
		getExonSupport(sampleLoopCnt, prevJuncNum + 1);
	}
#endif

	/************************************************************************/
	/* CALCULATE COVERAGE                                                   */
	/************************************************************************/
#ifdef COVERAGE
	//calculate read coverage
	for (long i = 1; i <= sortList_Junction_Num; i++)
	{
		if (sortList_Junction[i]->junc->type == frag_exon || sortList_Junction[i]->junc->type == frag_retained_intron)
		{
			for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; sampleLoopCnt++)
			{
				sortList_Junction[i]->junc->support[sampleLoopCnt] = sortList_Junction[i]->junc->support[sampleLoopCnt] / (sortList_Junction[i]->junc->end - sortList_Junction[i]->junc->start + 1);
			}
		}
	}
#endif

#ifdef NORMALIZE_COVERAGE
	//normalize exon coverage and junction support
	for (sampleLoopCnt = 0; sampleLoopCnt < SUPPORT_VECTOR_SIZE; sampleLoopCnt++)
	{
		if (DatasetReadCount_total[sampleLoopCnt] > 0)
		{
			DatasetNormalizationRatio[sampleLoopCnt] = DatasetReadCount_normalization / DatasetReadCount_total[sampleLoopCnt];
			for (long i = 1; i <= sortList_Junction_Num; i++)
			{
				sortList_Junction[i]->junc->support[sampleLoopCnt] = sortList_Junction[i]->junc->support[sampleLoopCnt] * DatasetNormalizationRatio[sampleLoopCnt];
			}
		}
	}
#endif

	//cout << "total " << sortList_Junction_Num << "\t";
#ifdef FILTER_FRAGMENTS
	filterIntron(prevJuncNum + 1); //filter introns
	//cout << "filtering intron... " << sortList_Junction_Num << "\t";
#endif

#ifdef COUTSTEPS
	cout << "get junctions and exons done." << endl;
#endif

	return;
}


void input_config(string filename)
{
	ifstream config_file;
	config_file.open(filename.c_str());

	string parameter, info;

	if (config_file.is_open() == true)
	{
		while (config_file >> parameter)
		{
			if (parameter.compare("SUPPORT_VECTOR_SIZE") == 0)
				config_file >> SUPPORT_VECTOR_SIZE;
			else if (parameter.compare("thresh_junction_filter_max_read_support") == 0)
				config_file >> thresh_junction_filter_max_read_support;
			else if (parameter.compare("thresh_junction_filter_mean_read_support") == 0)
				config_file >> thresh_junction_filter_mean_read_support;			
			else if (parameter.compare("thresh_junction_filter_num_samples_presence") == 0)
				config_file >> thresh_junction_filter_num_samples_presence;
			else if (parameter.compare("COUNT_MAJOR_PATHS_ONLY") == 0)
			{
				config_file >> info;
				if (info.compare("yes") == 0)
					COUNT_MAJOR_PATHS_ONLY = true;
				else
					COUNT_MAJOR_PATHS_ONLY = false;
			}
			else if (parameter.compare("coverageThreshold_exon") == 0)
				config_file >> coverageThreshold_exon;
			else if (parameter.compare("coverageThreshold_intron") == 0)
				config_file >> coverageThreshold_intron;

			getline(config_file, info);
		}
	} 
	else
	{
		cout << "Error: fail to open the config file for gtree" << endl;
		exit(1);
	}

	config_file.close();

	if (SUPPORT_VECTOR_SIZE <= 0)
	{
		cout << "Error: incomplete config file for gtree" << endl;
		exit(1);
	}

	return;
}



int main(int argc, char *argv[])
{
	RangeJunctionList* origList;
	string config_filename;

#ifdef COUTSTEPS
	cout << "GTree begins..." << endl;
#endif

#ifdef UNIX
	if (argc == 5)
	{
		inputPath = argv[1];
		resultPath = argv[2];
		resultNamePrefix = argv[3];
		config_filename = argv[4];

		cout << resultNamePrefix << "\t";
	}	
	else
	{
		cout << argv[0] << "\t<inputPath>\t<resultPath>\t<resultNamePrefix>\t<config_file>\n";
		exit(1);
	}
#else
	if (argc == 1)
	{
		inputPath = "";
		resultPath = "result\\";
		resultNamePrefix = "chr7";
		config_filename = "config_gtree";

		cout << resultNamePrefix << "\t";
	}
	else
	{
		exit(1);
	}
#endif	

	input_config(config_filename);
	initialization();	

	/************************************************************************/
	/* INPUT COVERAGE                                                       */
	/************************************************************************/
	inputData();


	/************************************************************************/
	/* CONSTRUCTION                                                         */
	/************************************************************************/

	//build tree
	gTree = new GenomeTree;
	origList = buildOrigList();
	cout << "total exon and junction = " << sortList_Junction_Num << "... " << flush;
	
#ifdef COUTSTEPS
	cout << "Constructing... " << endl;
#endif

	cout << "Decomposing the splice graph... " << flush;
	constructGTree(origList);
	
#ifdef COUTSTEPS
	cout << "Counting... " << endl;
#endif




	/************************************************************************/
	/* ESTIMATION                                                           */
	/************************************************************************/
		
	//find alternative start and end
//	preCountGTree(gTree->root);
#ifdef COUTSTEPS
	cout << "pre-count done" << endl;
#endif
//	alterStart(gTree->root);
	
#ifdef COUTSTEPS
	cout << " Alter start done" << endl;
#endif

	//recount
	cout << "Estimating the transcript abundance... " << flush;
	countGTree(gTree->root);

#ifdef COUTSTEPS
	cout << "Counting done." << endl;
#endif

	cout << "Writing output files... " << flush;
	/************************************************************************/
	/* DIFFERENTIAL ANALYSIS                                                */
	/************************************************************************/
	differential_analysis_expression_level();
	
	differential_analysis_asm();

	/************************************************************************/
	/* OUTPUT                                                               */
	/************************************************************************/

	//output
	GTree_output(gTree->root, 1, &outfile_gtree);
	
	cout << "Finished" << flush;
	cleanAll();

//	cin >> filename;
	return 0;
}


