// *******************************************************************
//      file with specific functions to solve the alpha-NpMP
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Global Variables
extern int n;                               // size of random-key vector

//------ DEFINITION OF TYPES OF PROBLEM ON HAND --------


//------ DEFINITION OF GLOBAL VARIABLES OF PROBLEM ON HAND --------

static std::vector <std::vector <int> > dist;	// matrix with Euclidean distance

static int numNodes;						    // number of nodes
static int numMedians;							// number of medians
static int alphaN;								// number of alpha neighbors


//-------------------------- FUNCTIONS OF PROBLEM OON HAND --------------------------

/************************************************************************************
 Method: ReadData
 Description: read the input data
*************************************************************************************/
void ReadData(char nameTable[])
{ 
    char name[200] = "../Instances/";
    strcat(name,nameTable);

    FILE *arq;
    arq = fopen(name,"r");

    if (arq == NULL)
    {
        printf("\nERROR: File (%s) not found!\n",name);
        getchar();
        exit(1);
    }

    // => read data
    int aux;
    if (fscanf(arq, "%d %d %d", &numNodes, &aux, &numMedians) != 3) {
        printf("Erro ao ler os dados do arquivo do problema.\n");
    }

    numMedians = 10;
    alphaN = 5;

    // numMedians = 20;
    // alphaN = 10;

    dist.clear();
    dist.resize(numNodes, std::vector<int>(numNodes, 9999999));
    for (int i = 0; i < numNodes; i++){
        dist[i][i] = 0;
    }

    // read distance informations
    while (!feof(arq)) //for (int k = 0; k < aux; k++)
    {
        int i, j, distance;
        if (fscanf(arq, "%d %d %d", &i, &j, &distance) != 3) {
            printf("Erro ao ler os dados do arquivo problema - p2.\n");
        }
    	// fscanf(arq, "%d %d %d", &i, &j, &distance);
    	dist[i-1][j-1] = distance;
        dist[j-1][i-1] = distance;

        // printf("\n%d \t %d \t %d", i, j, dist[i-1][j-1]);
    }
    fclose(arq);

    // Floyd-Warshall algorithm to compute the shortest paths between every pair of vertices
    /* Add all vertices one by one to the set of intermediate vertices.
      ---> Before start of a iteration, we have shortest distances between all
      pairs of vertices such that the shortest distances consider only the
      vertices in set {0, 1, 2, .. k-1} as intermediate vertices.
      ----> After the end of a iteration, vertex no. k is added to the set of
      intermediate vertices and the set becomes {0, 1, 2, .. k} */

    for (int k = 0; k < numNodes; k++){
        // Pick all vertices as source one by one
        for (int i = 0; i < numNodes; i++){
            // Pick all vertices as destination for the above picked source
            for (int j = 0; j < numNodes; j++){
                // If vertex k is on the shortest path from i to j, then update the value of dist[i][j]
                if (dist[i][k] + dist[k][j] < dist[i][j])
                    dist[i][j] = dist[i][k] + dist[k][j];
            }
        }
    }

    // print
    // for (int i = 0; i < numNodes; i++){
    //     printf("\n");
    //     for (int j = 0; j < numNodes; j++){
    //         printf("%d ", dist[i][j]);
    //     }
    // }
    // getchar();

    

    // define the size of the solution vector
    n = numMedians;
}

/************************************************************************************
 Method: CalculateFitness
 Description: calculate the fitness of a chromossom s
*************************************************************************************/
double CalculateFitness(std::vector <int> sol)
{    
    double ofv = 0;
    // std::vector <std::vector <int> > solution;
    // solution.resize(numNodes, std::vector<int>(alphaN));

    // for each node find the alpha nearest medians 
    for (int j=0; j<numNodes; j++)
    {
        std::vector<int> medians(numMedians,0);
        int minDist;
        int minIndex;

        // find the k nearest median
        for (int k=0; k<alphaN; k++)
        {
            minDist = 99999999;
            minIndex = 0;

            for (int i=0; i<numMedians; i++)
            {
                if ( (dist[sol[i]][j] < minDist) && (medians[i] == 0) )
                {
                    minDist = dist[sol[i]][j];
                    minIndex = i;
                }
            }

            // insert the alpha nearest median in the solution
            // solution[j][k] = s.vec[minIndex].sol;
            medians[minIndex] = 1;

            ofv += minDist;

            // printf("\n");
            // for (int i=0; i<numMedians; i++){
            //     printf("%d \t", medians[i]);
            // }
        }
        // getchar();
    }

    // printf("\n");
    // for (int i=0; i<numMedians; i++)
    // {
    //     printf("%d \t", s.vec[i].sol);
    // }
    // printf("\n");

    // for (int j=0; j<numNodes; j++)
    // {
    //     printf("\n");
    //     for (int k=0; k<alphaN; k++)
    //     {
    //         printf("%d \t", solution[j][k]);
    //     }
    // }
    // getchar();


    return ofv;
}

/************************************************************************************
 Method: Decoders
 Description: users need to implement at least one decoder, DecK (K = [1,2,3,4,5])
*************************************************************************************/
double Decoder(TSol s) 
{
    // create an initial solution of the problem
    std::vector<int> nodesList;

    // vector with p medians
    std::vector <int> sol(n);

    // list of candidate nodes
    for (int j = 0; j < numNodes; j++){
        nodesList.push_back(j);
	}

    // define the p medians
    for (int i=0; i<n; i++)
    {
        // define the index of the nodeList
        int index = 0;
        if (s.rk[i] >= 0 || s.rk[i] < 1.0)
            index = floor(s.rk[i] * (int)nodesList.size());

        // define the node that is a median
        sol[i] = nodesList[index];

        // remove this node of the nodeList
        nodesList.erase(nodesList.begin()+index);
    }

    s.ofv = CalculateFitness(sol);

    // print the solution in the screen
    if (debug && print)
    {
        for (int i=0; i<n; i++)
		    printf("%d ", sol[i]);
    }

    // print the solution in a file
    if (!debug && print)
    {
        for (int i=0; i<n; i++)
		    fprintf(arqSol,"%d ", sol[i]);
    }

    return s.ofv;
}


/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/

void FreeMemoryProblem()
{
    //specific problem
    dist.clear();
}

#endif