// *******************************************************************
//      file with specific functions to solve the HMP
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Global Variables
extern int n;                               // size of cromossoms

//---------------------- DEFINITION OF TYPES OF PROBLEM SPECIFIC --------------------


//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  --------

static std::vector <std::vector <int> > hij;	// matrix with edge weight
static std::vector <double> ti;				    // vector of nodes weigth
static int numNodes;							// number of nodes
static int numClusters;							// number of clusters
static double cap;								// cluster capacity
static int totalHij;                            // total edge weigth


//----------------------- IMPLEMENTATION OF FUNCTIONS  -------------------------------


/************************************************************************************
 Method: ReadData
 Description: read input data of the problem
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
    fscanf(arq, "%d %d %lf", &numNodes, &numClusters, &cap);

    // read graph informations
    ti.clear();
    ti.resize(numNodes);

    hij.clear();
    hij.resize(numNodes, std::vector<int>(numNodes));

    for (int i=0; i<numNodes; i++)
    {
    	fscanf(arq, "%lf", &ti[i]);
    }

    totalHij = 0;
    for (int i=0; i<numNodes; i++)
    {
        for (int j=0; j<numNodes; j++)
        {
    	    fscanf(arq, "%d", &hij[i][j]);
            totalHij += hij[i][j];
        }
    }

    fclose(arq);
    
    n = numNodes+1; // a ultima posicao (n-1) representa o numero de base que adicionamos inicialmente em rncs diferentes


    // imprimir
    // printf("\n\n%d %d %lf \n", numNodes, numClusters, cap);

    // for (int i=0; i<numNodes; i++)
    // {
    // 	printf("%lf\n", ti[i]);
    // }

    // for (int i=0; i<numNodes; i++)
    // {   
    //     printf("\n");
    //     for (int j=0; j<numNodes; j++)
    //     {
    // 	    printf("%d\t", hij[i][j]);
    //     }
    // }
    // printf("\n\n");
    // getchar();
}

/************************************************************************************
 Method: CalculateFitness
 Description: calculate the fitness of a chromossom s
*************************************************************************************/
int CalculateFitness(std::vector< std::vector<int> > RNC)
{
    // calculate objective function
    int soma = 0;

    // somar os pesos de cada RNC
    // printf("\nRNC = %d; \t", RNC.size());
    for (int k=0; k<(int)RNC.size(); k++)
    {
        // printf("%d \t", RNC[k].size());
        // total handover between base stations assigned to the same RNC
        for (int i=0; i<(int)RNC[k].size(); i++)
        {
            for (int j=0; j<(int)RNC[k].size(); j++)
            {
                soma += hij[RNC[k][i]][RNC[k][j]];
            }
        }
    }
    return totalHij - soma;

    // for (int k=0; k<(int)RNC.size(); k++)
    // {
    //     for (int i=0; i<(int)RNC[k].size(); i++)
    //     {
    //         for (int k1=0; k1<(int)RNC.size(); k1++)
    //         {
    //             if (k1 != k)
    //             {
    //                 for (int j=0; j<(int)RNC[k1].size(); j++)
    //                 {
    //                     soma += hij[RNC[k][i]][RNC[k1][j]];
    //                 }
    //             }
    //         }
    //     }
    //     // soma += ((double)(rand()%10000)/10000.0)*(100-10)+10;
    // }

    // return soma;
}

/************************************************************************************
 Method: Decoders
 Description: Convert a random key solution in a real problem solution
*************************************************************************************/
double Decoder(TSol s)
{
    // best RNC with numRNC initial alocation

    // definir o numero de RNC usados inicialmente
    int numRNC = 0;

    // create a initial solution of the problem
    std::vector <int> sol(n); 
    for (int j = 0; j < n; j++)
	{
        if (j < n-1)
		    sol[j] = j;
        else
        if (j == n-1)
            numRNC = (s.rk[n-1] * (numClusters)) + 1; // define o numero de bases alocadas em diferentes RNCs
            // s.vec[j].sol = (s.vec[n-1].rk * (numClusters/2)) + numClusters/2; 
            // s.vec[j].sol = numClusters;
	}

    // sort the problem vector based on the values in the rk vector
    std::sort(sol.begin(), sol.end()-1, [&s](int i1, int i2) {
        return s.rk[i1] < s.rk[i2];
    });

    // criar as estruturas dos RNCs
    std::vector< std::vector<int> > RNC;
    std::vector<double> capRNC;

    RNC.resize(numClusters);
    capRNC.resize(numClusters);

    //empacotar as bases nos rncs
    int penalty = 0;
    for(int i=0; i<numNodes; i++)
    {
        int bestRNC = -1;
        int bestCost = -1;

        // alocar os primeiros numRNC bases em RNCs separados
        if(i<numRNC)
        {
            capRNC[i] += ti[sol[i]];
            RNC[i].push_back(sol[i]);
        }
        else
        {
            //encontrar o melhor rnc k que possui capacidade para atender o node i
            int k=0;
            while(k < numClusters)
            {
                //verificar a capacidade disponivel do rnc k
                if (capRNC[k] + ti[sol[i]] <= cap)
                {
                    // calcular o custo de inserir o ponto i no rnc k
                    int cost = 0;
                    for (int j=0; j<(int)RNC[k].size(); j++)
                    {
                        cost += hij[RNC[k][j]][sol[i]] + hij[sol[i]][RNC[k][j]];
                    }

                    if (cost > bestCost)
                    {
                        bestRNC = k;
                        bestCost = cost;
                    }
                }
                k++;
            }

            //inserir o node i no melhor rnc aberto disponivel
            if (bestRNC >= 0)
            {
                // atualiza a capacidade 
                capRNC[bestRNC] += ti[sol[i]];

                // insere a base 
                RNC[bestRNC].push_back(sol[i]);
            }

            // penalizar a solucao se um node nao foi agrupado
            else
            {
                penalty += 100000 + ti[sol[i]] * 10; 
            }
        }
    }

    s.ofv = CalculateFitness(RNC) + penalty;

    // print the solution in the screen
    if (debug && print)
    {
        for (int i=0; i<n; i++){
		    printf("%d ", sol[i]);
        }

        for (int k=0; k<(int)RNC.size(); k++)
        {
            printf("\nRNC %d (%.3lf): ", k, capRNC[k]);
            // total handover between base stations assigned to the same RNC
            for (int i=0; i<(int)RNC[k].size(); i++)
            {
                printf("%d ", RNC[k][i]);
            }
        }
    }

    // print the solution in a file
    if (!debug && print)
    {
        for (int i=0; i<n; i++){
		    fprintf(arqSol,"%d ", sol[i]);
        }

        for (int k=0; k<(int)RNC.size(); k++)
        {
            fprintf(arqSol,"\nRNC %d (%.3lf): ", k, capRNC[k]);
            // total handover between base stations assigned to the same RNC
            for (int i=0; i<(int)RNC[k].size(); i++)
            {
                fprintf(arqSol, "%d ", RNC[k][i]);
            }
        }
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
    hij.clear();
    ti.clear();
}

#endif
