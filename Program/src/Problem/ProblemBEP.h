// ************************************************************************************
//      file with specific functions to solve the Balanced Edge Partition Problem
// ************************************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Variables declared in main.cpp
extern int n;                               // size of cromossoms

//----------------- DEFINITION OF PROBLEM SPECIFIC TYPES___ -----------------------
struct Edge {
    int node1;
    int node2;
};

struct GraphData {
    int numNodes;
    int numEdges;
    std::vector<Edge> edges;
};
GraphData graphData;

static std::vector <std::vector <int> > nodesEdges;	// matrix with index of edges of echa node

static int K = 0;    // numero de particoes

//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  --------



//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------

/************************************************************************************
 Method: ReadData
 Description: read the input data
*************************************************************************************/
void ReadData(char nameTable[])
{ 
     char name[200] = "../Instances/";
    strcat(name,nameTable);

    // std::ifstream para leitura de arquivos em C++
    std::ifstream file(name);
    std::string line;
    bool inGraphSection = false;

    // => read data
    if (file.is_open()) {
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string token;
            ss >> token;

            if (token == "SECTION" && ss >> token && token == "Graph") {
                inGraphSection = true;
            } else if (inGraphSection) {
                if (token == "Nodes") {
                    ss >> graphData.numNodes;
                } else if (token == "Edges") {
                    ss >> graphData.numEdges;
                } else if (token == "E") {
                    Edge edge;
                    ss >> edge.node1 >> edge.node2;
                    graphData.edges.push_back(edge);
                } else if (token == "END") {
                    inGraphSection = false;
                }
            }
        }
        file.close();
    } else {
        std::cerr << "Erro ao abrir o arquivo: " << name << std::endl;
    }

    nodesEdges.clear();
    nodesEdges.resize(graphData.numNodes);
    for (int i=0; i<graphData.numEdges; i++)
    {
        // corrigir o indice
        graphData.edges[i].node1 = graphData.edges[i].node1-1;
        graphData.edges[i].node2 = graphData.edges[i].node2-1;

        nodesEdges[graphData.edges[i].node1].push_back(i);
        nodesEdges[graphData.edges[i].node2].push_back(i);
    }

    // print
    // printf("\nEdges: %d \nNodes: %d", graphData.numEdges, graphData.numNodes);
    // for (int i=0; i<graphData.edges.size(); i++)
    // {
    //     printf("\n%d %d", graphData.edges[i].node1, graphData.edges[i].node2);
    // }
    // for (int i=0; i<graphData.numNodes; i++)
    // {
    //     printf("\nNode %d: ",i+1);
    //     for (unsigned int j=0; j<nodesEdges[i].size(); j++)
    //         printf("%d ", nodesEdges[i][j]);
    // }
    // getchar();

    // define the size of the solution vector
    // n = graphData.numEdges + graphData.numNodes + 1;
    // n = graphData.numEdges + graphData.numNodes; 
    // n = graphData.numEdges + 1;
    n = graphData.numEdges; 
}

double CalculateOFV(std::vector <int> X, std::vector <int> Y)
{
    int cost = 0;
    double total = 0;
    for (int i=0; i<graphData.numEdges; i++)
    {
        if (Y[graphData.edges[i].node1] == Y[graphData.edges[i].node2] && Y[graphData.edges[i].node1] == X[i]){
            cost = 0;
        }

        else 
        if (Y[graphData.edges[i].node1] != Y[graphData.edges[i].node2] && Y[graphData.edges[i].node1] == X[i]){
            cost = 1;
        }

        else 
        if (Y[graphData.edges[i].node1] != Y[graphData.edges[i].node2] && Y[graphData.edges[i].node2] == X[i]){
            cost = 1;
        }

        else 
        if (Y[graphData.edges[i].node1] == Y[graphData.edges[i].node2] && Y[graphData.edges[i].node1] != X[i]){
            cost = 2;
        }

        else 
        if (Y[graphData.edges[i].node1] != Y[graphData.edges[i].node2] && Y[graphData.edges[i].node1] != X[i] && Y[graphData.edges[i].node2] != X[i]){
            cost = 3;
        }

        total += cost;

        if (debug && print && cost>0) printf("\n%d: %d %d [%d]", i, graphData.edges[i].node1, graphData.edges[i].node2, cost);
    }

    return total;
}

/************************************************************************************
 Method: Decoders
 Description: mapping the random-key solution into a problem solution
*************************************************************************************/
double Decoder(TSol s) 
{
    // create an initial list of edges
    std::vector <int> edges(graphData.numEdges);  
    for (int j = 0; j < graphData.numEdges; j++){ edges[j] = j;}

    // sort the problem vector based on the values in the rk vector
    std::sort(edges.begin(), edges.end(), [&s](int i1, int i2) {
        return s.rk[i1] < s.rk[i2];
    }); 

    // problem solution
    std::vector <int> X(graphData.numEdges,-1);     // edges allocation
    std::vector <int> Y(graphData.numNodes,-1);     // nodes allocation

    
    std::vector <int> capP(K,0);                    // capacidade das particoes
    int limite = graphData.numEdges/K;              // limite inferior da capacidade de cada particao
    
    // alocar cada edge em uma particao
    for(int i=0; i<graphData.numEdges; i++)
    {
        // alocar as primeiras K arestas em K particoes
        // if (i < K || i < graphData.numEdges * 0.3)
        // if (i < K || i < graphData.numEdges * s.rk[n-1])
        // if (i < 0) 
        // {
        //     X[edges[i]] = i%K;
        //     capP[i%K]++;
        // }

        // alocar as arestas na melhor particao que tenha capacidade
        // else
        // {
            // armazenar o numero de nos em comum em cada particao
            std::vector <int> partitions(K,0);

            // aresta
            int e = edges[i];

            // nos da aresta e
            int n1 = graphData.edges[e].node1;
            int n2 = graphData.edges[e].node2;
        
            // verificar n1
            for (unsigned int j=0; j<nodesEdges[n1].size(); j++)
            {
                // outras arestas que tem o no n1
                int eAux = nodesEdges[n1][j];

                // particao da aresta eAux ja foi alocada
                if (X[eAux] >= 0)
                {
                    partitions[X[eAux]]++;
                }
            }

            // verificar n2
            for (unsigned int j=0; j<nodesEdges[n2].size(); j++)
            {
                // outras arestas que tem o no n2
                int eAux = nodesEdges[n2][j];

                // particao da aresta eAux ja foi alocada
                if (X[eAux] >= 0)
                {
                    partitions[X[eAux]]++;
                }
            }

            // inserir na melhor particao que tenha capacidade
            int bestPartition = -1;
            int maior = -1;            
            for (int k=0; k<K; k++)
            {
                if ((partitions[k] > maior) && (capP[k] < limite))
                {
                    maior = partitions[k];
                    bestPartition = k;
                }

                // desempatar considerando a menor capacidade ocupada
                else
                if ((partitions[k] == maior) && (capP[k] < limite))
                {
                    if (capP[k] < capP[bestPartition])
                    {
                        maior = partitions[k];
                        bestPartition = k;
                    }
                }
            }

            // nao ha particao com capacidade
            if (bestPartition == -1)
            {
                // aumentar o limite e inserir na melhor particao
                limite++;

                bestPartition = -1;
                maior = -1;
                
                for (int k=0; k<K; k++)
                {
                    if ((partitions[k] > maior) && (capP[k] < limite))
                    {
                        maior = partitions[k];
                        bestPartition = k;
                    }
                }
            }

            X[edges[i]] = bestPartition;
            capP[bestPartition]++; 
        // }
    }

    // alocar os nodes
    capP.clear();
    capP.resize(K,0);
    for(int i=0; i<graphData.numNodes; i++)
    {
        // encontrar a particao mais comum entre as arestas e nos conectados ao no i
        std::vector <int> partitions(K,0);
    
        for (unsigned int j=0; j<nodesEdges[i].size(); j++)
        {
            // aresta
            int e = nodesEdges[i][j];

            // particao da aresta
            partitions[X[e]]++;

            // particao do no 1 (se ja alocado)
            if (Y[graphData.edges[e].node1] >= 0)
                partitions[Y[graphData.edges[e].node1]]++;

            // particao do no 2 (se ja alocado)
            if (Y[graphData.edges[e].node2] >= 0)
                partitions[Y[graphData.edges[e].node2]]++;
        }

        int bestPartition = 0;
        int maior = -1;
        for (unsigned int k=0; k<partitions.size(); k++)
        {
            if (partitions[k] > maior)
            {
                maior = partitions[k];
                bestPartition = k;
            }
            else
            if (partitions[k] == maior)
            {
                // desempatar considerando a menor capacidade ocupada em termos de nos
                if (capP[k] < capP[bestPartition])
                {
                    maior = partitions[k];
                    bestPartition = k;
                }
            }
        }

        Y[i] = bestPartition;
        capP[bestPartition]++; 
    }


    // printf("\n\nSolucao:\nX = ");
    // for(unsigned int i=0; i<edges.size(); i++)
    //     printf("%d ", X[i]);
    // printf("\nY = ");
    // for(unsigned int i=0; i<nodes.size(); i++)
    //     printf("%d ", Y[i]);
    // getchar();

    // calculate fitness
    s.ofv = CalculateOFV(X,Y);

    // tentar trocar a posicao das ultimas arestas
    // double ofvLine = s.ofv;
    // for (int i=(graphData.numEdges/K * K); i<graphData.numEdges; i++)
    // // for (int i=0; i<graphData.numEdges-1; i++)
    // {
    //     for (int j=0; j<graphData.numEdges/K * K; j++)
    //     // for (int j=0; j<graphData.numEdges; j++)
    //     {
    //         int oldKi = X[i];
    //         int oldKj = X[j];

    //         if (X[i] != X[j])
    //         {
    //             // trocar a particao
    //             X[i] = oldKj;
    //             X[j] = oldKi;

    //             // calcular fitness
    //             ofvLine = CalculateOFV(X,Y);

    //             // atualizar ofv se melhorou
    //             if (ofvLine < s.ofv)
    //                 s.ofv = ofvLine;

    //             // retornar o valor antigo caso contrario
    //             else
    //             {
    //                 X[i] = oldKi;
    //                 X[j] = oldKj;
    //             }
    //         }
    //     }
    // }

    // busca local - tentar trocar a posicao dos vertices
    // double ofvLine = s.ofv;
    // for (int i=0; i<graphData.numNodes; i++)
    // {
    //     for (int k=0; k<K; k++)
    //     {
    //         int oldK = Y[i];
    //         if (Y[i] != k)
    //         {
    //             // trocar a particao
    //             Y[i] = k;

    //             // calcular fitness
    //             ofvLine = CalculateOFV(X,Y);

    //             // atualizar ofv se melhorou
    //             if (ofvLine < s.ofv)
    //                 s.ofv = ofvLine;

    //             // retornar o valor antigo caso contrario
    //             else
    //                 Y[i] = oldK;
    //         }
    //     }
    // }


    // print the solution in the screen
    if (debug && print)
    {
        // para cada particao
        for (int k=0; k<K; k++)
        {
            int capOcupada = 0;
            // imprimir as arestas alocadas na particao k
            printf("\nk = %d: ",k);
            for (int j=0; j<graphData.numEdges; j++)
            {
                if (X[j] == k)
                {
                    printf("%d ", j);    
                    capOcupada++;
                }
            }
            printf("\t #E = %d ",capOcupada);
        }

        printf("\n");
        for (int i=0; i<graphData.numNodes; i++)
        {
            printf("%d ", Y[i]);
        }
    }

    // print the solution in a file
    if (!debug && print)
    {
        fprintf(arqSol,"\nSolucao\nX = ");
        for (int i=0; i<graphData.numEdges; i++)
        {
            fprintf(arqSol,"%d ", X[i]);
        }

        fprintf(arqSol,"\nY = ");
        for (int i=0; i<graphData.numNodes; i++)
        {
            fprintf(arqSol,"%d ", Y[i]);
        }

        // para cada particao
        fprintf(arqSol,"\nParticoes");
        for (int k=0; k<K; k++)
        {
            int capOcupada = 0;
            // imprimir as arestas alocadas na particao k
            fprintf(arqSol,"\nk = %d: ",k);
            for (int j=0; j<graphData.numEdges; j++)
            {
                if (X[j] == k)
                {
                    fprintf(arqSol,"%d ", j);    
                    capOcupada++;
                }
            }
            fprintf(arqSol,"\t #E = %d: ",capOcupada);
        }
    }

    return s.ofv;
}

/************************************************************************************
 Method: Decoders
 Description: mapping the random-key solution into a problem solution
*************************************************************************************/
double DecodeOld(TSol s) 
{
    // create an initial list of edges
    std::vector <int> edges(graphData.numEdges);  
    for (int j = 0; j < graphData.numEdges; j++){ edges[j] = j;}

    // sort the problem vector based on the values in the rk vector
    std::sort(edges.begin(), edges.end(), [&s](int i1, int i2) {
        return s.rk[i1] < s.rk[i2];
    });

    // create an initial list of nodes
    std::vector <int> nodes(graphData.numNodes);  
    for (int j = 0; j < graphData.numNodes; j++){ nodes[j] = j;}

    // sort the problem vector based on the values in the rk vector
    // std::sort(nodes.begin(), nodes.end(), [&s](int i1, int i2) {
    //     return s.rk[i1+graphData.numEdges] < s.rk[i2+graphData.numEdges];
    // });

    // problem solution
    std::vector <int> X(graphData.numEdges,-1);    // edges allocation

    std::vector <int> Y(graphData.numNodes,-1); // nodes allocation

    // alocar cada edge em uma particao
    std::vector <int> capP(K,0);                // capacidade das particoes
    // int kCurr = 0;
    int limite = graphData.numEdges/K;
    for(int i=0; i<graphData.numEdges; i++)
    {
        // alocar as primeiras K arestas em K particoes
        // if (i< K || i < graphData.numEdges * 0.3)
        if (i< K || i < graphData.numEdges * s.rk[n-1])
        // if (i< K)
        {
            X[edges[i]] = i%K;
            capP[i%K]++;
        }
        // alocar as demais arestas na melhor particao que tenha capacidade
        else
        {
            std::vector <int> partitions(K,0);

            // aresta
            int e = edges[i];

            // nos da aresta e
            int n1 = graphData.edges[e].node1;
            int n2 = graphData.edges[e].node2;
        
            // verificar n1
            for (unsigned int j=0; j<nodesEdges[n1].size(); j++)
            {
                // outras arestas que tem o no n1
                int eAux = nodesEdges[n1][j];

                // particao da aresta eAux ja foi alocada
                if (X[eAux] >= 0)
                {
                    partitions[X[eAux]]++;
                }
            }

            // verificar n2
            for (unsigned int j=0; j<nodesEdges[n2].size(); j++)
            {
                // outras arestas que tem o no n2
                int eAux = nodesEdges[n2][j];

                // particao da aresta eAux ja foi alocada
                if (X[eAux] >= 0)
                {
                    partitions[X[eAux]]++;
                }
            }

            // inserir na melhor particao que tenha capacidade
            int bestPartition = -1;
            int aux = -1;
            // int limite = graphData.numEdges/K;
            for (int k=0; k<K; k++)
            {
                if ((partitions[k] >= aux) && (capP[k] < limite))
                {
                    aux = partitions[k];
                    bestPartition = k;
                }
            }

            // nao ha particao com capacidade
            if (bestPartition == -1)
            {
                // aumentar o limite e inserir na melhor particao
                limite++;

                bestPartition = -1;
                aux = -1;
                
                for (int k=0; k<K; k++)
                {
                    if ((partitions[k] > aux) && (capP[k] < limite))
                    {
                        aux = partitions[k];
                        bestPartition = k;
                    }
                }
            }

            X[edges[i]] = bestPartition;
            capP[bestPartition]++;
        }
    }

    // alocar os nodes
    for(int i=0; i<graphData.numNodes; i++)
    {
        // encontrar a particao mais comum entre as arestas e nodes conectados a este ponto
        std::vector <int> partitions(K,0);
    
        for (unsigned int j=0; j<nodesEdges[nodes[i]].size(); j++)
        {
            // aresta
            int e = nodesEdges[nodes[i]][j];

            // particao da aresta
            partitions[X[e]]++;

            // particao do no 1 (se ja alocado)
            // if (Y[graphData.edges[e].node1] >= 0)
            //     partitions[Y[graphData.edges[e].node1]]++;

            // particao do no 2 (se ja alocado)
            // if (Y[graphData.edges[e].node2] >= 0)
            //     partitions[Y[graphData.edges[e].node2]]++;
        }

        int bestPartition = 0;
        int aux = -1;
        for (unsigned int k=0; k<partitions.size(); k++)
        {
            if (partitions[k] > aux)
            {
                aux = partitions[k];
                bestPartition = k;
            }
        }

        Y[nodes[i]] = bestPartition;
    }


    // printf("\n\nSolucao:\nX = ");
    // for(unsigned int i=0; i<edges.size(); i++)
    //     printf("%d ", X[i]);
    // printf("\nY = ");
    // for(unsigned int i=0; i<nodes.size(); i++)
    //     printf("%d ", Y[i]);
    // getchar();

    // calculate fitness
    s.ofv = CalculateOFV(X,Y);

    // tentar trocar a posicao das ultimas arestas
    // double ofvLine = s.ofv;
    // for (int i=(graphData.numEdges/K * K); i<graphData.numEdges; i++)
    // // for (int i=0; i<graphData.numEdges-1; i++)
    // {
    //     for (int j=0; j<graphData.numEdges/K * K; j++)
    //     // for (int j=0; j<graphData.numEdges; j++)
    //     {
    //         int oldKi = X[i];
    //         int oldKj = X[j];

    //         if (X[i] != X[j])
    //         {
    //             // trocar a particao
    //             X[i] = oldKj;
    //             X[j] = oldKi;

    //             // calcular fitness
    //             ofvLine = CalculateOFV(X,Y);

    //             // atualizar ofv se melhorou
    //             if (ofvLine < s.ofv)
    //                 s.ofv = ofvLine;

    //             // retornar o valor antigo caso contrario
    //             else
    //             {
    //                 X[i] = oldKi;
    //                 X[j] = oldKj;
    //             }
    //         }
    //     }
    // }

    // // tentar trocar a posicao dos vertices
    // double ofvLine = s.ofv;
    // for (int i=0; i<graphData.numNodes; i++)
    // {
    //     for (int k=0; k<K; k++)
    //     {
    //         int oldK = Y[i];
    //         if (Y[i] != k)
    //         {
    //             // trocar a particao
    //             Y[i] = k;

    //             // calcular fitness
    //             ofvLine = CalculateOFV(X,Y);

    //             // atualizar ofv se melhorou
    //             if (ofvLine < s.ofv)
    //                 s.ofv = ofvLine;

    //             // retornar o valor antigo caso contrario
    //             else
    //                 Y[i] = oldK;
    //         }
    //     }
    // }


    // print the solution in the screen
    if (debug && print)
    {
        // para cada particao
        for (int k=0; k<K; k++)
        {
            int capOcupada = 0;
            // imprimir as arestas alocadas na particao k
            printf("\nk = %d: ",k);
            for (int j=0; j<graphData.numEdges; j++)
            {
                if (X[j] == k)
                {
                    printf("%d ", j);    
                    capOcupada++;
                }
            }
            printf("\t #E = %d ",capOcupada);
        }

        printf("\n");
        for (int i=0; i<graphData.numNodes; i++)
        {
            printf("%d ", Y[i]);
        }
    }

    // print the solution in a file
    if (!debug && print)
    {
        fprintf(arqSol,"\nSolucao\nX = ");
        for (int i=0; i<graphData.numEdges; i++)
        {
            fprintf(arqSol,"%d ", X[i]);
        }

        fprintf(arqSol,"\nY = ");
        for (int i=0; i<graphData.numNodes; i++)
        {
            fprintf(arqSol,"%d ", Y[i]);
        }

        // para cada particao
        fprintf(arqSol,"\nParticoes");
        for (int k=0; k<K; k++)
        {
            int capOcupada = 0;
            // imprimir as arestas alocadas na particao k
            fprintf(arqSol,"\nk = %d: ",k);
            for (int j=0; j<graphData.numEdges; j++)
            {
                if (X[j] == k)
                {
                    fprintf(arqSol,"%d ", j);    
                    capOcupada++;
                }
            }
            fprintf(arqSol,"\t #E = %d: ",capOcupada);
        }
    }

    return s.ofv;
}

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem(){
    graphData.edges.clear();
    nodesEdges.clear();
}

#endif