// *******************************************************************
//      file with specific functions to solve the VRPTWMD2R
// (Problema de Roteamento de Veıculos com Janelas de Tempo, Multiplos 
//      Entregadores, Roteamento em Dois Niveis e Clusterizacao)
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"
#include <vector>
#include <bits/algorithmfwd.h>
#include <math.h>
#include <string.h>
#include <stdio.h>

// Variables declared in main.cpp
extern int n;                               // size of cromossoms

//----------------- DEFINITION OF PROBLEM SPECIFIC TYPES___ -----------------------
int nV;                                                   /* number of vehicles */
int nClu;                                                 /* number of clusters */
int nC;                                                   /* number of customers */
int nN;                                                   /* number of cluster + 1 */
int Q;                                                    /* maximum capacity of each vehicle */
double vd_speed_ratio = 3;                                /* we assume that vehicles speed is 1 and deliverymen speed is 1 / vd_speed_ratio */

int ML = 3;                                               /* Maximum number of deliverymen in each vehicle */
int fv = 1000;                                            /* Fixed cost associated with each vehicle */ 
int fd = 100;                                             /* Fixed cost associated with each deliveryman */
int cv = 10;                                              /* Unitary distance cost of first-level routes (vehicles) */
int cd = 1;                                               /* Unitary distance cost of second-level routes (deliverymen) */

typedef struct {
    int id;
    int xCoord;                                           /* x coordinate of cluster */
    int yCoord;                                           /* y coordinate of cluster */
    int demand;                                           /* demand at each cluster */
    int readyTime;                                        /* instant when the time window opens at each cluster */
    int dueDate;                                          /* instant when the time window closes at each cluster */
    int serviceTime;                                      /* time spent in a cluster */
} Cluster;
static std::vector <Cluster> clusters;					  /* vector of Clusters */

typedef struct {
    int id;
    int xCoord;                                           /* x coordinate of customers */
    int yCoord;                                           /* y coordinate of customers */
    int demand;                                           /* demand at each customer */
    int readyTime;                                        /* instant when the time window opens at each customer */
    int dueDate;                                          /* instant when the time window closes at each customer */
    int serviceTime;                                      /* time spent in a customer */
    int cluster;                                          /* cluster of a customer */
} Customer;
static std::vector <Customer> customers;				  /* vector of Customers */

typedef struct {
    int cluster_index;
    int dist_total;
    std::vector<int> dist_dl;
    std::vector<int> time_dl;
    std::vector<std::vector<int>> route_dl;
} TClu;

typedef struct {
    int num_dl;
    int dist_vehicle;
    int time_vehicle;
    int q;
    std::vector <TClu> clusters;
} Route;

std::vector<std::vector<int>> distance;                   /* distance between any two nodes, (nClu + 1 + nC) x (nClu + 1 + nC) */
std::vector<std::vector<int>> travel_time;                /* travel time between any two nodes, (nClu + 1 + nC) x (nClu + 1 + nC) */

std::vector<int> nC_per_Clu;                              /* number of customers in each cluster (nClu) */
std::vector<std::vector<int>> customers_per_cluster;      /* list of customers ids that belong to each cluster */


//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  --------



//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------

/* Run Floyd-Warshall Algorithm to ensure triangle inequality */
void ensure_triangle_inequality(std::vector<std::vector<int>> &dist)
{
    bool changed = true;
    while (changed)
    {
        changed = false;
        for (unsigned int k = 0; k < dist.size(); k++)
        {
            for (unsigned int i = 0; i < dist.size(); i++)
            {
                if (i != k)
                {
                    for (unsigned int j = 0; j < dist.size(); j++)
                    {
                        if (j != i && j != k)
                        {
                            if (dist.at(i)[k] + dist.at(k)[j] < dist.at(i)[j] - 1.E-4)
                            {
                                dist.at(i)[j] = dist.at(i)[k] + dist.at(k)[j];
                                changed = true;
                            }
                        }
                    }
                }
            }
        }
    }
}

/************************************************************************************
 Method: ReadData
 Description: read the input data
*************************************************************************************/
void ReadData(char nameTable[])
{
    /*
    We defined travel time and distance as the euclidian distance truncated to integer.
    We also ensured that the triangle inequality was valid.
    For deliverymen travel time, the speed scaling was also considered.
    */

    char filename[200] = "../Instances/";
    strcat(filename,nameTable);

    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("Error opening file.\n");
        return;
    }

    char problemName[20];
    
    // Read problem name and general information
    fscanf(file, "%s", problemName);
    fscanf(file, "%d %d", &nClu, &nC);

    nN = nClu + 1;
    
    // printf("Problem Name: %s\n", problemName);
    // printf("Clusters: %d, Customers: %d\n", nClu, nC);

    // Read vehicle information
    fscanf(file, "%*s %*s %*s %d %d", &nV, &Q);
    // printf("Vehicle Number: %d, Capacity: %d\n", nV, Q);

    // Read depot and customer location data
    clusters.clear();
    customers.clear();
    customers.resize(nC);

    // Read cluster section
    fscanf(file, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
    for (int i = 0; i < nClu+1; i++) {
        Cluster cAux;
        fscanf(file, "%d %d %d %d %d %d %d",
               &cAux.id, &cAux.xCoord, &cAux.yCoord,
               &cAux.demand, &cAux.readyTime, &cAux.dueDate,
               &cAux.serviceTime);

        clusters.push_back(cAux);

        // printf("Cluster %d: X = %d, Y = %d, Demand = %d, Ready Time = %d, Due Date = %d, Service Time = %d\n",
        //        clusters[i].id, clusters[i].xCoord, clusters[i].yCoord,
        //        clusters[i].demand, clusters[i].readyTime, clusters[i].dueDate,
        //        clusters[i].serviceTime);
    }

    // Read customer data
    fscanf(file, "%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s %*s");
    for (int i = 0; i < nC; i++) {
        fscanf(file, "%d %d %d %d %d %d %d %d",
               &customers[i].id, &customers[i].xCoord, &customers[i].yCoord,
               &customers[i].demand, &customers[i].readyTime, &customers[i].dueDate,
               &customers[i].serviceTime, &customers[i].cluster);

        // printf("Customer %d: X = %d, Y = %d, Demand = %d, Ready Time = %d, Due Date = %d, Service Time = %d, Cluster = %d\n",
        //        customers[i].id, customers[i].xCoord, customers[i].yCoord,
        //        customers[i].demand, customers[i].readyTime, customers[i].dueDate,
        //        customers[i].serviceTime, customers[i].cluster);
    }



    /* Create number of customers in each cluste */
    nC_per_Clu.clear();
    nC_per_Clu.resize(nClu + 1);

    nC_per_Clu[0] = 1;        // Cluster 0 is equivalent to the depot

    for (int i = 1; i < nClu + 1; i++)
    {
        for (int j = 0; j < nC; j++)
        {
            if (customers[j].cluster == i)
            {
                nC_per_Clu[i]++;
            }
        }
    }

    // for (int i = 0; i < nClu+1; i++)
    // {
    //     printf("\nCluster %d: %d ", i, nC_per_Clu[i]);
    // }

    /* Define list of customers ids in each cluster */
    customers_per_cluster = std::vector<std::vector<int>>(nClu + 1, std::vector<int>());
    customers_per_cluster[0].push_back(0);

    for (int i = 1; i < nClu + 1; i++)
    {
        for (int j = 0; j < nC; j++)
        {
            if (customers[j].cluster == i)
            {
                customers_per_cluster[i].push_back(j);
            }
        }
    }

    // for (int i = 1; i < nClu + 1; i++)
    // {
    //     printf("\nCluster %d [%d, %d]: \t", i, clusters[i].readyTime, clusters[i].dueDate);
    //     for (unsigned int j = 0; j < customers_per_cluster[i].size(); j++)
    //     {
    //         int pos = customers_per_cluster[i][j];
    //         printf("%d [%d, %d] \t", customers_per_cluster[i][j], customers[pos].readyTime, customers[pos].dueDate);
    //     }
    // }

    distance.resize(nN + nC, std::vector<int>(nN + nC, 0));
    travel_time.resize(nN + nC, std::vector<int>(nN + nC, 0));

    /* Fill the matrix of travel time */
    for (int i = 0; i < nN + nC; i++)
    {
        for (int j = 0; j < nN + nC; j++)
        {
            double x;
            double y;

            if (i < nN && j < nN) // i e j sao clusters
            {
                x = clusters[i].xCoord - clusters[j].xCoord;
                y = clusters[i].yCoord - clusters[j].yCoord;
            }
            else
            if (i >= nN && j < nN) // i eh customer e j eh cluster
            {
                x = customers[i-nN].xCoord - clusters[j].xCoord;
                y = customers[i-nN].yCoord - clusters[j].yCoord;
            }
            else
            if (i < nN && j >= nN) // i eh cluster e j eh customer
            {
                x = clusters[i].xCoord - customers[j-nN].xCoord;
                y = clusters[i].yCoord - customers[j-nN].yCoord;
            }
            else    // ambos sao customers
            {
                x = customers[i-nN].xCoord - customers[j-nN].xCoord;
                y = customers[i-nN].yCoord - customers[j-nN].yCoord;
            }

            int aux = (int)(sqrt(x * x + y * y));
            distance[i][j] = aux;
        }
    }

    // **** preciso verificar para cada sub matriz (clusters, customers of cluster 1, ...) ****
    std::vector<std::vector<int>> mat;
    std::vector<int> matID;

    // matriz clusters
    mat.clear();
    mat.resize(nN, std::vector<int>(nN, 0));
    mat.resize(nN);

    matID.clear();
    matID.resize(nN);

    for (int i=0; i<nN; i++){
        matID[i] = i;
    }

    for (int i=0; i<nN; i++){
        for (int j=0; j<nN; j++){
            mat[i][j] = distance[matID[i]][matID[j]];
        }
    }

    ensure_triangle_inequality(mat);

    for (int i=0; i<nN; i++){
        // printf("\n");
        for (int j=0; j<nN; j++){
            distance[matID[i]][matID[j]] = mat[i][j];
            // printf("%d ", mat[i][j]);
        }
    }

    // matrizes secundarias
    for (int h = 1; h < nN; h++)
    {
        // dimension = numero de customers por cluster + o deposito (cluster i)
        int dim = (int)customers_per_cluster[h].size() + 1;

        mat.clear();
        mat.resize(dim, std::vector<int>(dim, 0));
        mat.resize(dim);

        matID.clear();
        matID.resize(dim);

        matID[0] = h; // cluster h
        for (unsigned int j=0; j<customers_per_cluster[h].size(); j++){
            matID[j+1] = customers_per_cluster[h][j] + nN;
        }

        for (int i=0; i<dim; i++){
            for (int j=0; j<dim; j++){
                mat[i][j] = distance[matID[i]][matID[j]];
                // printf("(%d, %d) ", matID[i], matID[j]);
            }
        }

        ensure_triangle_inequality(mat);

        // printf("\n\nCluster %d:", h);
        for (int i=0; i<dim; i++){
            // printf("\n");
            for (int j=0; j<dim; j++){
                distance[matID[i]][matID[j]] = mat[i][j];
                // printf("%d ", mat[i][j]);
            }
        }
    }

    // for (int i = 0; i < (int)distance.size(); i++)
    // {
    //     // printf("\ni = %d: %d | ", i, (int)distance[i].size());
    //     printf("\n");
    //     for (int j = 0; j < (int)distance[i].size(); j++)
    //     {
    //         printf("%d ", distance[i][j]); 
    //     }
    // }

    for (int i = 0; i < nN + nC; i++)
    {
        for (int j = 0; j < nN + nC; j++)
        {
            // ambos sao clusters
            if (i < nN && j < nN)
                travel_time[i][j] = distance[i][j];
            // um deles eh customer
            else
                travel_time[i][j] = distance[i][j] * vd_speed_ratio;
        }
    }

    // for (int i = 0; i < (int)distance.size(); i++)
    // {
    //     // printf("\ni = %d: %d | ", i, (int)distance[i].size());
    //     printf("\n");
    //     for (int j = 0; j < (int)travel_time[i].size(); j++)
    //     {
    //         printf("%d ", travel_time[i][j]); 
    //     }
    // }

    


    // update the number of vehicles
    if (nClu < nV)
        nV = nClu;
    
    // define n (clusters + customers + vehicles (limitado ao numero de clusters))
    n = nClu + nC + nV;

    fclose(file);
    // getchar();
    
}

/************************************************************************************
 Method: Decoders
 Description: mapping the random-key solution into a problem solution
*************************************************************************************/
double Decoder(TSol s) 
{
    //horario de chegada A_i = D_i-1 + t_i-1,i
    std::vector <int> A_pri(nClu+1);
    std::vector <int> A_sec(nC);
    int A0;

    //horario de inicio de servico B_i = max {e_i, A_i}
    std::vector <int> B_pri(nClu+1);
    std::vector <int> B_sec(nC);
    // int B0, Bn;

    //horario de partida D_i = B_i + s_i
    std::vector <int> D_pri(nClu+1); 
    std::vector <int> D_sec(nC); 

    //tempo de espera antes de iniciar o servico W_i = B_i - A_i
    std::vector <int> W_pri(nClu+1);
    std::vector <int> W_sec(nC);

    //capacidade dos veiculos
    std::vector <int> Q_pri(nClu+1);



    /*create a solution of the problem */ 
    // printf("\n\nRK: ");
    // for (int i = 0; i < n; i++)
    // {
    //     printf("%.3lf, ", s.rk[i]);
    // }

    // create an initial list of clusters
    std::vector <int> sClusters(nClu);  
    for (int j = 0; j < nClu; j++){ sClusters[j] = j+1;}

    // sort the cluster vector based on the values in the rk vector
    std::sort(sClusters.begin(), sClusters.end(), [&s](int i1, int i2) {
        return s.rk[i1-1] < s.rk[i2-1];
    });

    // printf("\nClusters: ");
    // for(int i=0; i<nClu; i++)
    // {
    //     printf("%d, ", sClusters[i]);
    // }

    // create an initial list of customers
    std::vector <int> sCustomers(nC);  
    for (int j = 0; j < nC; j++){ sCustomers[j] = j;}

    // sort the cluster vector based on the values in the rk vector
    std::sort(sCustomers.begin(), sCustomers.end(), [&s](int i1, int i2) {
        return s.rk[i1+nClu] < s.rk[i2+nClu];
    });

    // printf("\n\nCustomers: ");
    // for(int i=0; i<nC; i++)
    // {
    //     printf("%d, ", sCustomers[i]);
    // }


    /* rotas dos veiculos */
    std::vector <Route> routes;
    routes.clear();
    routes.resize(nV);

    // define number of deliveryman per vehicle
    for (int i = 0; i < nClu; i++)
    {
        routes[i].num_dl = (int)ceil(s.rk[i+nClu+nC]*ML);
        routes[i].q = 0;
        routes[i].dist_vehicle = 0;
        routes[i].time_vehicle = 0;

        //criar a rota com o cluster 0 (depot)
        TClu aux;
        aux.cluster_index = 0;
        aux.dist_total = 0;
        aux.dist_dl.resize(routes[i].num_dl,0);
        aux.time_dl.resize(routes[i].num_dl,0);
        aux.route_dl.clear();

        routes[i].clusters.push_back(aux);
    }

    // printf("\n\nDeliveryman per vehicle: ");
    // for (int i = 0; i < nV; i++)
    // {
    //     printf("%d, ", routes[i].num_dl);
    // }
    
    //horario de inicio de servico na garagem => janela de tempo da garagem
    // B0 = veiculos[k].e;
    // Bn = B0;

    // //horario de partida da garagem
    // D0 = B0;

    // //carga do veiculo k ao sair da garagem
    // Q0 = 0;

    // //calcular a maior carga do veiculo k
    // maiorCap = 0;

    /* construir as rotas incluindo cada um dos clusters */
    for (int i=0; i<nClu; i++)
    {
        int cluCurrent = sClusters[i];

        // encontrar o primeiro veiculo para inserir o cluster corrente (respeitando capacidade e janela de tempo)
        for (int k=0; k<nV; k++)
        {
            // ultimo cluster visitado
            int lastCluster = routes[k].clusters[routes[k].clusters.size()-1].cluster_index;

            // rota tem capacidade disponivel e chega no cluster corrente antes da janela fechar
            if (routes[k].q + clusters[cluCurrent].demand <= Q && routes[k].time_vehicle + travel_time[lastCluster][cluCurrent] <= clusters[cluCurrent].dueDate)
            {
                // atualizar a demanda atendida
                routes[k].q += clusters[cluCurrent].demand;

                // inserir o cluster na rota
                TClu aux;
                aux.cluster_index = cluCurrent;
                aux.dist_dl.resize(routes[i].num_dl,0);
                aux.time_dl.resize(routes[i].num_dl,0);

                // criar a rota dos entregadores no cluster (menor distancia e janela de tempo)
                aux.route_dl.clear();
                aux.route_dl.resize(routes[k].num_dl, std::vector<int>(1, cluCurrent));
                
                int lastCustomer = 0;                                                   // ultimo customer visitado
                A0 = routes[k].time_vehicle += travel_time[lastCluster][cluCurrent];    // horario de chegada ao cluster

                int e = 0;
                for (int j = 0; j < nC; j++)
                {
                    if (customers[sCustomers[j]].cluster == cluCurrent)
                    {
                        // verificar qual o primeiro entregador para atender o customer j
                        for (e=0; e<routes[k].num_dl; e++)
                        {
                            aux.route_dl[e].push_back(sCustomers[j]);

                            if (aux.time_dl[e] == 0){
                                aux.time_dl[e] = travel_time[lastCluster][sCustomers[j]+nClu];
                                aux.dist_dl[e] = distance[lastCluster][sCustomers[j]+nClu];
                            }
                            else
                            {
                                aux.time_dl[e] += travel_time[lastCustomer][sCustomers[j]+nClu];
                                aux.dist_dl[e] = distance[lastCustomer][sCustomers[j]+nClu];
                            }

                            lastCustomer = sCustomers[j];
                            break;
                        }
                        
                    }
                }

                if (e < routes[k].num_dl){
                    aux.time_dl[e] += travel_time[lastCustomer][lastCluster];
                    aux.dist_dl[e] = distance[lastCustomer][lastCluster];
                }

                routes[k].clusters.push_back(aux);

                // atualizar o tempo da rota com a viagem ate o cluster, o maior *** tempo dos entregadores no cluster corrente e o tempo de servico
                routes[k].time_vehicle += travel_time[lastCluster][cluCurrent] + aux.time_dl[e] + clusters[cluCurrent].serviceTime;

                // interromper o loop
                break;
            }
        }
    }
    
    for (int k=0; k<nV; k++)
    {
        printf("\n\nRota %d: \n NumDL = %d \n q = %d \n DistVehicle = %d", k, routes[k].num_dl, routes[k].q, routes[k].dist_vehicle);
        for (unsigned int i = 0; i < routes[k].clusters.size(); i++)
        {
            printf("\n Cluster %d: \n  Dist_dl = %d \n  Route_dl: ", routes[k].clusters[i].cluster_index, routes[k].clusters[i].dist_total);
            for (unsigned int j=0; j< routes[k].clusters[i].route_dl.size(); j++)
            {
                printf("\n    dl %d: ", j);
                for (unsigned int e=0; e< routes[k].clusters[i].route_dl[j].size(); e++)
                    printf("%d ", routes[k].clusters[i].route_dl[j][e]);
            }
        }
        
    }

    // calculate fitness
    s.ofv = 0;

    getchar();

    return s.ofv;
}

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem(){
    clusters.clear();
    customers.clear();
    distance.clear();
    travel_time.clear();
    nC_per_Clu.clear();
    customers_per_cluster.clear();
}

#endif



/************************************************************************************
*** M�todo: CalcularFO(TSol s)                                                    ***
*** Fun��o: Calcula o valor da fun��o objetivo da solu��o s                       ***
// *************************************************************************************/
// float CalcularFO(TSol sol)
// {
//     unsigned int temp = 0;
//     for(int k=0; k<m; k++)
//         temp += sol.rota[k].size();
//     if (temp < n)
//         printf("\n*********************************Erro: ponto nao realocado (%d)",temp);

//     //pesos
//     int omega[5] = {8,0,3,1,0}; //(distancia, #V, rota, percurso, espera)
//     int teta[5] = {1500, 1500, 0, 1500, 1500}; //(rotaMax, percursoMax, esperaMax, capMax, janMax)

//     float distTotal = 0,
//           rotaTotal = 0,
//           percursoTotal = 0,
//           esperaTotal = 0,
//           inv_Rota = 0,
//           inv_Percurso = 0,
//           inv_Espera = 0,
//           inv_Janela = 0;

//     int numVeiculos = 0,
//         inv_Cap = 0,
//         maiorCap = 0;

//     //horario de chegada A_i = D_i-1 + t_i-1,i
//     float A[2*n], An;

//     //horario de inicio de servi�o B_i = max {e_i, A_i}
//     float B[2*n], B0, Bn;

//     //horario de partida D_i = B_i + s_i
//     float D[2*n], D0;

//     //tempo de espera antes de iniciar o servi�o W_i = B_i - A_i
//     float W[2*n];

//     //capacidade dos veiculos
//     int Q[2*n], Q0;

//     //tempo de viagem R_i = B_n+1 - D_i
//     float R[n];

//     //calcular os custos para cada veiculo
//     for (int k=0; k<m; k++)
//     {
//         //horario de inicio de servico na garagem => janela de tempo da garagem
//         B0 = veiculos[k].e;
//         Bn = B0;

//         //horario de partida da garagem
//         D0 = B0;

//         //ocupacao do veiculo k ao sair da garagem
//         Q0 = 0;

//         //calcular a maior ocupacao do veiculo k
//         maiorCap = 0;

//         //percorrer a rota do veiculo k
//         for (unsigned int i=0; i<sol.rota[k].size(); i++)
//         {
//             //calcular a distancia da garagem ao primeiro cliente
//             if (i == 0)
//             {
//                 distTotal += G[k][sol.rota[k][i]];
//                 numVeiculos++; //aumentar o numero de veiculos utilizados

//                 A[sol.rota[k][i]] = D0 + Gt[k][sol.rota[k][i]];
//                 B[sol.rota[k][i]] = MAX(pontos[sol.rota[k][i]].e, A[sol.rota[k][i]]);
//                 D[sol.rota[k][i]] = B[sol.rota[k][i]] + pontos[sol.rota[k][i]].s;
//                 W[sol.rota[k][i]] = B[sol.rota[k][i]] - A[sol.rota[k][i]];
//                 Q[sol.rota[k][i]] = Q0 + pontos[sol.rota[k][i]].q;

//                 if (Q[sol.rota[k][i]] > maiorCap)
//                     maiorCap = Q[sol.rota[k][i]];
//             }

//             //calcular a distancia de i-1 para i (at� o �ltimo cliente da rota)
//             else
//             {
//                 distTotal += Dist[sol.rota[k][i-1]][sol.rota[k][i]];

//                 A[sol.rota[k][i]] = D[sol.rota[k][i-1]] + Tempo[sol.rota[k][i-1]][sol.rota[k][i]];
//                 B[sol.rota[k][i]] = MAX(pontos[sol.rota[k][i]].e, A[sol.rota[k][i]]);
//                 D[sol.rota[k][i]] = B[sol.rota[k][i]] + pontos[sol.rota[k][i]].s;
//                 W[sol.rota[k][i]] = B[sol.rota[k][i]] - A[sol.rota[k][i]];
//                 Q[sol.rota[k][i]] = Q[sol.rota[k][i-1]] + pontos[sol.rota[k][i]].q;

//                 if (Q[sol.rota[k][i]] > maiorCap)
//                     maiorCap = Q[sol.rota[k][i]];
//             }

//             //calcular a distancia do ultimo cliente at� a garagem
//             if (i == sol.rota[k].size()-1)
//             {
//                 distTotal += G[k][sol.rota[k][i]];
//                 An = D[sol.rota[k][i]] + Gt[k][sol.rota[k][i]];
//                 Bn = An;
//             }
//         }

//         //calcular a duracao das viagens dos clientes do veiculo k
//         for (unsigned int i = 0; i < sol.rota[k].size(); i++)
//         {
//             if (sol.rota[k][i] < n)
//             {
//                 R[sol.rota[k][i]] = B[n+sol.rota[k][i]] - D[sol.rota[k][i]];
//                 percursoTotal += R[sol.rota[k][i]];

//                 //verificar se o percurso do cliente eh maior que o limite maximo
//                 if (R[sol.rota[k][i]] > Rmax)
//                     inv_Percurso += R[sol.rota[k][i]] - Rmax;
//             }

//             //calcular a espera do veiculo k em todos os pontos
//             esperaTotal += W[sol.rota[k][i]];

//             //verificar se estrapola as janelas de tempo
//             if (pontos[sol.rota[k][i]].e > B[sol.rota[k][i]])
//                 inv_Janela += pontos[sol.rota[k][i]].e - B[sol.rota[k][i]];
//             if (pontos[sol.rota[k][i]].l < B[sol.rota[k][i]])
//                 inv_Janela += B[sol.rota[k][i]] - pontos[sol.rota[k][i]].l;
//         }

//         //calcular a duracao da rota
//         rotaTotal += Bn - D0;

//         //verificar se a duracao da rota eh maior que a duracao maxima
//         if ((Bn - D0) > Tmax)
//             inv_Rota += (Bn - D0) - Tmax;

//         //verificar se estrapolou a capacidade do veiculo
//         if (maiorCap > veiculos[k].cap)
//             inv_Cap += maiorCap - veiculos[k].cap;

//     }


//     sol.fo = omega[0]*distTotal +
//              omega[1]*numVeiculos +
//              omega[2]*rotaTotal +
//              omega[3]*percursoTotal +
//              omega[4]*esperaTotal +
//               teta[0]*inv_Rota +
//               teta[1]*inv_Percurso +
//               teta[2]*inv_Espera +
//               teta[3]*inv_Cap +
//               teta[4]*inv_Janela;

//     if (debug == 1)
//     {
//          printf("\n\n< ---------------------- DARP ----------------------- >\n");
//          printf("Numero de veiculos disponiveis.................: %d\n",m);
//          printf("Capacidade dos veiculos........................: %d\n",veiculos[0].cap);
//          printf("Numero de localidades..........................: %d\n",n*2);
//          printf("Numero de clientes.............................: %d\n\n",n);

//          printf("Distancia das rota.............................: %.2f km\n",distTotal);
//          printf("Numero de veiculos utilizados..................: %d \n",numVeiculos);
//          printf("Duracao das rotas..............................: %.2f minutos\n",rotaTotal);
//          printf("Tempo total de viagem..........................: %.2f minutos\n",percursoTotal);
//          printf("Tempo total de espera..........................: %.2f minutos\n\n",esperaTotal);

//          printf("Violacoes nas cargas dos veiculos..............: %d\n",inv_Cap);
//          printf("Violacoes nas janelas de tempo.................: %.2f\n",inv_Janela);
//          printf("Violacoes nas duracoes das viagens.............: %.2f\n",inv_Percurso);
//          printf("Violacoes no tempo de espera...................: %.2f\n",inv_Espera);
//          printf("Violacoes nas duracoes das rotas...............: %.2f\n",inv_Rota);

//          for (int k=0; k<m; k++)
//          {
//              printf("\n\nVeiculo %d => ", k);
//              for (unsigned int i=0; i<sol.rota[k].size(); i++)
//              {
//                  printf("\n%2d (a: %.2f; \t w: %.2f; \t q: %d) ", sol.rota[k][i], A[sol.rota[k][i]],  W[sol.rota[k][i]], Q[sol.rota[k][i]]);
//              }
//          }
//     }

//     //aplicar a heuristica de atraso

//     return sol.fo;
// }