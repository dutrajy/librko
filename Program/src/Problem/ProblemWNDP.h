// ********************************************************************
// file with specific functions to solve a Water Network Design Problem
// ********************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// #include <epanet2.h>

// Variables declared in main.cpp
extern int n;                               // size of cromossoms

//----------------- DEFINITION OF PROBLEM SPECIFIC TYPES___ -----------------------


//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  --------
int N;
int A;
int K;
double alphaW;
double pr;
int roughness;
double gammaW;
std::vector<int> Le;
std::vector<double> Commercial_D;
std::vector<int> Cost_D;
std::vector<int> Demand_in;
std::vector<double> Demand;
std::vector<int> H_min;

//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------

/************************************************************************************
 Method: ReadData
 Description: read the input data
*************************************************************************************/
void ReadData(char nameTable[])
{ 
    char name[200] = "../Instances/";
    strcat(name,nameTable);

        FILE *arq;
    arq = fopen(name, "r");

    if (arq == NULL) {
        printf("\nERROR: File (%s) not found!\n", name);
        getchar();
        exit(1);
    }

    // => read data
    char temp[100];

    fscanf(arq, "%s %d", temp, &N);
    fscanf(arq, "%s %d", temp, &A);
    fscanf(arq, "%s %d", temp, &K);
    fscanf(arq, "%s %lf", temp, &alphaW);
    fscanf(arq, "%s %lf", temp, &pr);
    fscanf(arq, "%s %d", temp, &roughness);
    fscanf(arq, "%s %lf", temp, &gammaW);

    // Read Le
    // fscanf(arq, "%s", temp);
    Le.resize(A);
    for (int j = 0; j < A; j++) {
        fscanf(arq, "%d", &Le[j]);
    }

    // Read Commercial_D
    // fscanf(arq, "%s", temp);
    Commercial_D.resize(K);
    for (int k = 0; k < K; k++) {
        fscanf(arq, "%lf", &Commercial_D[k]);
    }

    // Read Cost_D
    // fscanf(arq, "%s", temp);
    Cost_D.resize(K);
    for (int k = 0; k < K; k++) {
        fscanf(arq, "%d", &Cost_D[k]);
    }

    // Read Demand_in
    // fscanf(arq, "%s", temp);
    Demand_in.resize(N);
    for (int i = 0; i < N; i++) {
        fscanf(arq, "%d", &Demand_in[i]);
    }

    // Read Demand
    // fscanf(arq, "%s", temp);
    Demand.resize(N);
    for (int i = 0; i < N; i++) {
        fscanf(arq, "%lf", &Demand[i]);
    }

    // Read H_min
    // fscanf(arq, "%s", temp);
    H_min.resize(N);
    for (int i = 0; i < N; i++) {
        fscanf(arq, "%d", &H_min[i]);
    }

    // Print values for verification
    printf("\nN: %d", N);
    printf("\nA: %d", A);
    printf("\nK: %d", K);
    printf("\nAlfa: %lf", alphaW);
    printf("\npr: %lf", pr);
    printf("\nRoughness: %d", roughness);
    printf("\nGamma: %lf", gammaW);

    printf("\nLe: ");
    for (int j = 0; j < A; j++) {
        printf("%d ", Le[j]);
    }

    printf("\nCommercial_D: ");
    for (int k = 0; k < K; k++) {
        printf("%lf ", Commercial_D[k]);
    }

    printf("\nCost_D: ");
    for (int k = 0; k < K; k++) {
        printf("%d ", Cost_D[k]);
    }

    printf("\nDemand_in: ");
    for (int i = 0; i < N; i++) {
        printf("%d ", Demand_in[i]);
    }

    printf("\nDemand: ");
    for (int i = 0; i < N; i++) {
        printf("%lf ", Demand[i]);
    }

    printf("\nH_min: ");
    for (int i = 0; i < N; i++) {
        printf("%d ", H_min[i]);
    }

    // getchar();
    fclose(arq);

    // define the size of the solution vector
    n = A;
}


/************************************************************************************
 Method: Decoders
 Description: mapping the random-key solution into a problem solution
*************************************************************************************/
double Decoder(TSol s) 
{
    // create a solution of the problem
    std::vector <int> sol(n);  

    double ofv = 0.0;

    // define the diameters of each arc
    for (int j=0; j<n; j++){
        sol[j] = floor(s.rk[j] * K);
    }

    // # One-objective function that calculates the cost of a solution
    for (int j=0; j<n; j++){
        ofv += Cost_D[sol[j]] * Le[j];
    }



    // // Inicialização do interpretador Python
    // Py_Initialize();

    // // Código Python para importar e usar o pacote WNTR
    // const char *command = "import wntr\n"
    //                       "from wntr import network\n"
    //                       "wn = network.WaterNetworkModel()\n"  // cria um rede vazia
    //                       //# --- Nodes --- #
    //                       "wn.add_reservoir('res', base_head=210, head_pattern=None,coordinates=(2000,2000))\n"
    //                       "wn.add_junction('node2', base_demand=Demand[0], elevation=150.0, coordinates=(1000,2000))\n"
    //                       "wn.add_junction('node3', base_demand=Demand[1], elevation=160.0, coordinates=(0,2000))\n"
    //                       "wn.add_junction('node4', base_demand=Demand[2], elevation=155.0, coordinates=(1000,1000))\n"
    //                       "wn.add_junction('node5', base_demand=Demand[3], elevation=150.0, coordinates=(0,1000))\n"
    //                       "wn.add_junction('node6', base_demand=Demand[4], elevation=165.0, coordinates=(1000,0))\n"
    //                       "wn.add_junction('node7', base_demand=Demand[5], elevation=160.0, coordinates=(0,0))\n"
    //                       //# --- Links --- #
    //                       "wn.add_pipe('pipe1', 'res', 'node2', length=1000, diameter=" sol[0] ", roughness=130, minor_loss=0.0)\n"
    //                       "wn.add_pipe('pipe2', 'node2', 'node3', length=1000, diameter=sol[1], roughness=130, minor_loss=0.0)\n"
    //                       "wn.add_pipe('pipe3', 'node2', 'node4', length=1000, diameter=sol[2], roughness=130, minor_loss=0.0)\n"
    //                       "wn.add_pipe('pipe4', 'node4', 'node5', length=1000, diameter=sol[3], roughness=130, minor_loss=0.0)\n"
    //                       "wn.add_pipe('pipe5', 'node3', 'node5', length=1000, diameter=sol[4], roughness=130, minor_loss=0.0)\n"
    //                       "wn.add_pipe('pipe6', 'node6', 'node4', length=1000, diameter=sol[5], roughness=130, minor_loss=0.0)\n"
    //                       "wn.add_pipe('pipe7', 'node7', 'node6', length=1000, diameter=sol[6], roughness=130, minor_loss=0.0)\n"
    //                       "wn.add_pipe('pipe8', 'node5', 'node7', length=1000, diameter=sol[7], roughness=130, minor_loss=0.0)\n"

    //                       "epanet_sim = wntr.sim.EpanetSimulator(wn)\n"
    //                       "epanet_sim_results = epanet_sim.run_sim()\n"
    //                       "for i in range(8): epanet_sim_results.node['pressure'].values[0][i]\n";


    // // Executa o código Python
    // PyRun_SimpleString(command);

    // // Finalização do interpretador Python
    // Py_Finalize();


    double penalty = randomico(1000, 19000000);
    ofv += penalty;
    

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
    
    // calculate fitness
    return ofv;
}

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem(){}

#endif