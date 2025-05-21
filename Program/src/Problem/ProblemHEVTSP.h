// *******************************************************************
//      file with specific functions to solve a HEVTSP
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Variables declared in main.cpp
extern int n;                               // size of cromossoms

//----------------- DEFINITION OF TYPES OF PROBLEM SPECIFIC -----------------------

#define INFINITO 999999999
#define BAD_SOLUTION 100000
#define FACTOR_MULT 100
#define NUM_MODES 4
#define C 1
#define CC 2
#define E 3
#define B 4

// #include <vector>
#include <iostream>
#include <stdio.h>
#include<stdlib.h>
#include <string>
#include <string.h>
#include <sstream>
// #include "Data.h"

// using namespace std;

// extern int n;



struct TNode {
    int node;
    int serviceTime;
    int timeWindow_start;
    int timeWindow_end;
    int departure_time;
};

struct TEdge {
    int nodeStart; // maybe not needed
    int nodeEnd;
    int mode;
    double cost;
    double time;
};

struct TInfo {
    int node;
    int pos;
    int mode;
    double arrival_time;
    double beginning_service;
    double end_service;
    double wait_time;
    double travel_time;
    double travel_cost;
    double service_time;
};



//------ DEFINITION OF GLOBAL CONSTANTS AND VARIABLES OF SPECIFIC PROBLEM  --------

///////////////////////////////////////////////////////
//              HEV-TSPTW given values              //
//////////////////////////////////////////////////////
static int nCustomers;
// Battery values adapted do Watt/min
static double initialCharge;
static double maxCharge;
static double chargeRate;
static double dischargeRate;
// Time values all given in minutes
static std::vector<int> serviceTime;
static std::vector<int> timeWindows;
static std::vector<std::vector<double>> c_timeMatrix;
static std::vector<std::vector<double>> cc_timeMatrix;
static std::vector<std::vector<double>> e_timeMatrix;
static std::vector<std::vector<double>> b_timeMatrix;
static std::vector<std::vector<double>> c_costMatrix;
static std::vector<std::vector<double>> cc_costMatrix;
static std::vector<std::vector<double>> e_costMatrix;
static std::vector<std::vector<double>> b_costMatrix;


//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------

void splitStringIntoIntVector(std::string line, std::vector<int> &vector);

void initializeMatrix(std::vector<std::vector<double>> &matrix, int n);

void splitAllModesMatrix(char line[], int n, std::vector<std::vector<double>> &m1, std::vector<std::vector<double>> &m2, std::vector<std::vector<double>> &m3, std::vector<std::vector<double>> &m4);

/////////////////////////////////////////////////////////////////////////////////////
//                      Utils auxiliar function                                    //
/////////////////////////////////////////////////////////////////////////////////////
bool sortInfoByCost(const TInfo &lhs, const TInfo &rhs) { return lhs.travel_cost < rhs.travel_cost; }

int getKeyInterval(float val) {
    if(val < 0.5) {
        if(val < 0.25) return C;
        else return CC;
    } else {
        if(val < 0.75) return E;
        else return B;
    }
}

void splitStringIntoIntVector(std::string line, std::vector<int> &vector) {
    std::stringstream ss(line);
    std::string sub;
    while(ss.good()) {
        getline(ss, sub, ',');
        if(!sub.empty())
            vector.push_back(stoi(sub));
    }
}

void initializeMatrix(std::vector<std::vector<double>> &matrix, int n) {
    for(int i=0; i<n; i++) {
        std::vector<double> aux(n);
        matrix.push_back(aux);
    }
}

void splitAllModesMatrix(char line[], int n, std::vector<std::vector<double>> &m1, std::vector<std::vector<double>> &m2, std::vector<std::vector<double>> &m3, std::vector<std::vector<double>> &m4) {
    std::stringstream ss(line);
    std::string sub;

    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            getline(ss, sub, ',');
            m1[i][j] = stod(sub);
        }
    }
    
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            getline(ss, sub, ',');
            m2[i][j] = stod(sub);
        }
    }
    
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            getline(ss, sub, ',');
            m3[i][j] = stod(sub);
        }
    }

    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            getline(ss, sub, ',');
            m4[i][j] = stod(sub);
        }
    }
}
/////////////////////////////////////////////////////////////////////////////////////

void ReadData(char nameFile[]) {
    char name[200] = "../Instances/";
    FILE* file;
    int lineSize = 100000, lineCount = 1, mSize = 0;
    char line[lineSize];

    strcat(name, nameFile);
    file = fopen(name, "r");

    if (!file) {
        printf("\nERROR: File (%s) not found!\n",name);
        getchar();
        exit(1);
    } 
    // read all values from file that will be used
    while(fgets(line, lineSize, file)) {
        switch (lineCount) {
            case 5:
                nCustomers = std::stoi(line);
                mSize = nCustomers + 1;
            break;
        
            case 8:
                splitStringIntoIntVector(line, serviceTime);
            break;

            case 11:
                initialCharge = std::stoi(line);
            break;

            case 14:
                maxCharge = std::stoi(line);
            break;

            case 17:
                chargeRate = std::stoi(line) / 60.0;
            break;

            case 20:
                dischargeRate = std::stoi(line) / 60.0;
            break;

            case 29:
                initializeMatrix(c_timeMatrix, mSize);
                initializeMatrix(cc_timeMatrix, mSize);
                initializeMatrix(e_timeMatrix, mSize);
                initializeMatrix(b_timeMatrix, mSize);
                splitAllModesMatrix(line, mSize, c_timeMatrix, cc_timeMatrix, e_timeMatrix, b_timeMatrix);
            break;

            case 32:
                initializeMatrix(c_costMatrix, mSize);
                initializeMatrix(cc_costMatrix, mSize);
                initializeMatrix(e_costMatrix, mSize);
                initializeMatrix(b_costMatrix, mSize);
                splitAllModesMatrix(line, mSize, c_costMatrix, cc_costMatrix, e_costMatrix, b_costMatrix);
            break;

            case 35:
                splitStringIntoIntVector(line, timeWindows);
            break;

            default:
            break;
        }

        lineCount++;
    }

    // The keys vector doesn't include the depot node, but the path to return to depot needs to be included, making nCustomers nodes and nCustomers+1 edges
    n = (2 * nCustomers) + 1; 


    // printf("\n\nLeitura concluida..."); getchar();

    fclose(file);
}

double getCost(int nodeA, int nodeB, int mode) {
    switch (mode) {
        case C:
            return c_costMatrix[nodeA][nodeB];         
        break;
            
        case CC:
            return cc_costMatrix[nodeA][nodeB];
        break;

        case E:
            return e_costMatrix[nodeA][nodeB];
        break;

        case B:
            return b_costMatrix[nodeA][nodeB];
        break;

        default:
            return INFINITO;
        break;
    }
}

double getTime(int nodeA, int nodeB, int mode) {
        switch (mode) {
        case C:
            return c_timeMatrix[nodeA][nodeB];         
        break;
            
        case CC:
            return cc_timeMatrix[nodeA][nodeB];
        break;

        case E:
            return e_timeMatrix[nodeA][nodeB];
        break;

        case B:
            return b_timeMatrix[nodeA][nodeB];
        break;

        default:
            return INFINITO;
        break;
    }
}


/************************************************************************************
 Method: Decoders
 Description: users need to implement at least one decoder, DecK (K = [1,2,3,4,5])
*************************************************************************************/

double Decoder(TSol s) 
{
    // Constructed Path Variables 
    static std::vector<TNode> nodes;
    static std::vector<TEdge> edges;

    s.ofv = 0;

    std::vector <int> sol(n);  
    for(int i=0; i<n; i++) {
        if(i < nCustomers) {
            sol[i] = i + 1;    
        } else {
            sol[i] = getKeyInterval(s.rk[i]);
        }
    }

    // sort the problem vector based on the values in the rk vector
    std::sort(sol.begin(), sol.begin() + nCustomers, [&s](int i1, int i2) {
        return s.rk[i1] < s.rk[i2];
    });

    nodes.clear();
    nodes.resize(0);
    edges.clear();
    edges.resize(0);

    TNode auxNode;
    TEdge auxEdge;
    TInfo auxPosInfo;
    std::vector<TNode>::iterator nodePos;
    std::vector<TEdge>::iterator edgePos;
    std::vector<TInfo> infoDataVector;
    std::vector<double> auxVector;
    int index, nodeA, nodeB;
    double auxTime;
   
    // Add depot node
    auxNode.node = 0;
    auxNode.serviceTime = serviceTime[0];
    auxNode.timeWindow_start = timeWindows[0];
    auxNode.timeWindow_end = timeWindows[1];
    auxNode.departure_time = 0;
    nodes.push_back(auxNode);

    // Add customers by the order given
    for(int i=0; i<nCustomers; i++) {
        // For each i that will be inserted, calculate the Info distance/costs to all possible positions to be inserted
        auxPosInfo.node = sol[i];
        auxPosInfo.mode = sol[i + nCustomers];
        auxPosInfo.service_time = serviceTime[auxPosInfo.node];
        
        for(int j=0; j<nodes.size(); j++) {
            auxPosInfo.pos = j+1;
            auxPosInfo.travel_cost = getCost(nodes[j].node, auxPosInfo.node, auxPosInfo.mode);            
            auxPosInfo.travel_time = getTime(nodes[j].node, auxPosInfo.node, auxPosInfo.mode);
            auxPosInfo.arrival_time = nodes[j].departure_time + auxPosInfo.travel_time;
            
            // Arrival (plus service) needs to be smaller than the end of the time window, or else the insertion is for sure not viable
            if ((auxPosInfo.arrival_time + auxPosInfo.service_time) < timeWindows[(auxPosInfo.node*2)+1]) {
                if(auxPosInfo.arrival_time < timeWindows[auxPosInfo.node*2]) {
                    auxPosInfo.beginning_service = timeWindows[auxPosInfo.node*2];
                    auxPosInfo.wait_time = auxPosInfo.beginning_service - auxPosInfo.arrival_time;
                } else {
                    auxPosInfo.beginning_service = auxPosInfo.arrival_time;
                    auxPosInfo.wait_time = 0.0;
                }
                auxPosInfo.end_service = auxPosInfo.beginning_service + auxPosInfo.service_time;

                infoDataVector.push_back(auxPosInfo);
            }
        }
        
        // Select first viable path with smallest cost
        sort(infoDataVector.begin(), infoDataVector.end(), sortInfoByCost);
        int x = 0;
        do {
            auxPosInfo = infoDataVector[x];
            if (auxPosInfo.wait_time > 0.0) {
                x = -1;
            } else if ((auxPosInfo.beginning_service > timeWindows[(auxPosInfo.node*2)+1]) || (auxPosInfo.end_service > timeWindows[(auxPosInfo.node*2)+1])) {
                x++;
            } else {
                x = -1;
            }
        } while(x >= 0 && x < infoDataVector.size());
    
        // Fill nodes info - (check later if these are necessary)
        auxNode.node = auxPosInfo.node;
        auxNode.timeWindow_start = timeWindows[auxNode.node*2];
        auxNode.timeWindow_end = timeWindows[(auxNode.node*2)+1];
        auxNode.serviceTime = serviceTime[auxNode.node];
        auxNode.departure_time = auxPosInfo.end_service;
        
        nodePos = (nodes.begin() + auxPosInfo.pos);
        nodes.insert(nodePos, auxNode);

        auxEdge.mode = auxPosInfo.mode;        
        edgePos = (edges.begin() + auxPosInfo.pos - 1);
        edges.insert(edgePos, auxEdge);

        infoDataVector.clear();
    }

    // Add path back to depot
    auxEdge.mode = sol[n-1];
    edges.push_back(auxEdge);

    // printf("\nTest edges: ");
    // for(int k=0; k<edges.size(); k++) printf("(%d - %d (%d): %.2f|%.2f)", edges[k].nodeStart, edges[k].nodeEnd, edges[k].mode, edges[k].cost, edges[k].time);
    // getchar(); 

    int n = nodes.size();
    double level = initialCharge;
    double factor = 0.0;

    for(int i=0; i<edges.size(); i++) {
        edges[i].nodeStart = nodes[i].node;
        edges[i].nodeEnd = nodes[(i+1)%n].node;
        edges[i].cost = getCost(edges[i].nodeStart, edges[i].nodeEnd, edges[i].mode);
        edges[i].time = getTime(edges[i].nodeStart, edges[i].nodeEnd, edges[i].mode);

        switch (edges[i].mode) {
            case CC:
                level += chargeRate * edges[i].time;
            break;

            case E:
                level -= dischargeRate * edges[i].time;
            break;

            case B:
                level -= dischargeRate * edges[i].time;
            break;
        
            default:
            break;
        }

        if(level > maxCharge) {
            level = maxCharge;
        } else if (level < 0.0) {
            factor += FACTOR_MULT * (level* -1);
            level = 0.0;
        }

        s.ofv += edges[i].cost;
//        printf("\n(%d - %d)(%d) c:%.2f  t:%.2f  || B_level:%.2f  ||  current objF = %.2f", edges[i].nodeStart, edges[i].nodeEnd, edges[i].mode, edges[i].cost, edges[i].time, level, s.ofv);
//        getchar();        
    }

    s.ofv += factor;    
//    printf("\nFactor = %.2f  |  Final Cost = %.2f", factor, s.ofv);
//    getchar();

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
    // dist.clear();
    // node.clear();
}



#endif