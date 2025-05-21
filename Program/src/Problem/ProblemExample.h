// *******************************************************************
//      file with specific functions to solve a Problem
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Variables declared in main.cpp
extern int n;                               // size of cromossoms

//----------------- DEFINITION OF PROBLEM SPECIFIC TYPES___ -----------------------


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

    FILE *arq;
    arq = fopen(name,"r");

    if (arq == NULL)
    {
        printf("\nERROR: File (%s) not found!\n",name);
        getchar();
        exit(1);
    }

    // => read data
    

    fclose(arq);

    // define the size of the solution vector
    n = 1;
}


/************************************************************************************
 Method: Decoders
 Description: mapping the random-key solution into a problem solution
*************************************************************************************/
double Decoder(TSol s) 
{
    // create a solution of the problem

    // calculate fitness
    s.ofv = 0;
    return s.ofv;
}

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem(){}

#endif