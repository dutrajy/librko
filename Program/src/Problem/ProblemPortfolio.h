// *******************************************************************
// file with specific functions to solve a MIXED-INTEGER LIMITED-ASSET MARKOWITZ MODEL
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"

// Variables declared in main.cpp
// extern int n;                                           // size of random-key vector

//----------------- DEFINITION OF PROBLEM SPECIFIC TYPES___ -----------------------

struct TProblemData {
    int n;                                           // size of the RKO vector 

    int nAssets;                                     // number of assets
    std::vector <double> mean_return;                // expected return
    std::vector <double> std_dev;                    // standard deviation of return
    std::vector <std::vector <double> > covariance;  // covariance between assets

    double lambda;                                   // risk-aversion parameter
    double li;                                       // minimum proportion that must be held of asset i
    double ui;                                       // maximum proportion that can be held of asset i
    int KAssets;                                     // desired number of assets
};

//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------

/************************************************************************************
 Method: ReadData
 Description: read the input data
*************************************************************************************/
TProblemData ReadData(char nameTable[])
{ 
    TProblemData data;

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
    data.lambda = 0.3;                             // risk-aversion parameter
    data.li = 0.01;                                // minimum proportion that must be held of asset i
    data.ui = 0.25;                                // maximum proportion that can be held of asset i
    data.KAssets = 10;                             // desired number of assets

    // number of assets
    fscanf(arq, "%d", &data.nAssets);

    // for each asset i: mean return, standard deviation of return
    data.mean_return.clear();
    data.mean_return.resize(data.nAssets);

    data.std_dev.clear();
    data.std_dev.resize(data.nAssets);
    for (int i = 0; i < data.nAssets; i++)
    {
        fscanf(arq, "%lf %lf", &data.mean_return[i], &data.std_dev[i]);
    }
    
    // for all possible pairs of assets: i, j, correlation between asset i and asset j
    data.covariance.clear();
    data.covariance.resize(data.nAssets, std::vector<double>(data.nAssets));
    while (!feof(arq))
    {
        int i, j; 
        double correlation;
    	fscanf(arq, "%d %d %lf", &i, &j, &correlation);
    	
        // calculate covariance between assets i and j 
        data.covariance[i - 1][j - 1] = correlation * data.std_dev[i - 1] * data.std_dev[j - 1];
        data.covariance[j - 1][i - 1] = data.covariance[i - 1][j - 1];
    }
    fclose(arq);

    // define the size of the solution vector
    // n = 2 * nAssets;

    data.n = 2 * data.KAssets;

    // // print
    // for (int i = 0; i < data.nAssets; i++)
    // {
    //     printf("%lf \t %lf\n", data.mean_return[i], data.std_dev[i]);
    // }
    // // for all possible pairs of assets: i, j, covariance between asset i and asset j
    // for (int i = 0; i < data.nAssets; i++)
    // {
    //     printf("\n");
    //     for (int j = 0; j < data.nAssets; j++)
    //     {
    // 	    printf("%lf \t", data.covariance[i][j]);
    //     }
    // }
    // getchar();

    return data;
}

/************************************************************************************
 Method: Decoder for Markowitz MIP
 Description: mapping the random-key solution into a problem solution
*************************************************************************************/
double Decoder(TSol &s, const TProblemData &data) 
{
    std::vector <int> z(data.KAssets, 0);            // vector with K selected assets
    std::vector <double> w(data.KAssets, 0);         // vector with weights of the K selected assets
    std::vector<bool> selected(data.nAssets, false); // vector with already selected assets

    // Define the K assets and
    // Update weights to comply with bounds 
    double sumW = 0.0;
    for (int k = 0; k < data.KAssets; k++)
    {
        // define the index of the assetsList
        int index = 0;
        if (s.rk[k] >= 0 && s.rk[k] < 1.0)
            index = floor(s.rk[k] * data.nAssets);

        // search the next asset not selected
        while (selected[index])
            index = (index + 1) % data.nAssets;

        // define the selected asset
        z[k] = index;
        selected[index] = true;

        // w = lowerBounds x zi + (upperBounds x zi − lowerBounds x zi) ∗ rk[i+nAssets]
        w[k] = data.li + (data.ui - data.li) * s.rk[k+data.KAssets];
        sumW += w[k];
    }

    // Enforce budget constraint, sum(w) = 1
    double penalty = 0;
    for (int k=0; k<data.KAssets; k++)
    {
        w[k] = w[k]/sumW;

        if (w[k] > data.ui)
            penalty += 10000000 * (w[k] - data.ui);

        else if (w[k] < data.li)
            penalty += 10000000 * (data.li - w[k]);
    }                                 


    double risk_term   = 0.0, 
           return_term = 0.0;

    // Calculate of the risk term wT Σ w
    for (int i = 0; i < data.KAssets; i++) {
        for (int j = 0; j < data.KAssets; j++) {
            risk_term += w[i] * data.covariance[z[i]][z[j]] * w[j];
        }

        // Calculate of the return term μT w
        return_term += data.mean_return[z[i]] * w[i];
    }

    // calculate objective function
    double ofv = (data.lambda * risk_term - (1 - data.lambda) * return_term) + penalty; 

    // print
    // if (print)
    // {
    //     if (debug)
    //     {
    //         printf("\nSolution: %lf (%lf, %lf)\n", s.ofv, risk_term, return_term);
    //         for (int i = 0; i < data.KAssets; i++)
    //         {
    //             printf("z[%d] = %d \t w[%d] = %lf\n", i, z[i], i, w[i]);
    //         }
    //     }

    //     if (!debug)
    //     {
    //         char name[256]="../Results/Results_Pareto";
    //         strcat(name,".csv");

    //         FILE *arq;
    //         arq = fopen(name,"a");

    //         if (!arq)
    //         {
    //             printf("\n\nFile not found %s!!!",name);
    //             getchar();
    //             exit(1);
    //         }

    //         fprintf(arq,"\n%d, %.2lf, %lf, %lf, %lf", data.KAssets, data.lambda, s.ofv, return_term, risk_term);

    //         fclose(arq);
    //     }
    // }
    
    z.clear();
    w.clear();
    selected.clear();
    
    return ofv;
}


/*double DecoderAssets(TSol &s) 
{
    // create a solution of the problem
    std::vector <int> z(nAssets, 0);  
    std::vector <double> w(nAssets, 0);  

    // Binarize indicators: z = 1 for KAssets small random-keys
    for (int k=0; k<KAssets; k++)
    {
        double valueRK = 1.0;
        int idRK = 0;
        for (int i = 0; i < nAssets; i++)
        {
            // find the kth small random-key not selected
            if (s.rk[i] < valueRK && z[i] == 0)
            {
                valueRK = s.rk[i];
                idRK = i;
            }
        }
        // select the idRK asset
        z[idRK] = 1;
    }
      
    // Update weights to comply with bounds w = lowerBounds x zi + (upperBounds x zi − lowerBounds x zi) ∗ rk[i+nAssets]
    double sumW = 0.0;
    for (int i=0; i<nAssets; i++)
    {
        if (z[i] == 1)
            w[i] = li + (ui - li) * s.rk[i+nAssets];

        // w[i] = li * z[i] + (ui * z[i] - li * z[i]) * s.rk[i+nAssets];

        sumW += w[i];
    }

    // Enforce budget constraint, sum(w) = 1
    double penalty = 0;
    for (int i=0; i<nAssets; i++)
    {
        w[i] = w[i]/sumW;

        if (w[i] > ui)
            penalty += 10000000 * (w[i] - ui);

        if (w[i] < li && z[i] == 1)
            penalty += 10000000 * (li - w[i]);
    }                                 


    // calculate objective function
    double risk_term = 0.0, return_term = 0.0;

    // Cálculo do termo de risco wT Σ w
    for (int i = 0; i < nAssets; i++) {
        for (int j = 0; j < nAssets; j++) {
                risk_term += w[i] * covariance[i][j] * w[j];
        }

        // Cálculo do termo de retorno μT w
        return_term += mean_return[i] * w[i];
    }

    // Cálculo final da função objetivo
    double ofv = (lambda * risk_term - (1 - lambda) * return_term) + penalty; 

    // // print
    if (print)
    {
        printf("\nSolution: %lf (%lf, %lf)\n", s.ofv, risk_term, return_term);
        for (int i = 0; i < nAssets; i++)
        {
            // printf("z[%d] = %d \t w[%d] = %lf\n", i, z[i], i, w[i]);
            printf("%d %lf \n", i+1, w[i]);
        }
    }
    
    return ofv;
}*/


/************************************************************************************
 Method: Decoder for Markowitz MIP
 Description: mapping the random-key solution into a problem solution
*************************************************************************************/
/*void DecoderBIOBJ(TSol s) 
{
    // create a solution of the problem
    std::vector <int> z(nAssets, 0);  
    std::vector <double> w(nAssets, 0);  

    // Binarize indicators: z = 1 for KAssets small random-keys
    for (int k=0; k<KAssets; k++)
    {
        double valueRK = 1.0;
        int idRK = 0;
        for (int i = 0; i < nAssets; i++)
        {
            // find the kth small random-key not selected
            if (s.rk[i] < valueRK && z[i] == 0)
            {
                valueRK = s.rk[i];
                idRK = i;
            }
        }
        // select the idRK asset
        z[idRK] = 1;
    }
      
    // Update weights to comply with bounds w = lowerBounds x zi + (upperBounds x zi − lowerBounds x zi) ∗ rk[i+nAssets]
    double sumW = 0.0;
    for (int i=0; i<nAssets; i++)
    {
        if (z[i] == 1)
            w[i] = li + (ui - li) * s.rk[i+nAssets];

        // w[i] = li * z[i] + (ui * z[i] - li * z[i]) * s.rk[i+nAssets];

        sumW += w[i];
    }

    // Enforce budget constraint, sum(w) = 1
    double penalty = 0;
    for (int i=0; i<nAssets; i++)
    {
        w[i] = w[i]/sumW;

        if (w[i] > ui)
            penalty += 10000000 * (w[i] - ui);

        if (w[i] < li && z[i] == 1)
            penalty += 10000000 * (li - w[i]);
    }                                 


    // calculate objective function
    double risk_term = 0.0, return_term = 0.0;

    // Cálculo do termo de risco wT Σ w
    for (int i = 0; i < nAssets; i++) {
        for (int j = 0; j < nAssets; j++) {
                risk_term += w[i] * covariance[i][j] * w[j];
        }

        // Cálculo do termo de retorno μT w
        return_term += mean_return[i] * w[i];
    }

    // Cálculo final da função objetivo
    s.ofv = (lambda * risk_term - (1 - lambda) * return_term) + penalty; 

    char name[256]="../Results/Results_Pareto";
	strcat(name,".csv");

	FILE *arq;
    arq = fopen(name,"a");

	if (!arq)
	{
		printf("\n\nFile not found %s!!!",name);
		getchar();
		exit(1);
	}

    fprintf(arq,"\n%d, %.2lf, %lf, %lf, %lf", KAssets, lambda, s.ofv, return_term, risk_term);

	fclose(arq);
}*/

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem(TProblemData &data)
{
    data.mean_return.clear();
    data.std_dev.clear();
    data.covariance.clear();
}

#endif