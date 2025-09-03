#include "../Problem.h"

#include <vector>
#include <iostream>

class KnapsackProblem final: public Problem
{
private:
    int numberOfItems;
    int capacity;
    std::vector<int> weights;
    std::vector<int> values;

public:
    void loadFromFile(char *path)
    {
        FILE *file;
        file = fopen(path, "r");

        if (file == NULL)
        {
            printf("\nERROR: File (%s) not found!\n", path);
            getchar();
            exit(1);
        }

        fscanf(file, "%d", &this->numberOfItems);
        fscanf(file, "%d", &this->capacity);

        this->weights.clear();
        this->weights.resize(this->numberOfItems);

        this->values.clear();
        this->values.resize(this->numberOfItems);

        for (int k = 0; k < this->numberOfItems; k++)
        {
            fscanf(file, "%d", &this->values[k]);
            fscanf(file, "%d", &this->weights[k]);
        }

        this->dimension = this->numberOfItems;
    }

    const double decode(TSol &solution) const override
    {
        std::vector<int> itemsSelected(this->dimension, 0);
        for (int i = 0; i < this->dimension; i++)
        {
            if (solution.rk[i] > 0.5)
                itemsSelected[i] = 1;
        }

        // calculate the objective function value
        int totalValue = 0;
        int totalWeight = 0;
        for (int i = 0; i < this->dimension; i++)
        {
            if (itemsSelected[i] == 1)
            {
                totalValue += this->values[i];
                totalWeight += this->weights[i];
            }
        }

        // Penalty infeasible solutions
        int infeasible = this->capacity < totalWeight ? totalWeight - this->capacity : 0;
        totalValue = totalValue - (100000 * infeasible);

        // Change to minimization problem
        totalValue = totalValue * -1;

        return totalValue;
    }
};
