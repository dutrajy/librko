// *******************************************************************
//      file with specific functions to solve the VSBPP
// *******************************************************************
#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "../Main/Data.h"
#include <queue>

//----------------- DEFINITION OF PROBLEM SPECIFIC TYPES -----------------------
struct TProblemData
{
    int n;                                      // size of the RKO vector 

    std::vector <int> ck;	                    // vector with cost of bins
    std::vector <int> Wk;	                    // vector with capacity of bins
    std::vector <int> wi;	                    // vector with weigth of items

    std::vector <int> itemsDec;                 // itens em ordem decrescente
    std::vector <int> itemsCre;                 // itens em ordem crescente                

    int nItems;                                  // numero de itens
    int mBins;                                   // numero de bins
    int set;                                     // conjunto da instancia
    int numInstance;                             // numero da instancia dentro do conjunto
    int UBbins;                                  // limitante superior do numero de bins necessarios
    int Wmin;                                    // capacidade do bin com menor capacidade
    int Wmax;                                    // capacidade do bin com maior capacidade
    int Wmean;                                   // capacidade media dos bins
    int custoReducao = 0;                        // custo com a reducao inicial dos itens
};

struct TBin
{
    std::vector <int> items;                   // itens alocados no bin
    int cap;                                   // capacidade ocupada do bin
    int type;                                  // tipo do bin
};

//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------


void Reduction(TProblemData &data)
{
    // create a solution of the problem
    TBin bins;
    std::vector <int> items(data.nItems);
    std::vector<bool> selected(data.nItems, false); // vetor para marcar itens já escolhidos

    // criar os itens
    items.resize(data.nItems);
    for (int j = 0; j < data.nItems; j++){
        items[j] = j;                                  
    }

    // ordenar em ordem decrescente em relacao ao peso
    std::sort(items.begin(), items.end(), [&](int i1, int i2) {
        return data.wi[i1] > data.wi[i2];
    });

    // reduction
    for (int i=0; i<data.nItems; i++)
    {
        // item eh mais pesado que o segundo bin com a maior capacidade
        if (data.wi[items[i]] > data.Wk[data.mBins-2] && data.wi[items[i]] <= data.Wk[data.mBins-1])
        {
            // inserir o item em um bin com a maior capacidade
            bins.cap = data.wi[items[i]];
            bins.type = data.mBins-1;
            bins.items.clear();
            bins.items.push_back(items[i]);

            // tentar preencher o restante deste bin
            for (int j=0; j<data.nItems; j++)
            {
                if (!selected[items[j]])
                {
                    // verificar se posso inserir um item neste bin
                    if (bins.cap + data.wi[items[j]] <= data.Wk[bins.type])
                    {
                        bins.cap += data.wi[items[j]];
                        bins.items.push_back(items[j]);
                    }
                }
            }

            // so marcar como selecionado se a cap for igual ao limite do bin
            if (bins.cap == data.Wk[bins.type])
            // if (bins.cap > 0)
            {
                for (int j = 0; j < (int)bins.items.size(); j++)
                {
                    selected[bins.items[j]] = true;
                }
                data.custoReducao += data.ck[bins.type];  
            }
            else
            {
                bins.cap = 0;
                bins.items.clear();
            }
        }
        else
        {
            break;
        }

        // imprimir o bin
        // printf("\nBin: [%d, %d]: ", bins.cap, bins.type);
        // for (int i = 0; i < (int)bins.items.size(); i++)
        // {
        //     printf("%d ", bins.items[i]);
        // }
        
    }

    // remover os itens selecionados do problema e considerar o custo da reducao
    int sumWeight = 0;
    for (int i=data.nItems-1; i>= 0; i--)
    {
        if (selected[i])
        {
            data.wi.erase(data.wi.begin() + i);
        }
        else
        {
            sumWeight += data.wi[i];
        }
    }

    // printf("\n\nCusto reducao: %d | Itens reduzidos: %d\n\n", data.custoReducao, data.nItems - (int)data.wi.size());

    data.nItems = data.wi.size();
    data.UBbins = (ceil((float)sumWeight/data.Wmin) + 1);

    // printf("\n\nItens: %d, Bins: %d, SomaPeso: %d, PesoMin: %d\n\n", data.nItems, data.UBbins, sumWeight, data.Wmin);
}

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
    fscanf(arq, "%d", &data.nItems);
    fscanf(arq, "%d", &data.mBins);
    fscanf(arq, "%d", &data.set);
    fscanf(arq, "%d", &data.numInstance);
    
    //  cost os bins types
    data.ck.clear();
    data.ck.resize(data.mBins);

    // weigth of bins types
    data.Wk.clear();
    data.Wk.resize(data.mBins);

    data.Wmin = 9999999;                          // capacidade do bin com menor capacidade
    data.Wmax = -9999999;                         // capacidade do bin com maior capacidade
    data.Wmean = 0;

    for (int k=0; k<data.mBins; k++)
    {
        fscanf(arq, "%d", &data.Wk[k]);
        fscanf(arq, "%d", &data.ck[k]);

        if (data.Wk[k] < data.Wmin)
            data.Wmin = data.Wk[k];

        if (data.Wk[k] > data.Wmax)
            data.Wmax = data.Wk[k];

        data.Wmean += data.Wk[k];
    }

    data.Wmean = data.Wmean / data.mBins;

    // weigth of items
    data.wi.clear();
    data.wi.resize(data.nItems);

    int sumWeight = 0;

    for (int i=0; i<data.nItems; i++)
    {
        fscanf(arq, "%d", &data.wi[i]);
        sumWeight += data.wi[i];
    }
    fclose(arq);

    // reduzir os itens maiores que completem um bin
    Reduction(data);
    
    // criar os itens
    data.itemsCre.resize(data.nItems);
    data.itemsDec.resize(data.nItems);
    for (int j = 0; j < data.nItems; j++){
        data.itemsCre[j] = j;                                  
        data.itemsDec[j] = j;                                  
    }

    // ordenar em ordem crescente em relacao ao peso
    std::sort(data.itemsCre.begin(), data.itemsCre.end(), [&](int i1, int i2) {
        return data.wi[i1] < data.wi[i2];
    });

    // ordenar em ordem decrescente em relacao ao peso
    std::sort(data.itemsDec.begin(), data.itemsDec.end(), [&](int i1, int i2) {
        return data.wi[i1] > data.wi[i2];
    });

    // imprimir
    // if(0)
    // {
    //     printf("\nItems: %d \nBins types: %d \nUBbins: %d\n", data.nItems, data.mBins, data.UBbins);
    //     for (int k=0; k<data.mBins; k++)
    //     {
    //         printf("%d \t", data.Wk[k]);
    //         printf("%d \t", data.ck[k]);
    //     }
    //     printf("\n");
    //     for (int i=0; i<data.nItems; i++)
    //     {
    //         printf("%d \t", data.wi[i]);
    //     }
    //     printf("\n n = %d", n);
    //     getchar();
    // }
    // if(1)
    // {
    //     // Encontrar o último '/'
    //     char* lastSlash = strrchr(nameTable, '/');
    //     char name2[200]; // Buffer para armazenar o novo nome  
    //     if (lastSlash) {
    //         // Pular o '/'
    //         std::string filename = lastSlash + 1;
    //         // Substituir a extensão "H&S" por "HS"
    //         size_t pos = filename.find("&");
    //         if (pos != std::string::npos) {
    //             filename.replace(pos, 1, "e");  // Substitui "H&S" por "HS"
    //         }
    //         // Substituir a extensão ".txt" por ".dat"
    //         size_t pos2 = filename.find(".txt");
    //         if (pos2 != std::string::npos) {
    //             filename.replace(pos2, 4, ".dat");  // Substitui ".txt" por ".dat"
    //         }
    //         // Copiar o resultado para name2
    //         strncpy(name2, filename.c_str(), sizeof(name2) - 1);
    //         name2[sizeof(name2) - 1] = '\0'; // Garantir terminação nula
    //         printf("\n%s\n", name2);
    //     }       
    //     FILE *arq2;
    //     arq2 = fopen(name2,"w");
    //     if (arq2 == NULL)
    //     {
    //         printf("\nERROR: File (%s) not found!\n",name2);
    //         getchar();
    //         exit(1);
    //     }
    //     // => read data
    //     fprintf(arq2, "param n := %d;   # Número de itens\n", nItems);
    //     fprintf(arq2, "param m := %d;    # Número de bins\n", mBins);      
    //     fprintf(arq2, "\nparam W := \n");
    //     for (int k=0; k<mBins; k++)
    //     {
    //         fprintf(arq2, "%d %d\n", k+1, Wk[k]);
    //     }
    //     fprintf(arq2, ";  # Capacidades dos bins\n");
    //     fprintf(arq2, "\nparam c := \n");
    //     for (int k=0; k<mBins; k++)
    //     {
    //         fprintf(arq2, "%d %d\n", k+1, ck[k]);
    //     }
    //     fprintf(arq2, ";  # custo dos bins\n");
    //     fprintf(arq2, "\nparam w :=\n");
    //     for (int i=0; i<nItems; i++)
    //     {
    //         fprintf(arq2, "%d %d \n", i+1, wi[i]);
    //     }
    //     fprintf(arq2, ";  # Pesos dos itens\n");
    //     fclose(arq2);
    // }


    // executar o FFD para definir o numero de bins necessarios
    std::vector <TBin> bins(data.UBbins);
    for (int j = 0; j < data.UBbins; j++)
    {            
        double aux = std::uniform_real_distribution<double>(0, 1000)(rng);
        int bin = (int)(aux)%data.mBins;

        bins[j].type = bin;    // tipos de bins
    }

    int numBins=0;
    for(int i=0; i<data.nItems; i++)
    {
        int primeiroBin = -1;

        //encontrar o primeiro bin j que possui capacidade para atender o item i
        int j=0;
        while(j < (int)bins.size())
        {
            //verificar a capacidade disponivel do bin j
            if (bins[j].cap + data.wi[data.itemsDec[i]] <= data.Wk[bins[j].type])
            {
                primeiroBin = j;
                break;
            }
            j++;
        }

        //inserir o item i no primeiro bin aberto disponivel
        if (primeiroBin >= 0)
        {
            // acrescentar o custo do tipo do bin para inserir itens nele
            if (bins[primeiroBin].cap == 0)
                numBins++;  

            bins[primeiroBin].cap += data.wi[data.itemsDec[i]];
            bins[primeiroBin].items.push_back(data.itemsDec[i]);
        }

        // criar um novo bin e alocar o item
        else
        {
            // criar um novo bin que possa inserir o item
            for (int k=0; k<data.mBins; k++)
            {
                if (data.wi[data.itemsDec[i]] <= data.Wk[k])
                {
                    TBin aux;
                    aux.cap = data.wi[data.itemsDec[i]];
                    aux.items.push_back(data.itemsDec[i]);
                    aux.type = k;
                    bins.push_back(aux);

                    numBins++;
                }
            }

        }
    }
    
    
    // definir o tamanho do vetor de random-key
    data.n = numBins;
    data.UBbins = numBins;

    // printf("\n\n n = %d\n\n", data.n);
    
    return data;
}

/************************************************************************************
 Method: Decoder minimum-cost VSBPP
 Description: mapping the random-key solution into a problem solution
*************************************************************************************/
typedef struct {
    int from;
    int to;
    float cost;
} Arc;

void print_shortest_path(const std::vector<int>& sequence, std::vector<Arc> arcs, std::vector<int> path, float min_cost, const TProblemData &data) {
    if (path.empty()) {
        std::cout << "No path found from " << path[0] << " to " << path[path.size()-1] << std::endl;
    } else {
        std::cout << "\n\nShortest path: ";
        for (size_t i = 0; i < path.size(); i++) {
            std::cout << path[i] << (i < path.size() - 1 ? " -> " : "\n");
        }
        std::cout << "Minimum cost: " << min_cost << std::endl;
    }

    std::vector <TBin> bins;
    bins.resize(path.size() - 1);

    // bin corrente
    int total_custo  = 0;
    int k = 0;
    std::cout << "\n\nBins: \n" << k << " = ";
    for (int i = 0; i < (int)sequence.size()-1; i++) 
    {
        // adicionar o ith item da sequencia no bin k
        bins[k].items.push_back(sequence[i]);
        bins[k].cap += data.wi[sequence[i]];

        std::cout << sequence[i] <<  ", ";

        // passar para o proximo bin quando chegar no proximo item do path
        if (sequence[i+1] == path[k+1])
        {
            // descobrir o tipo do bin de acordo com o custo do arco path[k], path[k+1]
            int cost = 0;
            for (const auto& arc : arcs) {
                if (arc.from == path[k] && arc.to == path[k+1]) {
                    cost = arc.cost;
                    break;
                }
            }
            for (int j=0; j<data.mBins; j++)
            {
                if (data.ck[j] == cost)
                {
                    bins[k].type = j;
                    break;
                }
            }


            // // descobrir o tipo do bin de acordo com o peso dos itens
            // for (int j=0; j<data.mBins; j++)
            // {
            //     if (data.Wk[j] >= bins[k].cap)
            //     {
            //         bins[k].type = j;
            //         break;
            //     }
            // }



            std::cout << " | cap = " << bins[k].cap;
            std::cout << " [" << data.Wk[bins[k].type] << "  " << data.ck[bins[k].type] << "]";
            std::cout << " type = " << bins[k].type;

            if (bins[k].cap > data.Wk[bins[k].type])
                std::cout << " ********** Erro = " << data.Wk[bins[k].type] - bins[k].cap;

            total_custo += data.ck[bins[k].type];

            k++;
            if (k < (int)path.size() - 1)
                std::cout << "\n" << k << " = ";
        }
    }
    std::cout << "\n\ncusto total = " << total_custo;
}

void print_graph(std::vector<Arc> arcs, const TProblemData &data) {
    std::cout << "Arcs in the Graph:" << std::endl;
    for (const auto& arc : arcs) {
        std::cout << "(" << arc.from << " -> " << arc.to << ") Cost: " << arc.cost << std::endl;
    }
}

float find_min_bin(int total_weight, const TProblemData &data) {
    float min_cost = 999999;
    for (int i = 0; i < data.mBins; i++) {
        if ((total_weight <= data.Wk[i]) && (data.ck[i] <= min_cost)) {
            min_cost = data.ck[i];
        }
    }
    return min_cost;
}

void generate_graph(const std::vector<int>& sequence, std::vector<Arc> &arcs, const TProblemData &data) {
    int total_weight = 0;
    int N = sequence.size()-1;
    // for (int j = 0; j < nItems; j++) {
    for (int j = 0; j < N; j++) {
        // Conecta cada nó ao próximo na sequência com um custo correspondente ao menor bin que pode armazenar o item.
        arcs.push_back({sequence[j], sequence[j + 1], find_min_bin(data.wi[sequence[j]], data)});
        
        // Calcula pesos acumulados para verificar se múltiplos itens podem ser colocados no mesmo bin.
        total_weight = data.wi[sequence[j]];
        // for (int k = j + 1; k < nItems; k++) {
        for (int k = j + 1; k < N; k++) {
            total_weight += data.wi[sequence[k]];
            // Se couberem, adiciona um arco direto do primeiro item até o último item possível com custo correspondente
            // ao melhor bin que acomoda todos esses itens juntos.
            // Para de tentar adicionar itens quando o peso ultrapassa a capacidade do maior bin disponível.
            if (total_weight <= data.Wmax) {
                arcs.push_back({sequence[j], sequence[k+1], find_min_bin(total_weight, data)});
            } else {
                break;
            }
        }
    }
    // print_graph(arcs, data); getchar();
}

/*int dijkstra(const std::vector<int>& sequence, std::vector<Arc> &arcs) {
    int start = sequence[0]; 
    int end   = sequence[nItems]; 
    int num_nodes = sequence.size();

    std::vector<int> dist(num_nodes, 999999);
    std::vector<int> prev(num_nodes, -1);
    std::vector<bool> visited(num_nodes, false);
    dist[start] = 0;
    
    for (int i = 0; i < num_nodes; i++) {
        int u = -1;
        for (int j = 0; j < num_nodes; j++) {
            if (!visited[j] && (u == -1 || dist[j] < dist[u])) {
                u = j;
            }
        }
        
        if (dist[u] == 9999999) break;
        visited[u] = true;
        
        for (const auto& arc : arcs) {
            if (arc.from == u) {
                int v = arc.to;
                int cost = arc.cost;
                if (dist[u] + cost < dist[v]) {
                    dist[v] = dist[u] + cost;
                    prev[v] = u;
                }
            }
        }
    }
    
    if (print)
    {
        std::vector<int> path;
        for (int at = end; at != -1; at = prev[at]) {
            path.insert(path.begin(), at);
        }
        print_shortest_path(sequence, arcs, path, dist[end], data);
    }

    return dist[end];
}*/

/*float dijkstra(const std::vector<Arc>& arcs, const std::vector<int>& sequence) {
    int N = sequence.size();
    int source = sequence[0];
    int target = sequence[N-1];

    // Construir lista de adjacência
    std::vector<std::vector<std::pair<int, float>>> adj(N); // (vizinho, custo)
    for (const Arc& arc : arcs) {
        adj[arc.from].emplace_back(arc.to, arc.cost);
    }

    // Vetor de distâncias
    std::vector<float> dist(N, 99999);
    std::vector<int> prev(N, -1);
    dist[source] = 0;

    // Min-heap (distância, nó)
    using pii = std::pair<int, int>;
    std::priority_queue<pii, std::vector<pii>, std::greater<pii>> pq;
    pq.emplace(0, source);

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        // Se já passou pelo destino, pode parar
        if (u == target) break;

        // Ignorar entradas desatualizadas
        if (d > dist[u]) continue;

        for (auto [v, cost] : adj[u]) {
            if (dist[v] > dist[u] + cost) {
                dist[v] = dist[u] + cost;
                prev[v] = u; // salva o caminho
                pq.emplace(dist[v], v);
            }
        }
    }
    printf("\nReconstruir o caminho..."); getchar();
    if (1)
    {
        // Reconstruir o caminho
        std::vector<int> path;
        for (int at = target; at != -1; at = prev[at]) {
            path.push_back(at);
        }
        std::reverse(path.begin(), path.end());

        //print_shortest_path(sequence, arcs, path, dist[target], data);
        path.clear();
    }

    printf("\n\n%lf", dist[target]); getchar();

    return dist[target] == 99999 ? 99999 : dist[target]; // 99999 se não alcançável
}*/

float dijkstra(const std::vector<Arc>& arcs, const std::vector<int>& sequence, const TProblemData &data) {
    int N = sequence.size();
    int source = sequence[0];
    int target = sequence[N - 1]; // Corrigido para o último nó da sequência

    // Encontrar o maior nó para definir o tamanho da lista de adjacência
    int maxNode = source;
    for (const Arc& arc : arcs) {
        maxNode = std::max(maxNode, std::max(arc.from, arc.to));
    }

    // Construir lista de adjacência com o tamanho correto
    std::vector<std::vector<std::pair<int, float>>> adj(maxNode + 1);
    for (const Arc& arc : arcs) {
        adj[arc.from].emplace_back(arc.to, arc.cost);
    }

    // Vetor de distâncias e predecessores
    std::vector<float> dist(maxNode + 1, 99999);
    std::vector<int> prev(maxNode + 1, -1);
    dist[source] = 0;

    // Min-heap (distância, nó)
    using pii = std::pair<float, int>;
    std::priority_queue<pii, std::vector<pii>, std::greater<pii>> pq;
    pq.emplace(0.0f, source);

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        if (u == target) break;
        if (d > dist[u]) continue;

        for (auto [v, cost] : adj[u]) {
            if (dist[v] > dist[u] + cost) {
                dist[v] = dist[u] + cost;
                prev[v] = u;
                pq.emplace(dist[v], v);
            }
        }
    }

    // if(print)
    // {
    //     // Verificar se o caminho é alcançável
    //     if (dist[target] == 99999) {
    //         printf("\nCaminho não encontrado.");
    //         return 9999999;
    //     }

    //     // Reconstruir o caminho
    //     print_graph(arcs, data);
    //     std::vector<int> path;
    //     for (int at = target; at != -1; at = prev[at]) {
    //         path.push_back(at);
    //     }
    //     std::reverse(path.begin(), path.end());

    //     // Exibir o caminho
    //     printf("\nCaminho mais curto: ");
    //     for (int node : path) {
    //         printf("%d ", node);
    //     }
    //     printf("\nCusto: %.2f", dist[target]);
    //     getchar();
    // }

    if (dist[target] == 99999) {
        return 9999999;
    }

    return dist[target];
}

/************************************************************************************
 Method: Decoder Bins Type
 Description: mapping the random-key solution into a problem solution
*************************************************************************************/
double Decoder(TSol &s, const TProblemData &data)
{   
    // custo dos bins utilizados
    int custo = 0;

    // create a solution of the problem
    std::vector <TBin> bins(data.UBbins);
    std::vector <int> items(data.nItems);
    
    // transformar o vetor de rk em uma solucao do vsbpp
    for (int j = 0; j < data.UBbins; j++)
    {            
        bins[j].type = floor(s.rk[j] * data.mBins);    // tipos de bins
    }

    items = data.itemsDec;

    // ** empacotar os items nos bins

    // usar a insercao first bin
    // if (s.rk[n-1] < 0.5)
    // {
        for(int i=0; i<data.nItems; i++)
        {
            int primeiroBin = -1;

            //encontrar o primeiro bin j que possui capacidade para atender o item i
            int j=0;
            // while(j < UBbins)
            while(j < (int)bins.size())
            {
                //verificar a capacidade disponivel do bin j
                if (bins[j].cap + data.wi[items[i]] <= data.Wk[bins[j].type])
                {
                    primeiroBin = j;
                    break;
                }
                j++;
            }

            //inserir o item i no primeiro bin aberto disponivel
            if (primeiroBin >= 0)
            {
                // acrescentar o custo do tipo do bin para inserir itens nele
                if (bins[primeiroBin].cap == 0)
                    custo += data.ck[bins[primeiroBin].type];  

                bins[primeiroBin].cap += data.wi[items[i]];
                bins[primeiroBin].items.push_back(items[i]);
            }

            // criar um novo bin e alocar o item
            else
            {
                // criar um novo bin que possa inserir o item
                for (int k=0; k<data.mBins; k++)
                {
                    if (data.wi[items[i]] <= data.Wk[k])
                    {
                        TBin aux;
                        aux.cap = data.wi[items[i]];
                        aux.items.push_back(items[i]);
                        aux.type = k;
                        bins.push_back(aux);

                        custo += data.ck[k];  
                        break;
                    }
                }
            }
        }

    // usar a insercao best bin
    /*else
    {
        for(int i=0; i<nItems; i++)
        {
            int melhorBin = -1;
            float melhorTaxa = 0.0;

            //encontrar o melhor bin j que possui capacidade para atender o item i, aquele que ficara mais cheio
            int j=0;
            while(j < UBbins)
            {
                //verificar a capacidade disponivel do bin j
                if (bins[j].cap + wi[items[i]] <= Wk[bins[j].type])
                {
                    // verificar se eh a melhor capacidade ocupada
                    if ((float)(bins[j].cap + wi[items[i]])/Wk[bins[j].type] > melhorTaxa)
                    {
                        melhorBin = j;
                        melhorTaxa = (float)(bins[j].cap + wi[items[i]])/Wk[bins[j].type];
                            if (melhorTaxa == 1)
                                break;
                    }
                }
                j++;
            }

            //inserir o item i no melhor bin disponivel
            if (melhorBin >= 0)
            {
                // acrescentar o custo do tipo do bin para inserir itens nele
                if (bins[melhorBin].cap == 0)
                    custo += ck[bins[melhorBin].type];  

                bins[melhorBin].cap += wi[items[i]];
                bins[melhorBin].items.push_back(items[i]);
            }

            // penalizar a solucao se um item nao foi empacotado
            else
            {
                penalty += 10000000; 
            }
        }
    }*/

    // ajustar o tipo de cada bin usado
    for(int j=0; j<(int)bins.size(); j++)
    {
        // verificar se eh possivel realizar um ajuste nos tipos de bins
        if (bins[j].cap > 0){
            int melhorCusto = 99999999;
            int melhorK = 0;

            // buscar um custo menor para o bin j
            for (int k=0; k<data.mBins; k++){
                if (bins[j].cap <= data.Wk[k] && data.ck[k] < melhorCusto){
                    melhorCusto = data.ck[k];
                    melhorK = k;
                }
            }

            // atualizar o tipo do bin j se houver melhora no custo
            if (melhorCusto < data.ck[bins[j].type]){
                //atualizar o custo
                custo = custo + data.ck[melhorK] - data.ck[bins[j].type];

                // atualizar o tipo do bin
                bins[j].type = melhorK;

                // gerar uma chave no intervalo do tipo j para atualizar a solucao
                // float fator = 1.0/mBins;
                // double min = melhorK*fator;
                // double max = (melhorK+1)*fator;
                // s.rk[j] = ((double)(rand()%10000)/10000.0)*(max-min)+min;
            }
        }
    }

    // calculate fitness
    s.ofv = custo + data.custoReducao;

    // if (print)
    // {
    //     std::vector<int> verifica(data.nItems,0);
    //     printf("\nSolution VSBPP: \n");
    //     for (int i=0; i<(int)bins.size(); i++)
    //     {
    //         printf("\nBin[%d] [%d: - %d %d] [%d] = ", i, bins[i].type, data.Wk[bins[i].type], data.ck[bins[i].type], bins[i].cap);
    //         for (int j = 0; j < (int)bins[i].items.size(); j++)
    //         {
    //             printf(" %d ", bins[i].items[j]);

    //             verifica[bins[i].items[j]] = 1;

    //             if (bins[i].cap > data.Wk[bins[i].type])
    //             {
    //                 printf("\n\nERRO Cap\n\n");
    //                 getchar();
    //             }
    //         }
    //     }
    //     for (int i = 0; i < data.nItems; i++)
    //     {
    //         if (verifica[i] == 0)
    //         {
    //             printf("\n\nERRO item nao alocado %d\n\n", i);
    //             getchar();
    //         }
    //     }
        
    //     printf("\n\nSolution VSBPP: %.0lf\n\n", s.ofv);
    // }

    return s.ofv;
}

/************************************************************************************
 Method: Decoder (First fit )
 Description: mapping the random-key solution into a problem solution
*************************************************************************************/
/*double DecoderFirst(TSol &s)
{
    // custo dos bins utilizados
    int custo = 0;

    // create a solution of the problem
    std::vector <TBin> bins;
    std::vector <int> items;

    items.resize(nItems);
    bins.resize(UBbins);
    
    // transformar o vetor de rk em uma solucao do vsbpp
    for (int j = 0; j < n; j++)
    {
        if (j < nItems){
            items[j] = j;                                    // parte 1 (items)
        }
        else{              
            bins[j-nItems].type = floor(s.rk[j] * mBins);    // parte 2 (tipos de bins)
        }
    }

    // ordenar as n primeiras chaves (itens)
    std::sort(items.begin(), items.end(), [&s](int i1, int i2) {
        return s.rk[i1] < s.rk[i2];
    });

    //empacotar os items nos bins
    int penalty = 0;
    for(int i=0; i<nItems; i++)
    {
        int primeiroBin = -1;

        //encontrar o primeiro bin j que possui capacidade para atender o item i
        int j=0;
        while(j < UBbins)
        {
            //verificar a capacidade disponivel do bin j
            if (bins[j].cap + wi[items[i]] <= Wk[bins[j].type])
            {
                primeiroBin = j;
                break;
            }
            j++;
        }

        //inserir o item i no primeiro bin aberto disponivel
        if (primeiroBin >= 0)
        {
            // acrescentar o custo do tipo do bin para inserir itens nele
            if (bins[primeiroBin].cap == 0)
                custo += ck[bins[primeiroBin].type];  

            bins[primeiroBin].cap += wi[items[i]];
            bins[primeiroBin].items.push_back(items[i]);
        }

        // penalizar a solucao se um item nao foi empacotado
        else
        {
            penalty += 10000000; 
        }
    }

    // calculate fitness
    s.ofv = custo + penalty;
    return s.ofv;
}*/

/************************************************************************************
 Method: Decoders (Best fit)
 Description: mapping the random-key solution into a problem solution
*************************************************************************************/
/*double DecoderBest(TSol &s)
{
    // custo dos bins utilizados
    int custo = 0;

    // create a solution of the problem
    std::vector <TBin> bins;
    std::vector <int> items;

    items.resize(nItems);
    bins.resize(UBbins);
    
    // transformar o vetor de rk em uma solucao do vsbpp
    for (int j = 0; j < n; j++)
    {
        if (j < nItems){
            items[j] = j;                                    // parte 1 (items)
        }
        else{              
            bins[j-nItems].type = floor(s.rk[j] * mBins);    // parte 2 (tipos de bins)
        }
    }

    // ordenar as n primeiras chaves (itens)
    std::sort(items.begin(), items.end(), [&s](int i1, int i2) {
        return s.rk[i1] < s.rk[i2];
    });

    //empacotar os items nos bins
    int penalty = 0;
    for(int i=0; i<nItems; i++)
    {
        int melhorBin = -1;
        float melhorCapOcupada = 0.0;

        //encontrar o melhor bin j que possui capacidade para atender o item i, aquele que ficara mais cheio
        int j=0;
        while(j < UBbins)
        {
            //verificar a capacidade disponivel do bin j e se eh a melhor capacidade ocupada
            if (bins[j].cap + wi[items[i]] <= Wk[bins[j].type])
            {
                // printf("\nj: %d, ", j);
                // printf("capOcup: %lf, ", (float)(bins[j].cap + wi[items[i]])/Wk[bins[j].type]);
                // printf("Melhor: %lf, ", melhorCapOcupada);
                // printf("Bin: %d, ", melhorBin);

                if ((float)(bins[j].cap + wi[items[i]])/Wk[bins[j].type] > melhorCapOcupada)
                {
                    melhorBin = j;
                    melhorCapOcupada = (float)(bins[j].cap + wi[items[i]])/Wk[bins[j].type];
                    // printf(" >> Entrei");
                }
            }
            j++;
        }

        // printf("\tMelhor bin: %d, ", melhorBin); getchar();
        //inserir o item i no melhor bin disponivel
        if (melhorBin >= 0)
        {
            // acrescentar o custo do tipo do bin para inserir itens nele
            if (bins[melhorBin].cap == 0)
                custo += ck[bins[melhorBin].type];  

            bins[melhorBin].cap += wi[items[i]];
            bins[melhorBin].items.push_back(items[i]);
        }

        // penalizar a solucao se um item nao foi empacotado
        else
        {
            penalty += 10000000; 
        }
    }

    // criterio de desempate
    // double menorTaxa = 999999;
    // for(int j=0; j<mBins; j++)
    // {
    //     if (bins[j].cap > 0)
    //     {
    //         if ((double)bins[j].cap/Wk[bins[j].type] < menorTaxa)
    //             menorTaxa = (double)bins[j].cap/Wk[bins[j].type];
    //     }
    // }

    // ajustar o tipo de bin
    for(int j=0; j<mBins; j++)
    {
        // verificar se eh possivel realizar um ajuste nos tipos de bins
        if (bins[j].cap > 0){
            int melhorCusto = 99999999;
            int melhorK = 0;

            // buscar um custo menor para o bin j
            for (int k=0; k<mBins; k++){
                if (bins[j].cap <= Wk[k] && ck[k] < melhorCusto){
                    melhorCusto = ck[k];
                    melhorK = k;
                }
            }

            // atualizar o tipo do bin j se houver melhora no custo
            if (melhorCusto < ck[bins[j].type]){
                // printf("\nmelhor custo: %d, melhorK: %d, custo atual: %d, tipo atual: %d", melhorCusto, melhorK, ck[bins[j].type], bins[j].type);
                
                //atualizar o custo
                custo = custo + ck[melhorK] - ck[bins[j].type];

                // atualizar o tipo do bin
                bins[j].type = melhorK;

                // gerar uma chave no intervalo do tipo j para atualizar a solucao
                // printf(", rk atual: %lf, rk nova: ", s.rk[nItems+j]);
                
                float fator = 1.0/mBins;
                double min = melhorK*fator;
                double max = (melhorK+1)*fator;
                s.rk[nItems+j] = ((double)(rand()%10000)/10000.0)*(max-min)+min;
                
                // printf("%lf", s.rk[nItems+j]);
                // getchar();
            }
        }
    }

    // calculate fitness
    s.ofv = custo + penalty; //  + menorTaxa;
    return s.ofv;
}*/

/*double DecoderGrafo(TSol &s)
{
    // custo dos bins utilizados
    double custo = 0.0;

    // create a solution of the problem
    std::vector <int> sequence(nItems);
    std::vector<Arc> arcs;

    // transformar o vetor de rk em uma solucao do vsbpp
    for (int j = 0; j < nItems; j++){
        sequence[j] = j;
    }

    // ordenar as n chaves (itens)
    std::sort(sequence.begin(), sequence.end(), [&s](int i1, int i2) {
        return s.rk[i1] < s.rk[i2];
    });

    // ordenar as nItems * s.rk[n-1] chaves (itens) de acordo com o peso dos itens da sequencia
    // int cutoff = static_cast<int>(nItems * 0.25);
    // std::sort(sequence.begin(), sequence.begin() + cutoff, [&](int i1, int i2) {
    //     return wi[i1] > wi[i2];
    // });

    // dummy node n + 1 
    sequence.push_back(nItems);
    
    // generate_graph
    generate_graph(sequence, arcs, data);

    // find the mininum cost of the VSBPP
    // custo = dijkstra(sequence, arcs, data);
    custo = dijkstra(arcs, sequence, data);
    
    arcs.clear();
    sequence.clear();

    // return fitness
    return custo;
}*/

/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem(TProblemData &data){
    data.ck.clear();
    data.Wk.clear();
    data.wi.clear();
}

#endif