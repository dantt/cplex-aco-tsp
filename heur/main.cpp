#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cfloat>
#include <cmath>
#include <random>
#include <chrono>
#include <cstdlib>

using namespace std;

class Graph {
    
public: 
    Graph(int n_nodes);
    void add_edge(int, int, double);
    void update_trail(int, int, double);
    double get_trail(int, int);
    double get_edge(int, int);
    double get_metric(int, int);
    void print();
    void calculate_metric(int);
    void update_metric(int, int, double);
    
private:
    int size;
    vector<vector<double>*>* edges;
    vector<vector<double>*>* trails;
    vector<vector<double>*>* metric;
};

Graph::Graph(int n_nodes) {
    this->size = n_nodes;
    this->edges = new vector<vector<double>*>();
    this->trails = new vector<vector<double>*>();
    this->metric = new vector<vector<double>*>();
    for (int i=0; i<n_nodes; i++) {
        this->edges->push_back(new vector<double>(n_nodes, 0.0));
        this->trails->push_back(new vector<double>(n_nodes, 1));
        this->metric->push_back(new vector<double>(n_nodes, 0.0));
    }
    
}

void Graph::add_edge(int from, int to, double weight) {
    (*(*this->edges)[from])[to] = weight;
    return;
}

void Graph::update_trail(int from, int to, double amount) {
    (*(*this->trails)[from])[to] = amount;
    (*(*this->trails)[to])[from] = amount;
}

double Graph::get_trail(int from, int to) {
    return (*(*this->trails)[from])[to];
}

double Graph::get_edge(int from, int to) {
    return (*(*this->edges)[from])[to];
}

double Graph::get_metric(int from, int to) {
    return (*(*this->metric)[from])[to];
}

void Graph::print() {
    for (int i=0; i<this->edges->size(); i++) {
        for (int j=0; j<this->edges->size(); j++) {
            cout << (*(*this->edges)[i])[j] << " ";
        }
        cout << endl;
    }
};

void Graph::calculate_metric(int beta) {
    for (int i=0; i<this->size; i++) {
        for (int j=i; j<this->size; j++) {
            this->update_metric(i, j, this->get_trail(i, j) * pow( pow(this->get_edge(i, j), -1), beta));
        }
    }
}

void Graph::update_metric(int from, int to, double amount) {
    (*(*this->metric)[from])[to] = amount;
    (*(*this->metric)[to])[from] = amount;
}

class ant {
public:
    ant(int, int);
    void print_path();
    void reset();
    vector<int> visited = vector<int>();
    int current = 0;
    vector<int> tour = vector<int>();
    double cost = 0;
    int start = 0;
    int N = 0;
};

ant::ant(int N, int start = 0) {
    this->start = start;
    for (int i = 0; i<N; i++) {
        this->visited.push_back(0);
        this->tour.push_back(0);
    }
    this->visited[start] = 1;
    this->tour[0] = start;
    this->current = start;
    this->N = N;
}
void ant::print_path() {
    for (int i=0; i<this->tour.size(); i++) {
        cout << this->tour[i] << "->";
    }
    cout << this->start << endl;
}
void ant::reset() {
    this->cost = 0;
    for (int i=0; i<N; i++) {
        this->visited[i] = 0;
        this->tour[i] = 0;
    }
    this->visited[start] = 1;
    this->tour[0] = start;
    this->current = this->start;
}
    

int main (int argc, char const *argv[])
{
    //TODO pls check argv[1] exists or die
    string problem_file = argv[1];
    ifstream data(problem_file.c_str());
    string line;
    
    int N = 0;
    
    getline(data, line);
    stringstream lineStream(line);
    string cell;
    while(getline(lineStream,cell,','))
    {
        N++;
    }
    
    Graph* graph = new Graph(N);
    
    ifstream data2(problem_file.c_str());
    
    // Leggi i pesi, crea il grafo
    int node_n = 0;
    while(getline(data2,line))
    {
        stringstream lineStream(line);
        string cell;
        int neighbor_n = 0;
        while(getline(lineStream,cell,','))
        {
            double weight = (atof(cell.c_str()) == 0.0) ? FLT_MAX : atof(cell.c_str());
            graph->add_edge(node_n, neighbor_n, weight);
            neighbor_n++;
        }
        node_n++;
    }
    
    // Params    
    float alpha = 0.1;
    int beta = 2;
    float q0 = 0.9;
    float t0 = pow(N, -1);

    
    // Formiche
    int n_ants = (N >= 36) ? 10 : 5;
    vector<ant*>* ants = new vector<ant*>();
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 generator (seed);
    uniform_int_distribution<int> rndint(0,N-1);
    uniform_real_distribution<double> rnd01(0.0,1.0);
    for (int i=0; i<n_ants; i++) {
        ants->push_back(new ant(N, rndint(generator)));
    }
    
    
    // Tour 
    vector<int> best_tour(N,0);
    double best_val = FLT_MAX;
    graph->calculate_metric(beta);
    int g_its = 0;
    int max_its = 10;
    bool stopping_criterion = false;
    
    while (g_its < max_its && !stopping_criterion) {
        int iteration = 1;
        vector<double> costs(N,0);
        for (int i=0; i<n_ants; i++) {
            ((*ants)[i])->reset();
        }
        while (iteration < N) {
            for (int i=0; i<ants->size(); i++) {
                ant* ant = (*ants)[i];
                int maxindex = 0;
                double max = 0.0;
                double total = 0.0;
                for (int k=0; k<N; k++) {
                    if (ant->visited[k] != 1) {
                        double metric = graph->get_metric(ant->current, k);
                        total += metric;
                        if (metric > max) {
                            max = metric;
                            maxindex = k;
                        }
                    }
                }
                int next;
                if ( rnd01(generator) <= q0 ) {
                    next = maxindex;
                }
                else {
                    int i = 0;
                    double sums = 0.0;
                    double rval = rnd01(generator);
                    do {
                        if (ant->visited[i] != 1) {
                            sums += graph->get_metric(ant->current, i);
                        }
                        i++;
                    }
                    while (rval - sums / total > 0 );
                    next = i - 1;
                }
                
                graph->update_trail(ant->current, next, (1-alpha)*(graph->get_trail(ant->current, next) + alpha*t0));
                
                ant->visited[next] = 1;
                ant->tour[iteration] = next;
                ant->cost += graph->get_edge(ant->current, next);
                costs[i] = ant->cost;
                ant->current = next;
            }
            iteration++;
        }
        
        int ant_i;
        double best_it_val = FLT_MAX;
        for (int i=0; i<ants->size(); i++) {
            ant* ant = (*ants)[i];
            ant->cost += graph->get_edge(ant->current, ant->start);
            if (ant->cost < best_it_val) {
                best_it_val = ant->cost;
                ant_i = i;
            }
        }
        if (best_it_val < best_val) {
            best_val = best_it_val;
            best_tour = vector<int>(((*ants)[ant_i])->tour);
        }
        
        for (int i=0; i<N-1; i++) {
            graph->update_trail(best_tour[i], best_tour[i+1], (1-alpha)*(graph->get_trail(best_tour[i], best_tour[i+1])) + alpha*pow(best_val, -1));
        }
        graph->update_trail(best_tour[N], best_tour[0], (1-alpha)*(graph->get_trail(best_tour[N], best_tour[0])) + alpha*pow(best_val, -1));
        g_its++;
    }
    
    string results_file = "../results/heur/";
    string filename(argv[1]);
    results_file =  results_file + filename.substr(filename.rfind('/')+1);
    ofstream res (results_file.c_str());
    
    if (res.is_open())
    {
        res << fixed << setprecision(8) << best_val;
        res.close();
    }
    
    
    
    return 0;
}
