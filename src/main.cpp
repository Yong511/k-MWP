#include "Graph.h"
#include "Utility.h"
#include "popl.hpp"
#pragma GCC optimize(3,"Ofast","inline")

using namespace std;
using namespace popl;

void print_usage() {
    printf("Example usage: ./MWPlex -g path_to_graph -k 3 -a MWPA -w 0\n");
}

int main(int argc, char *argv[]) {

    bool output = false;
    bool binary_input = false;
    bool mode = false;
    string alg;
    int weight_choice = -1;

    OptionParser op("Allowed options");
    auto help_option = op.add<Switch>("h", "help", "\'produce help message\'");
    auto graph_option = op.add<Value<string>>("g", "graph", "\'path to input graph file\'");
    auto alg_option = op.add<Value<string>>("a", "alg", "'algorithm name' (Base | MWPA- | MWPA)", "MWPA", &alg);
    auto k_option = op.add<Value<int>>("k", "k", "\'the value of k for k-plex\'");
    auto weight_option = op.add<Value<int>>("w", "weight", "'weight distribution id (0:degree,1:normal,2:exponential,3:powerlaw,4:lognormal,5:pareto,6:gamma), default uses weight.txt'", -1, &weight_choice);


    op.parse(argc, argv);
    if(help_option->is_set()||argc <= 1) {
        cout << op << endl;
        if(argc <= 1) {
            print_usage();
            return 0;
        }
    }
    if(!graph_option->is_set()) {
        printf("!!! The argument -g is required! Exit !!!\n");
        return 0;
    }
    if(!k_option->is_set()) {
        printf("!!! The argument -k is required! Exit !!!\n");
        return 0;
    }
    if(weight_option->is_set() && (weight_choice < 0 || weight_choice > 6)) {
        printf("!!! The argument -w must be in [0,6]! Exit !!!\n");
        return 0;
    }
    Graph *graph = new Graph(graph_option->value().c_str(),k_option->value(), 2 * k_option->value() - 1, weight_choice);
    if(weight_choice >= 0) {
        printf("Example usage: ./MWPlex -g %s -k  %d -a %s -w %d -q %d\n", graph_option->value().c_str(), k_option->value(), alg.c_str(), weight_choice, 2 * k_option->value() - 1);
    } else {
        printf("Example usage: ./MWPlex -g %s -k  %d -a %s -q %d\n", graph_option->value().c_str(), k_option->value(), alg.c_str(), 2 * k_option->value() - 1);
    }
    graph->read_graph();

    if(strcmp(alg.c_str(), "Base") == 0) graph->exact_base();
    if(strcmp(alg.c_str(), "MWPA-") == 0) graph->exact_lower();
    if(strcmp(alg.c_str(), "MWPA") == 0) graph->exact(0);
    if(strcmp(alg.c_str(), "MWPA-GR") == 0) graph->exact(1);
    if(strcmp(alg.c_str(), "MWPA-s012") == 0) graph->exact(2);
    if(strcmp(alg.c_str(), "MWPA-b1234") == 0) graph->exact(3);

    
    return 0;
}
