#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <utility>
#include "../core/graph.hpp"
#include "../core/circuit.hpp"

#ifndef OPTIMIZE_HPP
#define OPTIMIZE_HPP
//Cost Metric
//RCNOT_COST = DISTANCE * 4
//COUPLING_COST = DISTANCE * EDGE_WEIGHT WHERE EDGE_WEIGHT = BASE ^ (MAX_CNOT_DEPTH - CURRENT_GATE_DEPTH)
enum class Metric_Type {RCNOT_COST, COUPLING_COST};
//A1: Global ordering & RCNOT gates
//A2: Global ordering & SWAP gates
//A3: Global ordering & RCNOT + SWAP gates
enum class Approach_Type {A1, A2, A3};  
typedef struct {
  grph::node_t control;
  grph::node_t target;
  //unsigned window;
  double swap_cost;
  double swaps;
}mincost_t;

void displayPath(std::vector<grph::node_t>& path);
bool isInvalidPath(std::vector<grph::node_t> & path, std::map<rev::line_t, rev::line_t> & plmap);
void removeDuplicatePath(std::vector<std::vector<grph::node_t> > & paths);
void addInitNode(std::vector<std::vector<grph::node_t> > & paths, grph::node_t nd);
void displayAllPaths(std::vector<std::vector<grph::node_t> > & paths);
void getAllPaths(grph::node_t p1, grph::node_t p2, std::vector<std::vector<grph::node_t> > & paths, grph::Graph & pg, int i);

int getAllPaths(grph::node_t p1, grph::node_t p2, 
      std::vector<std::vector<grph::node_t> > & paths, std::pair<int, int> dim, int i); //hexagonal architecture
void getPath(grph::node_t p1, grph::node_t p2, std::vector<grph::node_t> & path, grph::Graph & pg, int i);

mincost_t minCostPathIndex(std::vector<grph::node_t>& path, std::map<grph::node_t, grph::node_t> & temp_phmap, grph::Graph & pg, grph::Graph & temp_g);

void updateMapping(std::vector<grph::node_t>& path, std::map<grph::node_t, grph::node_t> & temp_phmap, grph::Graph& pg, 
                                      std::vector<grph::node_t> & nodes, grph::node_t control, grph::node_t target);
void updateMapping(std::vector<grph::node_t>& path, std::map<grph::node_t, grph::node_t> & phmap, std::pair<int, int> dim, 
                             std::vector<grph::node_t> & nodes, grph::node_t control, grph::node_t target); //hexagonal 

//minimal estimation of coupling cost exploring all possible paths
double swapCost(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& cg, int base, unsigned depth);
double swapCost(std::map<grph::node_t, grph::node_t>& phmap, std::pair<int, int> dim, grph::Graph& cg, int base, unsigned depth); //hexagonal

//local ordering
//estimation using coupling cost metric
double mappingCost(std::map<grph::node_t, grph::node_t>& phmap, std::map<rev::line_t, rev::line_t>& lgmap, grph::Graph& pg, grph::Graph& cg); //IBM QX
double mappingCost(std::map<grph::node_t, grph::node_t>& phmap, std::map<rev::line_t, rev::line_t>& lgmap,  std::pair<int, int> dim, grph::Graph& cg); //hexagonal

//global ordering
//estimation using coupling cost metric
double mappingCost(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& cg); //ibm qx
double mappingCost(std::map<grph::node_t, grph::node_t>& phmap, std::pair<int, int> dim, grph::Graph& cg); //hexagonal

//local ordering
//estimation using remote CNOT cost metric
double remoteCNOTCost(std::map<grph::node_t, grph::node_t>& phmap, std::map<rev::line_t, rev::line_t>& lgmap, grph::Graph& pg, grph::Graph& cg); //IBM QX
double remoteCNOTCost(std::map<grph::node_t, grph::node_t>& phmap, std::map<rev::line_t, rev::line_t>& lgmap,  std::pair<int, int> dim, grph::Graph& cg); //hexagonal

//global ordering
//estimation using remote CNOT cost metric
double remoteCNOTCost(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& cg); //IBM QX
double remoteCNOTCost(std::map<grph::node_t, grph::node_t>& phmap, std::pair<int, int> dim, grph::Graph& cg); //hexagonal

//void optimize(rev::Circuit & ckt); //Need to redefine due angle type changed from double to string

double computeDepth(rev::Circuit & ckt);
long computeCNOTGates(rev::Circuit & ckt);
double move_cost(pair<grph::node_t, grph::node_t>& swapped_nodes, grph::Graph& p); 
void generateRemoteCNOT(std::vector<rev::Gate>& gates, rev::Gate & g, std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg);

#endif /* OPTIMIZE_HPP */
