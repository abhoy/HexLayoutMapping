#include <iostream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include "optimize.hpp"
#include "../core/graph.hpp"
#include "../core/circuit.hpp"

#ifndef LOCAL_HPP
#define LOCAL_HPP

extern rev::window_t window;
/*typedef struct{
  rev::line_t control;
  rev::line_t target;
  //unsigned window;
  double swap_cost;
  double swaps;
}mincost_t;*/



void insertSWAP(std::vector<rev::Gate>& swap_gates, std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& lgmap,
                             std::map<rev::line_t, rev::line_t>& glmap, std::map<rev::line_t, rev::line_t>& plmap,
                             rev::line_t control, rev::line_t target);  

//void insertRCNOT(std::vector<rev::Gate>& rcnot_gates, rev::Gate & g, std::vector<grph::node_t>& path, 
//                             rev::line_t control, rev::line_t target);
void insertRCNOT(std::vector<rev::Gate>& rcnot_gates, std::vector<grph::node_t>& path, 
                               rev::line_t control, rev::line_t target);

//Hexagonal Approach 2
mincost_t minCostPathIndex(std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& lgmap,
                             std::map<rev::line_t, rev::line_t>& glmap, std::map<rev::line_t, rev::line_t>& plmap, 
                             std::map<grph::node_t, grph::node_t> & phmap, std::vector<rev::Gate>& gates,
                             std::vector<rev::line_t>& lines, std::pair<int, int> dim, grph::Graph & lg, 
                             unsigned pos, int base, Metric_Type metric);

//Hexagonal Approach 3
mincost_t minCostPathIndexv2(std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& lgmap,
                             std::map<rev::line_t, rev::line_t>& glmap, std::map<rev::line_t, rev::line_t>& plmap, 
                             std::map<grph::node_t, grph::node_t> & phmap, std::vector<rev::Gate>& gates,
                             std::vector<rev::line_t>& lines, std::pair<int, int> dim, grph::Graph & lg, 
                             unsigned pos, int base, Metric_Type metric);

//Hexagonal Mapping Approach 1 : Metric = Coupling Cost, Global Ordering, Remote CNOT gate
void localorder(rev::Circuit & ckt, std::map<grph::node_t, grph::node_t> & phmap, std::pair<int, int> dim);

//Hexagonal Mapping Approach 2 : Metric =Coupling Cost, Global Ordering, Local Ordering, SWAP gate
void localorder(rev::Circuit & ckt, std::map<grph::node_t, grph::node_t> & phmap, std::pair<int, int> dim, 
                                                           grph::Graph & lg, int base, Metric_Type metric);
//Hexagonal Mapping Approach 3 : Metric =Coupling Cost, Global Ordering, Local Ordering, SWAP + RCNOT gate
void localorderv2(rev::Circuit & ckt, std::map<grph::node_t, grph::node_t> & phmap, std::pair<int, int> dim, 
                                                             grph::Graph & lg, int base, Metric_Type metric);

#endif /* LOCAL_HPP */
