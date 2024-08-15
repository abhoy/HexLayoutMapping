//Mapping qubits in a hexagonal architecture
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <cstdlib>
#include "core/graph.hpp"
#include "core/circuit.hpp"
#include "lib/local.hpp"
#include "lib/optimize.hpp"
#include "lib/global.hpp"
#define AVG 1
#define BASE 2   //edge weight of IG is power(BASE, depth - pos);
rev::window_t window = 0;
int main( int argc, char ** argv ) 
{ 
  std::string circuit;
  std::ofstream ofs1, ofs2, ofs3;
  std::pair<int, int> dim;
  int width, height;
  Metric_Type metric;  //RCNOT = 0 and Coupling Cost = 1 (DEFAULT = 0)
  //Approach_Type approach;//RCNOT = 0, SWAP = 1 and SWAP + RCNOT = 2 (DEFAULT = 0)
  if (argc != 6)
  {
    std::cout<<"Usage:\n\t 1. circuitname[*.qasm] \n\t 2. outputpath \n\t 3. width (hexagonal layout) \n\t 4. height (hexagonal layout) \n\t 5. metric (RCNOT = 0 and Coupling Cost = 1)"<<std::endl; //\n\t 6. approach (RCNOT = 0, SWAP = 1 and SWAP + RCNOT = 2)"
    std::exit(EXIT_FAILURE);
  }
  else{
    circuit = argv[1]; 
    string path = argv[2];     
    width = atoi(argv[3]);
    height = atoi(argv[4]);
    dim = std::pair<int, int> {width, height};
    metric = atoi(argv[5]) == 1 ?  Metric_Type::COUPLING_COST :  Metric_Type::RCNOT_COST;

    std::string path2 = path + "/M" + argv[5] + "A1/" + circuit.substr(circuit.find_last_of("/\\") + 1); //Approach 2
    std::string path3 = path + "/M" + argv[5] + "A2/" + circuit.substr(circuit.find_last_of("/\\") + 1); //Approach 3
    
    ofs2.open(path2.c_str());
    ofs3.open(path3.c_str());
  }
   
  rev::Circuit ckt = rev::Circuit(); 
  ckt.readQASM(circuit);
  
  window  = ckt.getCNOTDepth(); 
  unsigned CNOTs = ckt.twoQubitGateCount();
  
  if ( BASE == 2 && window > 10) window = 20;
  
  clock_t start1, end1, start2, end2, start3, end3, start4, end4;
  
  std::map<grph::node_t, grph::node_t> phmap_min = std::map<grph::node_t, grph::node_t> ( );
  
  double A1_time = 0., A2_time = 0., A3_time = 0.;
  double average_map_time = 0.,  min_A1_time = 0., min_A2_time = 0., min_A3_time = 0.;
  double average_bound = 0.;
  
  rev::window_t minw_a2, minw_a3;
  rev::Circuit min_ckt2, min_ckt3;
  
  start1 = clock(); 
  grph::Graph lg;
  if (metric == Metric_Type::COUPLING_COST)
    lg = grph::Graph(ckt, grph::NEW, window, BASE);
  else 
    lg = grph::Graph(ckt);
  //lg.displayGraph();
  std::map<grph::node_t, grph::node_t> phmap = std::map<grph::node_t, grph::node_t> ( );
  lg.mapLayout(phmap, std::pair<int, int> {width, height});
     
  ga_search(phmap, dim, lg, metric); //Global Ordering, Cost Metric: RCNOT, Coupling Cost
  
  end1 = clock();

  average_map_time += ( (end1 - start1) / (double) CLOCKS_PER_SEC );  
  average_bound += 2 * mappingCost(phmap, dim, lg);

  //with variable lookahead parameter
  for (auto i = 0; i <50; i+=10){
    window = 10 + i;
    rev::Circuit ckt2 = rev::Circuit(); //for approach 2
    ckt2.readQASM(circuit);
  
    phmap_min.clear();
    for(const auto &  l : phmap) {
        phmap_min[l.first] = l.second;
        //std::cout<<"logical qubit: "<<l.first<<" \tPhysical qubit: "<<l.second<<std::endl;
    }
  
    start3 = clock(); 
    //Approach 2 for both RCNOT cost and Coupling cost
    //Approach 2: Global Ordering, Local Ordering, SWAP gate
    localorder(ckt2, phmap_min, dim, lg, BASE, metric); 
    end3 = clock();

    A2_time = ( (end3 - start3) / (double) CLOCKS_PER_SEC ); 
    if (i == 0 || ckt2.getGates().size() < min_ckt2.getGates().size()){
      min_ckt2 = ckt2;
      min_A2_time = A2_time;
      minw_a2 = window;
    }
    rev::Circuit ckt3 = rev::Circuit(); //for approach 3
    ckt3.readQASM(circuit);

    phmap_min.clear();
    for(const auto &  l : phmap) {
      phmap_min[l.first] = l.second;
      //std::cout<<"logical qubit: "<<l.first<<" \tPhysical qubit: "<<l.second<<std::endl;
    }

    start4 = clock();
    //Approach 3 for both RCNOT cost and Coupling cost
    //Approach 3: Global Ordering, Local Ordering, SWAP + RCNOT gate      
    localorderv2(ckt3, phmap_min, dim, lg, BASE, metric);

    end4 = clock();  

    A3_time = ( (end4 - start4) / (double) CLOCKS_PER_SEC );
    if (i == 0 || ckt3.getGates().size() < min_ckt3.getGates().size()){
      min_ckt3 = ckt3;
      min_A3_time = A3_time;
      minw_a3 = window;
    }
    ckt2.clear();
    ckt3.clear();    
  }
  //average_map_time = average_map_time / 10;
  //average_bound = average_bound / 10;
  std::cout<<"#Circuit:\t"<<circuit.substr(circuit.find_last_of("/")+1, circuit.find_last_of(".") - circuit.find_last_of("/") - 1)<<"\t#Lines:\t"<<ckt.getLines().size()<<"\t#CNOTs:\t"<<CNOTs<<"\t#UB(c):\t"<<average_bound<<"\t#Map_Time:\t"<<average_map_time;//<<std::endl;
  min_ckt2.displayQASM(ofs2, dim.first*dim.second);
  unsigned CNOTs_II = min_ckt2.twoQubitGateCount() - CNOTs;
  std::cout<<"\t#A2_CNOTs:\t"<<CNOTs_II<<"\t#A2_Depth:\t"<<computeDepth(min_ckt2)<<"\t#A2_LA:\t"<<minw_a2<<"\t#A2_Time:\t"<<min_A2_time;  
  min_ckt2.clear();
  
  min_ckt3.displayQASM(ofs3, dim.first*dim.second); 
  unsigned CNOTs_III = min_ckt3.twoQubitGateCount() - CNOTs;
  std::cout<<"\t#A3_CNOTs:\t"<<CNOTs_III<<"\t#A3_Depth:\t"<<computeDepth(min_ckt3)<<"\t#A3_LA:\t"<<minw_a3<<"\t#A3_Time:\t"<<min_A3_time<<std::endl;  
  min_ckt3.clear();  
} 

