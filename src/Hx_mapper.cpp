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
//#define WIDTH 17      //hexagonal grid width
//#define HEIGHT 17     //hexagonal grid height
#define BASE 2   //edge weight of IG is power(BASE, depth - pos);
rev::window_t window = 0;
int main( int argc, char ** argv ) 
{ 
  //std::string layout;
  std::string circuit;
  std::ofstream ofs;
  //int base;
  //int i,j;
  int width, height;
  Metric_Type metric;  //RCNOT = 0 and Coupling Cost = 1 (DEFAULT = 0)
  Approach_Type approach;//RCNOT = 0, SWAP = 1 and SWAP + RCNOT = 2 (DEFAULT = 0)
  if (argc != 7)
  {
    std::cout<<"Usage:\n\t 1. circuitname[*.qasm] \n\t 2. outcircuitname[*.qasm] \n\t 3. width (hexagonal layout) \n\t 4. height (hexagonal layout) \n\t 5. metric (RCNOT = 0 and Coupling Cost = 1) \n\t 6. approach (RCNOT = 0, SWAP = 1 and SWAP + RCNOT = 2)"<<std::endl;
    std::exit(EXIT_FAILURE);
  }
  else{
    //layout = argv[1];
    circuit = argv[1];
    ofs.open(argv[2]);
    //base = atoi(argv[4]);
    width = atoi(argv[3]);
    height = atoi(argv[4]);
    metric = atoi(argv[5]) == 1 ?  Metric_Type::COUPLING_COST :  Metric_Type::RCNOT_COST;
    approach = atoi(argv[6]) == 2 ?  Approach_Type::A3 : atoi(argv[6]) == 1 ? Approach_Type::A2 : Approach_Type::A1;
  }
   
  //grph::Graph pg = grph::Graph();
  //pg.readLayout(layout);
  //std::cout<<circuit<<circuit.size()<<std::endl;
  rev::Circuit ckt = rev::Circuit();
  //std::cout<<"Hello"<<circuit<<std::endl;
  ckt.readQASM(circuit);
  //std::cout<<"Hello"<<std::endl;

  //unsigned long gates = ckt.getGates().size();
  //ckt.displayQASM(std::cout, ckt.getLines().size());  
  //ckt.displayQisKit(std::cout, width*height);
  window  = ckt.getCNOTDepth(); 
  //std::cout<<"Hello1"<<window<<std::endl; 
  if ( BASE == 2 && window > 20) window = 20;
  //if ( BASE == 3 && window > 80) window = 80;
  //if ( BASE == 4 && window > 60) window = 60;
  clock_t start, end;
  //clock_t start1, end1, start2, end2, start3, end3, start4, end4;
  //rev::Circuit min_ckt;
  //rev::Circuit tmp_ckt2;
  std::map<grph::node_t, grph::node_t> phmap_min = std::map<grph::node_t, grph::node_t> ( );
  //{{0, 1}, {1, 7}, {2, 12}, {3, 13}, {4, 19}}; //Mapping 2
  //{{0, 1}, {1, 11}, {2, 7}, {3, 6}, {4, 8}}; //->Mapping 1
  //{{0, 1}, {1, 7}, {2, 12}, {3, 18}}; //->Mapping 3
  
  double avgtime = 0.;
  start = clock(); 
  //start1 = clock(); 
  grph::Graph lg;
  if (metric == Metric_Type::COUPLING_COST)
    lg = grph::Graph(ckt, grph::NEW, window, BASE);
  else 
    lg = grph::Graph(ckt); //entire circuit
    //lg = grph::Graph(ckt, window); //fixed size window
    //lg = grph::Graph(ckt, 6, grph::window_type::variable); //variable size window, maxdeg for hexagonal architecture = 6
  //lg.displayGraph();
  std::map<grph::node_t, grph::node_t> phmap = std::map<grph::node_t, grph::node_t> ( );
  lg.mapLayout(phmap, std::pair<int, int> {width, height});

  
  /*std::cout<<grph::getDistance(0, 19, std::pair<int, int> {5, 4})<<std::endl;
  std::cout<<grph::getDistance(1, 19, std::pair<int, int> {5, 4})<<std::endl;
  std::cout<<grph::getDistance(5, 19, std::pair<int, int> {5, 4})<<std::endl;
  std::cout<<grph::getDistance(15, 4, std::pair<int, int> {5, 4})<<std::endl;
  std::cout<<grph::getDistance(15, 9, std::pair<int, int> {5, 4})<<std::endl;
 
  std::cout<<"Paths between 5 and 19"<<std::endl;
  std::vector<std::vector<grph::node_t> > paths;
  paths.push_back(std::vector<grph::node_t> ());
  getAllPaths(5, 19, paths, std::pair<int, int> {5, 4}, 0);
  for (auto nodes : paths){
    std::cout<<"Path: ";
    for (auto node : nodes){
      std::cout<<"Next Node: "<<node<<" ";
    }
    std::cout<<std::endl;
  }*/

  
  //lg.mapLayout(phmap, std::pair<int, int> {WIDTH, HEIGHT});
  

  for(const auto &  l : phmap) {
      phmap_min[l.first] = l.second;
      //std::cout<<"logical qubit: "<<l.first<<" \tPhysical qubit: "<<l.second<<std::endl;
  }
  //end1 = clock();
    
  //start2 = clock();*/
  //void  ga_search(std::map<grph::node_t, grph::node_t>& phmap, std::pair<int, int> dim, grph::Graph& lg);
  //need updation based on coupling cost or remote cnot cost
  //std::cout<<"global ordering"<<std::endl;
  ga_search(phmap_min, std::pair<int, int> {width, height}, lg, metric); //Global Ordering, Cost Metric: RCNOT, Coupling Cost
  
  //end2 = clock();

  /*for(const auto &  l : phmap_min) {
      std::cout<<"logical qubit: "<<l.first<<" \tPhysical qubit: "<<l.second<<std::endl;
  }*/

  //start3 = clock(); 
  //std::cout<<"local ordering"<<std::endl;
  switch (approach){
    case Approach_Type::A2:
      //Approach 2 for both RCNOT cost and Coupling cost
      //Approach 2: Global Ordering, Local Ordering, SWAP gate
      //std::cout<<"local ordering 2"<<std::endl;
      //lg = grph::Graph(ckt, window); //fixed size window
      //lg = grph::Graph(ckt, 6, grph::window_type::variable); //variable size window, maxdeg for hexagonal architecture = 6
      localorder(ckt, phmap_min, std::pair<int, int> {width, height}, lg, BASE, metric); 
      break;
    case Approach_Type::A3:
      //Approach 3 for both RCNOT cost and Coupling cost
      //Approach 3: Global Ordering, Local Ordering, SWAP + RCNOT gate      
      //std::cout<<"local ordering 3"<<std::endl;
      //lg = grph::Graph(ckt, window); //fixed size window
      //lg = grph::Graph(ckt, 6, grph::window_type::variable); //variable size window, maxdeg for hexagonal architecture = 6
      localorderv2(ckt, phmap_min, std::pair<int, int> {width, height}, lg, BASE, metric);
      break;
    default: 
      //Approach 1 for both RCNOT cost and Coupling cost
      //Approach 1: Global Ordering, RCNOT gate
      localorder(ckt, phmap_min, std::pair<int, int> {width, height}); 
  }
  end = clock();  
  //ckt.displayQASM(std::cout, width*height);
  //ckt.displayQisKit(std::cout, width*height);
  ckt.displayQASM(ofs, width*height); 
  //ckt.displayQisKit(ofs, width*height); 
  /*start4 = clock();
  optimize(tmp_ckt);
  end4 = clock();
    
  end = clock();
  if (k == 0 || tmp_ckt.getGates().size() < min_ckt.getGates().size()){
    min_ckt = tmp_ckt;
  }*/
  avgtime += ( (end - start) / (double) CLOCKS_PER_SEC );
   
  //min_ckt.display(std::cout, pg.getNodes().size());  
  //min_ckt.displayQisKit(ofs, pg.getNodes().size()); */
  
  std::cout<<"#Circuit:\t"<<circuit.substr(circuit.find_last_of("/")+1, circuit.find_last_of(".") - circuit.find_last_of("/") - 1)<<"\t#Lines:\t"<<ckt.getLines().size()<<"\t#Gates:\t"<<ckt.getGates().size()<<"\t#Depth:\t"<<computeDepth(ckt)<<"\t#Time:\t"<< avgtime / AVG <<std::endl;  
  //std::cout<<"#Circuit:\t"<<circuit.substr(circuit.find_last_of("/")+1, circuit.find_last_of(".") - circuit.find_last_of("/") - 1)<<"\t#Lines:\t"<<ckt.getLines().size()<<"\t#Gates:\t"<<ckt.getGates().size()<<"\t#Depth:\t"<<computeDepth(ckt)<<std::endl;
} 

