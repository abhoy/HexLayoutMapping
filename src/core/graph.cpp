/* ========================================================================== */
/*                                                                            */
/*   Graph.cpp                                                                */
/*   (c) 2015 Author                                                          */
/*                                                                            */
/*   Graph to represent quantum circuit                                       */
/*                                                                            */
/* ========================================================================== */
#include<iostream>
#include <sstream>
#include <fstream>
#include <cstdlib> 
#include "graph.hpp"
#include <cmath>

using namespace grph;

/*Coordinate::Coordinate(Coordinate& pos){
      this->x = pos.x;
      this->y = pos.y;
      this->z = pos.z;
}*/

int grph::roundNum(double num){
  if ( num + 0.5 >= (int) num + 1 ) return (int) num + 1;
  else return (int) num;
}


Coordinate::Coordinate(int x, int y, int z){
      this->x = x;
      this->y = y;
      this->z = z;
}

Coordinate::Coordinate(){
      this->x = 0;
      this->y = 0;
      this->z = 0;
}

bool Coordinate::isIdentical(Coordinate& pos){
  return this->x == pos.x && this->y == pos.y && this->z == pos.z;
}

void Graph::init(){
  this->nodes = std::vector<node_t> (); //all nodes
  this->edges = std::map<node_t, std::vector<node_t> > (); //node to edges mapping
  this->redges = std::map<node_t, std::vector<node_t> > (); //reverse node to edges mapping
  this->weight = std::map<node_t, std::map<node_t, weight_t> > (); //edge weight
  this->distance = std::map<node_t, std::map<node_t, distance_t> > (); //distance
  this->path = std::map<node_t, std::map<node_t, node_t> > (); //distance
  this->path2 = std::map<node_t, std::map<node_t, node_t> > (); //distance
  this->reverse = std::map<node_t, std::map<node_t, reverse_t> > (); //reverse directed edge
  this->degs = std::map<node_t, deg_t > (); //reverse directed edge
}

Graph::Graph(){
  init();  
}

Graph::Graph(const Graph & g){
  init();
  for(std::vector<grph::node_t>::const_iterator n = g.nodes.begin(); n != g.nodes.end(); ++n)
    this->nodes.push_back(*n);
  for (std::map<node_t, std::vector<node_t> >::const_iterator it = g.edges.begin(); it != g.edges.end(); ++it)
    for(std::vector<node_t>::const_iterator n = it->second.begin(); n != it->second.end(); ++n)
      this->edges[it->first].push_back(*n);
  for (std::map<node_t, std::vector<node_t> >::const_iterator it = g.redges.begin(); it != g.redges.end(); ++it)
    for(std::vector<node_t>::const_iterator n = it->second.begin(); n != it->second.end(); ++n)
      this->redges[it->first].push_back(*n);
  for (std::map<node_t, std::map<node_t, weight_t> >::const_iterator it = g.weight.begin(); it != g.weight.end(); ++it)
    for(std::map<node_t, weight_t>::const_iterator n = it->second.begin(); n != it->second.end(); ++n)
      this->weight[it->first][n->first] = n->second;
  for (std::map<node_t, std::map<node_t, distance_t> >::const_iterator it = g.distance.begin(); it != g.distance.end(); ++it)
    for(std::map<node_t, distance_t>::const_iterator n = it->second.begin(); n != it->second.end(); ++n)
      this->distance[it->first][n->first] = n->second;
  for (std::map<node_t, std::map<node_t, node_t> >::const_iterator it = g.path.begin(); it != g.path.end(); ++it)
    for(std::map<node_t, node_t>::const_iterator n = it->second.begin(); n != it->second.end(); ++n)
      this->path[it->first][n->first] = n->second;
  for (std::map<node_t, std::map<node_t, node_t> >::const_iterator it = g.path2.begin(); it != g.path2.end(); ++it)
    for(std::map<node_t, node_t>::const_iterator n = it->second.begin(); n != it->second.end(); ++n)
      this->path2[it->first][n->first] = n->second;
  for (std::map<node_t, std::map<node_t, reverse_t> >::const_iterator it = g.reverse.begin(); it != g.reverse.end(); ++it)
    for(std::map<node_t, reverse_t>::const_iterator n = it->second.begin(); n != it->second.end(); ++n)
      this->reverse[it->first][n->first] = n->second;
}

void Graph::operator = (const Graph & g){
  init();
  for(std::vector<grph::node_t>::const_iterator n = g.nodes.begin(); n != g.nodes.end(); ++n)
    this->nodes.push_back(*n);
  for (std::map<node_t, std::vector<node_t> >::const_iterator it = g.edges.begin(); it != g.edges.end(); ++it)
    for(std::vector<node_t>::const_iterator n = it->second.begin(); n != it->second.end(); ++n)
      this->edges[it->first].push_back(*n);
  for (std::map<node_t, std::vector<node_t> >::const_iterator it = g.redges.begin(); it != g.redges.end(); ++it)
    for(std::vector<node_t>::const_iterator n = it->second.begin(); n != it->second.end(); ++n)
      this->redges[it->first].push_back(*n);
  for (std::map<node_t, std::map<node_t, weight_t> >::const_iterator it = g.weight.begin(); it != g.weight.end(); ++it)
    for(std::map<node_t, weight_t>::const_iterator n = it->second.begin(); n != it->second.end(); ++n)
      this->weight[it->first][n->first] = n->second;
  for (std::map<node_t, std::map<node_t, distance_t> >::const_iterator it = g.distance.begin(); it != g.distance.end(); ++it)
    for(std::map<node_t, distance_t>::const_iterator n = it->second.begin(); n != it->second.end(); ++n)
      this->distance[it->first][n->first] = n->second;
  for (std::map<node_t, std::map<node_t, node_t> >::const_iterator it = g.path.begin(); it != g.path.end(); ++it)
    for(std::map<node_t, node_t>::const_iterator n = it->second.begin(); n != it->second.end(); ++n)
      this->path[it->first][n->first] = n->second;
  for (std::map<node_t, std::map<node_t, node_t> >::const_iterator it = g.path2.begin(); it != g.path2.end(); ++it)
    for(std::map<node_t, node_t>::const_iterator n = it->second.begin(); n != it->second.end(); ++n)
      this->path2[it->first][n->first] = n->second;
  for (std::map<node_t, std::map<node_t, reverse_t> >::const_iterator it = g.reverse.begin(); it != g.reverse.end(); ++it)
    for(std::map<node_t, reverse_t>::const_iterator n = it->second.begin(); n != it->second.end(); ++n)
      this->reverse[it->first][n->first] = n->second;
}

bool Graph::isEmpty(){
  if (this->nodes.size() == 0) return true;
  weight_t cum_w = 0.;
  for (auto & n : nodes){
    for (auto & m : nodes){
      if (n != m) cum_w += getWeight(n, m);
    }
  }
  return cum_w <= 0.;
}
//Creating general graph from a circuit where nodes are qubits and edge weight represents the number
//2-qubit gate operations
Graph::Graph(rev::Circuit& ckt){
  init();
  
  
  //creating vertex v for each qubit line present in the circuit
  //std::map<rev::line_t, double> CNOTDepth = std::map<rev::line_t, double> (); 
  std::vector<rev::line_t> lines = ckt.getLines(); 
  for(std::vector<rev::line_t>::iterator i=lines.begin(); i!=lines.end(); ++i){
    this->addNode(*i);  //Adding nodes to adjacency graph
    //CNOTDepth[*i] = 0;
  }
 
  
  //creating an weighted edge e for each two input gate present in the circuit
  std::vector<rev::Gate> gates = ckt.getGates();
  //int k=0;
  for(std::vector<rev::Gate>::iterator i=gates.begin(); i!=gates.end(); ++i){
    //std::cout<<"Hello"<<i->getType()<<std::endl;
    if(rev::CX == i->getType()){
      //std::cout<<"control:"<<i->getControls()[0]<<" target:"<<i->getTargets()[0]<<std::endl;
      this->addEdge2(i->getControls()[0], i->getTargets()[0], 0, 0, 0);       
    } 
  } 
}

//qubit interaction graph based on initial constant size window   
Graph::Graph(rev::Circuit& ckt, window_t w){
  init();
  
  
  //creating vertex v for each qubit line present in the circuit
  std::map<rev::line_t, window_t> CNOTDepth = std::map<rev::line_t, window_t> (); 
  std::vector<rev::line_t> lines = ckt.getLines(); 
  for(std::vector<rev::line_t>::iterator i=lines.begin(); i!=lines.end(); ++i){
    this->addNode(*i);  //Adding nodes to adjacency graph
    this->degs[*i] = 0;
    CNOTDepth[*i] = 0;
  }
 
  
  //creating an weighted edge e for each two input gate present in the circuit
  std::vector<rev::Gate> gates = ckt.getGates();
  //int k=0;
  for(std::vector<rev::Gate>::iterator i=gates.begin(); i!=gates.end(); ++i){
    //std::cout<<"Hello"<<i->getType()<<std::endl;
    if(rev::CX == i->getType()){
      //std::cout<<"control:"<<i->getControls()[0]<<" target:"<<i->getTargets()[0]<<std::endl;
      window_t cdepth = CNOTDepth[i->getControls()[0]] > CNOTDepth[i->getTargets()[0]] ? 
                                     CNOTDepth[i->getControls()[0]] + 1 : CNOTDepth[i->getTargets()[0]] + 1;
      if (cdepth > w) break;      
      this->degs[i->getControls()[0]] += 1;
      this->degs[i->getTargets()[0]] += 1;
      CNOTDepth[i->getControls()[0]] = cdepth;
      CNOTDepth[i->getTargets()[0]] = cdepth; 
      this->addEdge2(i->getControls()[0], i->getTargets()[0], 0, 0, 0);       
    } 
  } 
}

//qubit interaction graph based on initial variable size window   
Graph::Graph(rev::Circuit& ckt, deg_t maxdeg, window_type type){
  init();
  
  if (type == window_type::fixed){
    std::cout<<"Window type will be variable"<<std::endl;
    exit(0);
  }
  //creating vertex v for each qubit line present in the circuit
  std::map<rev::line_t, window_t> CNOTDepth = std::map<rev::line_t, window_t> (); 
  std::vector<rev::line_t> lines = ckt.getLines(); 
  for(std::vector<rev::line_t>::iterator i=lines.begin(); i!=lines.end(); ++i){
    this->addNode(*i);  //Adding nodes to adjacency graph
    this->degs[*i] = 0;
    CNOTDepth[*i] = 0;
  }
 
  
  //creating an weighted edge e for each two input gate present in the circuit
  std::vector<rev::Gate> gates = ckt.getGates();
  //int k=0;
  for(std::vector<rev::Gate>::iterator i=gates.begin(); i!=gates.end(); ++i){
    //std::cout<<"Hello"<<i->getType()<<std::endl;
    if(rev::CX == i->getType()){
      //std::cout<<"control:"<<i->getControls()[0]<<" target:"<<i->getTargets()[0]<<std::endl;
      deg_t cdeg = this->degs[i->getControls()[0]] > this->degs[i->getTargets()[0]] ? 
                                     this->degs[i->getControls()[0]] + 1 : this->degs[i->getTargets()[0]] + 1;
      window_t cdepth = CNOTDepth[i->getControls()[0]] > CNOTDepth[i->getTargets()[0]] ? 
                                     CNOTDepth[i->getControls()[0]] + 1 : CNOTDepth[i->getTargets()[0]] + 1;
      if (cdeg > maxdeg) break;      
      this->degs[i->getControls()[0]] += 1;
      this->degs[i->getTargets()[0]] += 1;
      CNOTDepth[i->getControls()[0]] = cdepth;
      CNOTDepth[i->getTargets()[0]] = cdepth; 
      this->addEdge2(i->getControls()[0], i->getTargets()[0], 0, 0, 0);       
    } 
  } 
}
Graph::Graph(rev::Circuit& ckt, nnc_t type, window_t w, base_t x){
  init();
  
  
  //creating vertex v for each qubit line present in the circuit
  std::map<rev::line_t, window_t> CNOTDepth = std::map<rev::line_t, window_t> (); 
  std::vector<rev::line_t> lines = ckt.getLines(); 
  for(std::vector<rev::line_t>::iterator i=lines.begin(); i!=lines.end(); ++i){
    this->addNode(*i);  //Adding nodes to adjacency graph
    CNOTDepth[*i] = 0;
  }
 
  
  //creating an weighted edge e for each two input gate present in the circuit
  std::vector<rev::Gate> gates = ckt.getGates();
  //int k=0;
  for(std::vector<rev::Gate>::iterator i=gates.begin(); i!=gates.end(); ++i)
    if(rev::CX == i->getType()){
      window_t cdepth = CNOTDepth[i->getControls()[0]] > CNOTDepth[i->getTargets()[0]] ? 
                                     CNOTDepth[i->getControls()[0]] + 1 : CNOTDepth[i->getTargets()[0]] + 1;
      if (cdepth > w) break;
      CNOTDepth[i->getControls()[0]] = cdepth;
      CNOTDepth[i->getTargets()[0]] = cdepth; 
      if(OLD == type) this->addEdge(i->getControls()[0], i->getTargets()[0]); 
      else this->addEdge2(i->getControls()[0], i->getTargets()[0], x, w, cdepth);//k++); 
      //if( k == w) break;  // break the loop when window size is reached
    }  
}


Graph::Graph(std::vector<rev::line_t>& lines, std::vector<rev::Gate>& gates, int index, nnc_t type, window_t w, base_t x){
  init();
  
  
  //creating vertex v for each qubit line present in the circuit
  std::map<rev::line_t, window_t> CNOTDepth = std::map<rev::line_t, window_t> (); 
  //std::vector<rev::line_t> lines = ckt.getLines(); 
  for(std::vector<rev::line_t>::iterator i=lines.begin(); i!=lines.end(); ++i){
    this->addNode(*i);  //Adding nodes to adjacency graph
    CNOTDepth[*i] = 0;
  }
  
  //creating an weighted edge e for each two input gate present in the circuit
  //std::vector<rev::Gate> gates = ckt.getGates();
  //std::cout<<"Hello<<"<<std::endl;
  //int k=0;
  //std::cout<<"window:"<<w<<std::endl;
  for(std::vector<rev::Gate>::iterator i=gates.begin() + index; i!=gates.end(); ++i){    
    //(*i).display(std::cout);
    if(rev::CX == i->getType()){  
      window_t cdepth = CNOTDepth[i->getControls()[0]] > CNOTDepth[i->getTargets()[0]] ? 
                                     CNOTDepth[i->getControls()[0]] + 1 : CNOTDepth[i->getTargets()[0]] + 1;
      if (cdepth > w) break;
      CNOTDepth[i->getControls()[0]] = cdepth;
      CNOTDepth[i->getTargets()[0]] = cdepth; 
      if(OLD == type) this->addEdge(i->getControls()[0], i->getTargets()[0]); 
      else this->addEdge2(i->getControls()[0], i->getTargets()[0], x, w, cdepth);
      //if( k == w) break;  // break the loop when window size is reached
    }  
  }
  //std::cout<<"Hello"<<std::endl;

}

// Constant size window
Graph::Graph(std::vector<rev::line_t>& lines, std::vector<rev::Gate>& gates, int index, window_t w){
  init();
  
  
  //creating vertex v for each qubit line present in the circuit
  std::map<rev::line_t, window_t> CNOTDepth = std::map<rev::line_t, window_t> (); 
  //std::vector<rev::line_t> lines = ckt.getLines(); 
  for(std::vector<rev::line_t>::iterator i=lines.begin(); i!=lines.end(); ++i){
    this->addNode(*i);  //Adding nodes to adjacency graph
    CNOTDepth[*i] = 0;
  }
  
  //creating an weighted edge e for each two input gate present in the circuit
  //std::vector<rev::Gate> gates = ckt.getGates();
  //std::cout<<"Hello<<"<<std::endl;
  //int k=0;
  //std::cout<<"window:"<<w<<std::endl;
  for(std::vector<rev::Gate>::iterator i=gates.begin() + index; i!=gates.end(); ++i){    
    //(*i).display(std::cout);
    if(rev::CX == i->getType()){  
      window_t cdepth = CNOTDepth[i->getControls()[0]] > CNOTDepth[i->getTargets()[0]] ? 
                                     CNOTDepth[i->getControls()[0]] + 1 : CNOTDepth[i->getTargets()[0]] + 1;
      if (cdepth > w) break;
      CNOTDepth[i->getControls()[0]] = cdepth;
      CNOTDepth[i->getTargets()[0]] = cdepth; 
      this->addEdge2(i->getControls()[0], i->getTargets()[0], 0, 0, 0);       
      //if( k == w) break;  // break the loop when window size is reached
    }  
  }
  //std::cout<<"Hello"<<std::endl;

}

// Variable size window
Graph::Graph(std::vector<rev::line_t>& lines, std::vector<rev::Gate>& gates, int index, deg_t maxdeg, window_type type){
  init();
  
  
  //creating vertex v for each qubit line present in the circuit
  std::map<rev::line_t, window_t> CNOTDepth = std::map<rev::line_t, window_t> (); 
  //std::vector<rev::line_t> lines = ckt.getLines(); 
  for(std::vector<rev::line_t>::iterator i=lines.begin(); i!=lines.end(); ++i){
    this->addNode(*i);  //Adding nodes to adjacency graph
    this->degs[*i] = 0;
    CNOTDepth[*i] = 0;
  }
  
  //creating an weighted edge e for each two input gate present in the circuit
  //std::vector<rev::Gate> gates = ckt.getGates();
  //std::cout<<"Hello<<"<<std::endl;
  //int k=0;
  //std::cout<<"window:"<<w<<std::endl;
  for(std::vector<rev::Gate>::iterator i=gates.begin() + index; i!=gates.end(); ++i){    
    //(*i).display(std::cout);
    if(rev::CX == i->getType()){  
      deg_t cdeg = this->degs[i->getControls()[0]] > this->degs[i->getTargets()[0]] ? 
                                     this->degs[i->getControls()[0]] + 1 : this->degs[i->getTargets()[0]] + 1;
      window_t cdepth = CNOTDepth[i->getControls()[0]] > CNOTDepth[i->getTargets()[0]] ? 
                                     CNOTDepth[i->getControls()[0]] + 1 : CNOTDepth[i->getTargets()[0]] + 1;
      if (cdeg > maxdeg) break;      
      this->degs[i->getControls()[0]] += 1;
      this->degs[i->getTargets()[0]] += 1;
      CNOTDepth[i->getControls()[0]] = cdepth;
      CNOTDepth[i->getTargets()[0]] = cdepth; 
      this->addEdge2(i->getControls()[0], i->getTargets()[0], 0, 0, 0);       
      //if( k == w) break;  // break the loop when window size is reached
    }  
  }
  //std::cout<<"Hello"<<std::endl;

}


void Graph::addNode(node_t nd){
  this->nodes.push_back(nd);
}

void Graph::addEdge(node_t nd1, node_t nd2){ //For Qubit Coupling Graph
  //std::map<node_t, map<node_t, weight_t>>:: const_iterator e;
  //e = this->weight.find(nd2);
  if(this->weight[nd1].end() == this->weight[nd1].find(nd2)){
    this->edges[nd1].push_back(nd2);   // nd2 is a neighbour of nd1
    this->edges[nd2].push_back(nd1);   // nd2 is a neighbour of nd2
    //this->redges[nd2].push_back(nd1);
    this->weight[nd1][nd2] = 1;
    //this->weight[nd2][nd1] = 1;        
  }
  else{
    this->weight[nd1][nd2]++;
    //this->weight[nd2][nd1]++;
  }
  this->distance[nd1][nd2] = FCOST;
  //std::cout<<"Path: "<<nd1<<" : "<<nd2<<std::endl;
  this->path[nd1][nd2] = nd2;
  this->distance[nd2][nd1] = FCOST;
}

void Graph::addDirEdge(node_t nd1, node_t nd2){
  if(this->weight[nd1].end() == this->weight[nd1].find(nd2)){
    this->edges[nd1].push_back(nd2);
    this->redges[nd2].push_back(nd1);
    this->weight[nd1][nd2] = 1;
    this->distance[nd1][nd2] = FCOST;
    //std::cout<<"Path: "<<nd1<<" : "<<nd2<<std::endl;
    this->path[nd1][nd2] = nd2;
    this->distance[nd2][nd1] = RCOST;
    //std::cout<<"Path: "<<nd2<<" : "<<nd1<<std::endl;
    this->path[nd2][nd1] = nd1;
    this->reverse[nd1][nd2] = false;
    this->reverse[nd2][nd1] = true;
  }
}

void Graph::updateEdgeWeight(node_t n, node_t m, double weight){
  this->weight[n][m] -= weight;
}


void Graph::addEdge2(node_t nd1, node_t nd2, base_t x, window_t y, int i){ //For Qubit Interaction Graph
  weight_t w;
  if (x == 0 && y == 0 && i == 0) //simple weight = number of gates
    w = 1;
  else
    w = std::pow(x, y - i);// - 1); //(y/2) - i
  //std::cout<<"w="<<w<<" x="<<x<<" y="<<y<<" i="<<i<<std::endl;
  //std::cout<<"Weight: "<<std::setprecision(5)<<w<<std::endl;
  if(this->weight[nd1].end() == this->weight[nd1].find(nd2) && this->weight[nd2].end() == this->weight[nd2].find(nd1)){
    this->edges[nd1].push_back(nd2);
     this->redges[nd2].push_back(nd1);
    //this->edges[nd2].push_back(nd1);
    this->weight[nd1][nd2] = w;
    //this->weight[nd2][nd1] = w;
    //std::cout<<"Here 1"<<std::endl;
  }
  else{
    if (this->weight[nd1].end() == this->weight[nd1].find(nd2))
      this->weight[nd2][nd1] += w;
    else 
      this->weight[nd1][nd2] += w;
    //this->weight[nd2][nd1] += w;
    //std::cout<<"Here 2"<<std::endl;
  }
  //std::cout<<nd1<<", "<<nd2<<": "<<this->weight[nd1][nd2]<<std::endl;
}

void Graph::removeCNOTEdge(rev::Gate & g, window_t w, base_t x){
  //node_t control = g.getControls()[0];
  //node_t target = g.getTargets()[0];
  //std::vector<node_t> controls = std::vector<node_t> ();
  //std::vector<nod_t> targets = std::vector<node_t> ();
  double cweight = std::pow(x, w-1);
  for(std::vector<node_t>::iterator i=this->nodes.begin(); i!=this->nodes.end(); ++i){
    for(std::vector<node_t>::iterator j=this->edges[*i].begin(); j!=this->edges[*i].end(); ++j){
      if (this->weight[*i][*j] >= cweight) {
        //controls.push_back(*i);
        //targets.push_back(*j);
        this->weight[*i][*j] -= cweight;
        if (this->weight[*i][*j] == 0) this->weight[*i].erase(*j);
      }
    }
  } 
}

weight_t Graph::getWeight(node_t nd1, node_t nd2){
  if (this->weight[nd1].end() == this->weight[nd1].find(nd2))
    return this->weight[nd2][nd1];
  else
    return this->weight[nd1][nd2];
}

/*weight_t Graph::getWeight(node_t nd1, node_t nd2){
  //if(this->weight.end() == this->weight.find(nd1))
  //  return -1;
  //else 
  return this->weight[nd1][nd2];
}*/

void Graph::setWeight(node_t nd1, node_t nd2, weight_t w){
  if (this->weight[nd1].end() == this->weight[nd1].find(nd2))
    this->weight[nd2][nd1] = w;
  else 
    this->weight[nd1][nd2] = w;
}

void Graph::displayGraph(){
  std::cout<<"Adjacent Nodes:"<<std::endl;
  std::cout<<"Node:"<<"\t"<<"Adjacent Nodes (Weight)"<<std::endl;
  for(std::vector<node_t>::iterator i=this->nodes.begin(); i!=this->nodes.end(); ++i){
    //std::cout<<"Node: "<<(*i)<<" Adjacent Nodes (Weight)";
    std::cout<<setw(3)<<(*i);
    for(std::vector<node_t>::iterator j=this->edges[*i].begin(); j!=this->edges[*i].end(); ++j)
      std::cout<<"\t"<<setw(3)<<(*j)<<"("<<setw(4)<<this->weight[*i][*j]<<")";
    std::cout<<std::endl;
  }
}
std::vector<node_t>& Graph::getNodes(){
  return this->nodes;
  }

std::map<node_t, std::vector<node_t> >& Graph::getEdges(){
  return this->edges;
}

std::map<node_t, std::vector<node_t> >& Graph::getRedges(){
  return this->redges;
}
      
std::map<node_t, std::map<node_t, weight_t> >& Graph::getEdgeWeights(){
  return this->weight;
}
node_t Graph::getNodeInPath(node_t p1, node_t p2){
  return this->path[p1][p2];
}
node_t Graph::getNodeInPath2(node_t p1, node_t p2){
  return this->path2[p1].end() == this->path2[p1].find(p2) ? -1 : this->path2[p1][p2];//this->path[p1][p2];
}
distance_t Graph::getDistance(node_t p1, node_t p2){
  return this->distance[p1][p2];
}
//distance for hexagonal architecture
distance_t grph::getDistance(node_t p1, node_t p2, std::pair<int, int> dim){
  // coordinate of node p1
  int y_p1 = p1 / dim.first;
  int x_p1 = (1 - (y_p1 % 2 )) + 2 * (p1 % dim.first);

  // coordinate of node p2
  int y_p2 = p2 / dim.first;
  int x_p2 = (1 - (y_p2 % 2 )) + 2 * (p2 % dim.first);
  
  int d1 = abs(y_p1 - y_p2);
  
  int d2 = (abs(x_p1 - x_p2) +  d1 ) / 2;
  //std::cout<<p1<<": ("<<x_p1<<", "<<y_p1<<"), "<<p2<<": ("<<x_p2<<", "<<y_p2<<"), Distance: ";//<<std::endl;
  return d1 > d2 ? d1 - 1 : d2 - 1; //max of d1 - 1 , d2 - 1
}

void Graph::computePath(){
  for(std::vector<node_t>::iterator k=this->nodes.begin(); k!=this->nodes.end(); ++k){
    //std::cout<<"Level:"<<*k<<std::endl;
    for(std::vector<node_t>::iterator i=this->nodes.begin(); i!=this->nodes.end(); ++i){
      for(std::vector<node_t>::iterator j=this->nodes.begin(); j!=this->nodes.end(); ++j){
        //std::cout<<(this->distance[*i].end() == this->distance[*i].find(*j) ? -1 : this->distance[*i][*j])<<"\t";
        if(*i == *j || *i == *k || *k == *j) continue;
        if (this->distance[*i][*k] <= RCOST){
          //if (*i == 0 && *j == 3) std::cout<<"Hello k="<<*k<<" and D="<<this->distance[*k][*j]<<std::endl;
          if (this->path[*i].end() == this->path[*i].find(*j))
            this->path[*i][*j] = *k;
          else if (this->distance[*k][*j] < this->distance[this->path[*i][*j]][*j]){
            if (this->distance[*k][*j] + RCOST >= this->distance[this->path[*i][*j]][*j])
              this->path2[*i][*j] = this->path[*i][*j]; 
            else if (this->path2[*i].end() != this->path2[*i].find(*j))
              this->path2[*i].erase(this->path2[*i].find(*j));          
            this->path[*i][*j] = *k;
          }
          else if (this->distance[*k][*j] - this->distance[this->path[*i][*j]][*j] <= RCOST && this->path[*i][*j] != *k){
            this->path2[*i][*j] = *k;           
          }   
        }
      }
    }
  }
}

void Graph::computeDistance(){
  for(std::vector<node_t>::iterator k=this->nodes.begin(); k!=this->nodes.end(); ++k){
    //std::cout<<"Level:"<<*k<<std::endl;
    for(std::vector<node_t>::iterator i=this->nodes.begin(); i!=this->nodes.end(); ++i){
      for(std::vector<node_t>::iterator j=this->nodes.begin(); j!=this->nodes.end(); ++j){
        //std::cout<<(this->distance[*i].end() == this->distance[*i].find(*j) ? -1 : this->distance[*i][*j])<<"\t";
        if(*i == *j || *i == *k || *k == *j) continue;
	if(this->distance[*i].end() != this->distance[*i].find(*k) && this->distance[*k].end() != this->distance[*k].find(*i) &&
           this->distance[*k].end() != this->distance[*k].find(*j) && this->distance[*j].end() != this->distance[*j].find(*k)){	
          double min_ik = this->distance[*i][*k] > this->distance[*k][*i] ? this->distance[*k][*i] : this->distance[*i][*k];
          double min_kj = this->distance[*k][*j] > this->distance[*j][*k] ? this->distance[*j][*k] : this->distance[*k][*j];
          double cur_distance = min_ik + min_kj;
          //if (*i == 11 && *j==9) std::cout<<"distance="<<cur_distance;
          if (this->distance[*k][*i] == min_ik && this->distance[*k][*j] > min_kj) //this->distance[*j][*k] == min_kj &&
            cur_distance += (1 + RCOST);
          else
            cur_distance += 1;
          //if (*i == 11 && *j==9) std::cout<<"distance="<<cur_distance<<std::endl;    
          /*if (*i == 11 && *j==9) { 
            std::cout<<"K="<<*k<<" min_ik="<<min_ik<<" min_kj="<<min_kj; 
            std::cout<<" d(i,k)="<<this->distance[*i][*k]<<" d(k,i)="<<this->distance[*k][*i]; 
            std::cout<<" d(k,j)="<<this->distance[*k][*j]<<" d(j,k)="<<this->distance[*j][*k]; 
          }*/
          if((this->distance[*i].end() == this->distance[*i].find(*j))){ //distance is infinity
            this->distance[*i][*j] = cur_distance; 
            //std::cout<<"Path: "<<*i<<" : "<<*j<<std::endl;  
            //this->path[*i][*j] = *k;
            this->reverse[*i][*j] = false;
            //if (*i == 11 && *j==9) std::cout<<" UPDATED HERE 1"<<std::endl;
	  }
          else if (this->distance[*i][*j] > cur_distance) {
            this->distance[*i][*j] = cur_distance;
            //std::cout<<"Path: "<<*i<<" : "<<*j<<std::endl;  
            //this->path[*i][*j] = *k;
            this->reverse[*i][*j] = false;
            //if (*i == 11 && *j==9) std::cout<<" UPDATED HERE 2"<<std::endl;
	  }       
          //else if (this->distance[*i][*j] == cur_distance && this->path[*i][*j] != *k) {
            //std::cout<<"Path2: "<<*i<<" : "<<*j<<std::endl;  
            //this->path2[*i][*j] = *k;
            //if (*i == 9 && *j==6) std::cout<<"Forward3 "<<*k<<std::endl;
	  //}
          //if (*i == 11 && *j==9) std::cout<<" Distance(11,9)="<<this->distance[*i][*j]<<std::endl; 
	}
      }
    }
  }
  for(std::vector<node_t>::reverse_iterator k=this->nodes.rbegin(); k!=this->nodes.rend(); ++k){
    //std::cout<<"Level:"<<*k<<std::endl;
    for(std::vector<node_t>::iterator i=this->nodes.begin(); i!=this->nodes.end(); ++i){
      for(std::vector<node_t>::iterator j=this->nodes.begin(); j!=this->nodes.end(); ++j){
        //std::cout<<(this->distance[*i].end() == this->distance[*i].find(*j) ? -1 : this->distance[*i][*j])<<"\t";
        if(*i == *j || *i == *k || *k == *j) continue;
	if(this->distance[*i].end() != this->distance[*i].find(*k) && this->distance[*k].end() != this->distance[*k].find(*i) &&
           this->distance[*k].end() != this->distance[*k].find(*j) && this->distance[*j].end() != this->distance[*j].find(*k)){	
          double min_ik = this->distance[*i][*k] > this->distance[*k][*i] ? this->distance[*k][*i] : this->distance[*i][*k];
          double min_kj = this->distance[*k][*j] > this->distance[*j][*k] ? this->distance[*j][*k] : this->distance[*k][*j];
          double cur_distance = min_ik + min_kj;
          if (this->distance[*k][*i] == min_ik && this->distance[*k][*j] > min_kj) //this->distance[*j][*k] == min_kj &&
            cur_distance += (1 + RCOST);
          else
            cur_distance += 1;

          /*if (*i == 11 && *j==9) { 
            std::cout<<"K="<<*k<<" min_ik="<<min_ik<<" min_kj="<<min_kj; 
            std::cout<<" d(i,k)="<<this->distance[*i][*k]<<" d(k,i)="<<this->distance[*k][*i]; 
            std::cout<<" d(k,j)="<<this->distance[*k][*j]<<" d(j,k)="<<this->distance[*j][*k]; 
          }*/
          if(this->distance[*i].end() == this->distance[*i].find(*j)){//distance is infinity
            this->distance[*i][*j] = cur_distance;  
            //std::cout<<"Path: "<<*i<<" : "<<*j<<std::endl;   
            //this->path[*i][*j] = *k;
            this->reverse[*i][*j] = false;
            //if (*i == 11 && *j==9) std::cout<<" UPDATED HERE 3"<<std::endl;
	  }
          else if (this->distance[*i][*j] > cur_distance) {
            this->distance[*i][*j] = cur_distance;
            //std::cout<<"Path: "<<*i<<" : "<<*j<<std::endl;  
            //this->path[*i][*j] = *k;
            //this->path2[*i].erase(*j);
            this->reverse[*i][*j] = false; 
            //if (*i == 11 && *j==9) std::cout<<" UPDATED HERE 4"<<std::endl;
	  }
          //else if (this->distance[*i][*j] == cur_distance && this->path[*i][*j] != *k) {
            //std::cout<<"Path2: "<<*i<<" : "<<*j<<std::endl;  
            //this->path2[*i][*j] = *k;
            //if (*i == 9 && *j==6) std::cout<<"Reverse3 "<<*k<<std::endl;
	  //}
          //if (*i == 11 && *j==9) std::cout<<" Distance(11,9)="<<this->distance[*i][*j]<<std::endl; 
	}
      }
    }
  }
}
void Graph::computeDistance2(){
  for(std::vector<node_t>::iterator k=this->nodes.begin(); k!=this->nodes.end(); ++k){
	//std::cout<<" Distance ("<<*k<<", "<<*k<<" ) "<<this->distance[*k][*k]<<std::endl;
    for(std::vector<node_t>::iterator i=this->nodes.begin(); i!=this->nodes.end(); ++i){
      for(std::vector<node_t>::iterator j=this->nodes.begin(); j!=this->nodes.end(); ++j){
        //std::cout<<(this->distance[*i].end() == this->distance[*i].find(*j) ? -1 : this->distance[*i][*j])<<"\t";
        if(*i == *j || *i == *k || *k == *j) continue;
        //std::cout<<"i:"<<*i<<" j:"<<*j<<" k:"<<*k<<" distance "<<this->distance[*i][*j]<<std::endl;
	if(this->distance[*i].end() != this->distance[*i].find(*k) && this->distance[*k].end() != this->distance[*k].find(*i) &&
           this->distance[*k].end() != this->distance[*k].find(*j) && this->distance[*j].end() != this->distance[*j].find(*k)){	
          //double min_ik = this->distance[*i][*k] > this->distance[*k][*i] ? this->distance[*k][*i] : this->distance[*i][*k];
	  //std::cout<<"i="<<*i <<" k="<<*k<<" distance="<<min_ik	<<std::endl;
          //double min_kj = this->distance[*k][*j] > this->distance[*j][*k] ? this->distance[*j][*k] : this->distance[*k][*j];
          double cur_distance = this->distance[*i][*k] + this->distance[*k][*j] + 1;//min_ik + min_kj;
          //if (*i == 11 && *j==9) 
          //std::cout<<" distance="<<cur_distance<<std::endl;
          //if (this->distance[*k][*i] == min_ik && this->distance[*k][*j] > min_kj) //this->distance[*j][*k] == min_kj &&
            //cur_distance += (1 + RCOST);
          //else
            //cur_distance += 1;
          //if (*i == 11 && *j==9) std::cout<<"distance="<<cur_distance<<std::endl;    
          /*if (*i == 11 && *j==9) { 
            std::cout<<"K="<<*k<<" min_ik="<<min_ik<<" min_kj="<<min_kj; 
            std::cout<<" d(i,k)="<<this->distance[*i][*k]<<" d(k,i)="<<this->distance[*k][*i]; 
            std::cout<<" d(k,j)="<<this->distance[*k][*j]<<" d(j,k)="<<this->distance[*j][*k]; 
          }*/
          if((this->distance[*i].end() == this->distance[*i].find(*j))){ //distance is infinity
            this->distance[*i][*j] = cur_distance; 
            //std::cout<<"Path: "<<*i<<" : "<<*j<<std::endl;  
            //this->path[*i][*j] = *k;
            //this->reverse[*i][*j] = false;
            //if (*i == 11 && *j==9) std::cout<<" UPDATED HERE 1"<<std::endl;
	  }
          else if (this->distance[*i][*j] > cur_distance) {
            this->distance[*i][*j] = cur_distance;
            //std::cout<<"Path: "<<*i<<" : "<<*j<<std::endl;  
            //this->path[*i][*j] = *k;
            //this->reverse[*i][*j] = false;
            //if (*i == 11 && *j==9) std::cout<<" UPDATED HERE 2"<<std::endl;
	  }       
          //else if (this->distance[*i][*j] == cur_distance && this->path[*i][*j] != *k) {
            //std::cout<<"Path2: "<<*i<<" : "<<*j<<std::endl;  
            //this->path2[*i][*j] = *k;
            //if (*i == 9 && *j==6) std::cout<<"Forward3 "<<*k<<std::endl;
	  //}
          //if (*i == 11 && *j==9) std::cout<<" Distance(11,9)="<<this->distance[*i][*j]<<std::endl; 
	}
      }
    }
  }
}
void Graph::displayDistance(){
  std::cout<<"Shortest Distance:"<<std::endl;
  for(std::vector<node_t>::iterator i=this->nodes.begin(); i!=this->nodes.end(); ++i){
    //print matrix label
    if (*i == 0){
      std::cout<<"   ";
      for(std::vector<node_t>::iterator j=this->nodes.begin(); j!=this->nodes.end(); ++j)
        std::cout<<setw(4)<<*j<<"  ";
      std::cout<<std::endl;
    }
    std::cout<<setw(2)<<*i<<" ";
    //end of print
    for(std::vector<node_t>::iterator j=this->nodes.begin(); j!=this->nodes.end(); ++j)
      std::cout<<setw(4)<<(this->distance[*i].end() == this->distance[*i].find(*j) ? -1 : this->distance[*i][*j])<<"  ";
    std::cout<<std::endl;
  }
}
void Graph::displayPath(){
  std::cout<<"Shortest Path details:"<<std::endl;
  for(std::vector<node_t>::iterator i=this->nodes.begin(); i!=this->nodes.end(); ++i){
    //print matrix label
    if (*i == 0){
      std::cout<<"   ";
      for(std::vector<node_t>::iterator j=this->nodes.begin(); j!=this->nodes.end(); ++j)
        std::cout<<setw(2)<<*j<<"  ";
      std::cout<<std::endl;
    }
    std::cout<<setw(2)<<*i<<" ";
    //end of print
    for(std::vector<node_t>::iterator j=this->nodes.begin(); j!=this->nodes.end(); ++j)
      std::cout<<setw(2)<<(this->path[*i].end() == this->path[*i].find(*j) ? 99 : this->path[*i][*j])<<"  ";
      //std::cout<<this->path[*i][*j]<<std::endl;
    std::cout<<std::endl;
  }
  std::cout<<"Alternate Path details:"<<std::endl;
  for(std::vector<node_t>::iterator i=this->nodes.begin(); i!=this->nodes.end(); ++i){
    //print matrix label
    if (*i == 0){
      std::cout<<"   ";
      for(std::vector<node_t>::iterator j=this->nodes.begin(); j!=this->nodes.end(); ++j)
        std::cout<<setw(2)<<*j<<"  ";
      std::cout<<std::endl;
    }
    std::cout<<setw(2)<<*i<<" ";
    //end of print
    for(std::vector<node_t>::iterator j=this->nodes.begin(); j!=this->nodes.end(); ++j)
      std::cout<<setw(2)<<(this->path2[*i].end() == this->path2[*i].find(*j) ? 99 : this->path2[*i][*j])<<"  ";
    std::cout<<std::endl;
  }
}
void Graph::readLayout(string file){
  char fileName[50];
  int i, ln = file.length();
  for(i=0; i < ln; i++)
    fileName[i] = file[i];
  fileName[i] = '\0';
  string line;
  ifstream myfile(fileName);
  std::map<string, node_t> q_map;
  if(myfile.is_open()){
    while(getline(myfile,line)){ 
      if(0 == line.size() || 1 == line.size()) continue; //skip empty lines
      if('#' == line[0])  continue; //skip comments  
      istringstream iss(line);
      string command;
      iss >> command;
      if("Qubit:" == command){ //Physical Qubits
        while(iss >> command) {
		q_map[command] = nodes.size();
       		addNode(nodes.size());
	}
        continue;
      }
      if("end" == command) break; //End of Layout 	
      if("Coupling:" == command) continue;//Qubits Coupling
      //read all coupling map information
      node_t nd1 = q_map[command];
      //std::cout<<nd1<<command;
      iss >> command;
      node_t nd2 = q_map[command];
      //std::cout<<" "<<nd2<<command<<std::endl;
      addDirEdge(nd1, nd2);
    }
    computeDistance();
    computePath();
  }
}		
void Graph::readLayout2(string file){  //Reading physical layout of undirected qubit coupling
  //std::cout<<"Hello "<<std::endl;
  char fileName[50];
  int i, ln = file.length();
  for(i=0; i < ln; i++)
    fileName[i] = file[i];
  fileName[i] = '\0';
  string line;
  ifstream myfile(fileName);
  std::map<string, node_t> q_map;
  if(myfile.is_open()){
    while(getline(myfile,line)){ 
      if(0 == line.size() || 1 == line.size()) continue; //skip empty lines
      if('#' == line[0])  continue; //skip comments  
      istringstream iss(line);
      string command;
      iss >> command;
      if("Qubit:" == command){ //Physical Qubits
        while(iss >> command) {
		q_map[command] = nodes.size();
       		addNode(nodes.size());
	}
        continue;
      }
      if("end" == command) break; //End of Layout 	
      if("Coupling:" == command) continue;//Qubits Coupling
      //read all coupling map information
      node_t nd1 = q_map[command];
      //std::cout<<nd1<<command;
      iss >> command;
      node_t nd2 = q_map[command];
      //std::cout<<" "<<nd2<<command<<std::endl;
      addEdge(nd1, nd2);
    }
    computeDistance2();
    computePath();
  }
}
/*
void Graph::mappLayout(Graph p) {
  std::vector<std::vector<map_t> > m = std::vector<std::vector<map_t> > ();
  std::map<node_t, bool > col = std::map<node_t, bool > ( ); //column mapped
  for(std::vector<node_t>::iterator j=p.nodes.begin(); j!=p.nodes.end(); ++j) col[*j] = false;
  int totcost = 0;
  for(std::vector<node_t>::iterator i=this->nodes.begin(); i!=this->nodes.end(); ++i) {
    int cost = this->edges[*i].size() ;
    node_t colp;
    std::vector<map_t> mi = std::vector<map_t> ();
    for(std::vector<node_t>::iterator j=p.nodes.begin(); j!=p.nodes.end(); ++j) {
      mi.push_back(INVALID);
      int cost1 = this->edges[*i].size() - p.edges[*j].size();//this->edges[*i].size() + this->redges[*i].size() - p.edges[*j].size() - p.redges[*j].size();
      //std::cout<<"Size: "<<this->edges[*i].size()<<" "<<this->redges[*i].size()<<" "<<p.edges[*j].size()<<" "<<p.redges[*j].size()<<" "<<cost1<<std::endl;
      cost1 = abs(cost1);
      if ( col[*j] == false && cost > cost1 ) {
           cost = cost1;
           colp = *j;
      }
    }
    col[colp] = true;
    mi[colp] = cost;
    m.push_back(mi);
    totcost += cost;
  }
  std::vector<std::vector<map_t> > m1 = m;

  std::cout<<"Initial mapping:"<<std::endl;
  for(std::vector<node_t>::iterator i=this->nodes.begin(); i!=this->nodes.end(); ++i) {
    for(std::vector<node_t>::iterator j=p.nodes.begin(); j!=p.nodes.end(); ++j) 
      std::cout<<setw(3)<<m[*i][*j]<<" ";
    std::cout<<"Total cost"<<totcost<<std::endl;
  }
}

void Graph::prune( std::vector<std::vector<map_t> >& m, Graph p){
  bool changed;
  do{
    changed = false;
    for(std::vector<node_t>::iterator i=this->nodes.begin(); i!=this->nodes.end(); ++i) {
      for(std::vector<node_t>::iterator j=p.nodes.begin(); j!=p.nodes.end(); ++j){
        if (m[*i][*j] != INVALID){
           for(std::vector<node_t>::iterator j=this->edges[*i].begin(); j!=this->edges[*i].end(); ++j)  
        }
      }
    }
  }while (changed);   
}*/
//physical qubits 0, 1 8
//logical qubits 0, 1, 2
//phmap[2] = 8
//Mapping circuit qubits to specific physical qubits
void Graph::mapQubits(Graph & pg, std::map<node_t, node_t> & phmap, std::vector<node_t> & pqubits){  
  std::map<node_t, distance_t> cdistance = std::map<node_t, distance_t> (); //cumulative distance of a node from other nodes
  std::vector<node_t> pnodes = std::vector<node_t> (); //all ordered physical nodes
  for(std::vector<node_t>::iterator i=pqubits.begin(); i!=pqubits.end(); ++i){
    cdistance[*i] = 0; //initialize distance
    for(std::vector<node_t>::iterator j=pqubits.begin(); j!=pqubits.end(); ++j){
      cdistance[*i] += pg.getDistance(*i, *j);
    }
    //std::cout<<"Node "<<*i<<" Distance "<<cdistance[*i]<<std::endl;
    if (pnodes.empty()) pnodes.push_back(*i); 
    else {
      std::vector<node_t>::iterator j=pnodes.begin();
      for(; j!=pnodes.end(); ++j){
        if(cdistance[*i] < cdistance[*j]) {
          pnodes.insert(j, *i);
          break;
        } 
      }
      if (j == pnodes.end()) pnodes.push_back(*i);
    }    
  }
  /*for(std::vector<node_t>::iterator j=pnodes.begin(); j!=pnodes.end(); ++j){
      std::cout<<"node"<<*j<<std::endl;        
  }*/
  std::vector<node_t>::iterator j=pnodes.begin();  //pick the logical and physical node     
  for(std::vector<node_t>::iterator i=this->nodes.begin(); 
                 i!=this->nodes.end() && j!=pqubits.end(); ++i, ++j){
    phmap[*i] = *j; //mapp the logical qubit *i to physical qubit *j
  }
}

//Mapping Hexagonal grid layout with specified width and height
//in x direction nodes are placed with distance 2, i.e. i, i+2, i+4
//with i start from 0 or 1 depending on y is odd or even.
void Graph::mapLayout(std::map<node_t, node_t> & phmap, std::pair<int, int> dim){
  //number of qubits in layout
  int tot_qubits = dim.first * dim.second;
  
  //coordinate of centered node
  int y = dim.second / 2;             
  int x = (1 - (y % 2 )) + 2 * (dim.first / 2);   
  
  //coordinate of selected node for mapping
  int x_new = x, y_new = y;
  int r = 0; //radious increase by 1
  int d = 1; //nodes that can be mapped at radious r d = d + r*6
  int n = 0; //nodes mapped
  int k = 0, l = 0;
  //std::vector<node_t> temp {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41};
  bool flag = false;
  for(const auto& lq: nodes){  // temp pick each logical qubit
    //std::cout<<"x:"<<x_new<<", y:"<<y_new<<std::endl;
    while (true) {
      //std::cout<<"K: "<<k<<", L:"<<l<<std::endl;
      x_new = x + k, y_new = y + l; //set the coordinate of next qubit to be mapped    
      node_t pq =  (y_new * dim.first) + (x_new / 2); //finding physical qubit at (x_new, y_new)      
      n += 1; //next node 
      //Compute coordinate of next node
      if ( n >= d){
        r += 1;  //increase radious
        d += r * 6; // in hexagonal architecture 6 neighbouring nodes        
        k = (-1) * r * 2; //-i;
        l = 0; //j;
        //std::cout<<"k:"<<k<<", l:"<<l<<std::endl;
      } 
      else if ( flag == false && l != 0) {
        //k = -k;
        l = -l;            
        flag = true;
      }
      else if (flag == true  &&  abs(k) <= r ){      
        if (k == r && l == -r) k += 1, l += 1;  //last coordinate having abs(k) =
        else k += 2, l = -l;
        flag = false;
      }
      else{
       k += 1;
       l = (k > r) ? abs(l) - 1 : abs(l) + 1;
       flag = false;
      }
      if ( pq < tot_qubits && x_new >= 0 && x_new < 2 * dim.first && y_new >= 0 && y_new < dim.second) { // Valid node from the layout 
        phmap[lq] = pq;       
        //std::cout<<"x:"<<x_new<<", y:"<<y_new<<", Q: "<<pq<<std::endl;
        break;
      }
    }        
  }
}

//Layout mapping based on degree of verticess
void Graph::mapLayout(Graph pg, std::map<node_t, node_t> & phmap){
  std::vector<node_t> pnodes = std::vector<node_t> (); //all ordered physical nodes
  std::map<node_t, bool> state = std::map<node_t, bool> (); //mapped status

  //ordering physical qubits based on their degree of association 
  for(std::vector<node_t>::iterator i=pg.nodes.begin(); i!=pg.nodes.end(); ++i){
    bool inserted = false;
    for(std::vector<node_t>::iterator j=pnodes.begin(); j!=pnodes.end(); ++j){
      if(pg.edges[*i].size() > pg.edges[*j].size()) {
        pnodes.insert(j, *i);
        inserted = true;
        break;
      }
      else if (pg.edges[*i].size() == pg.edges[*j].size() && pg.redges[*i].size() > pg.redges[*j].size()) {
        pnodes.insert(j, *i);
        inserted = true;
        break;
      }
    }
    if (!inserted) pnodes.push_back(*i); 
    state[*i] = false;
  }
  //ordering circuit qubits based on their corresponding degree of interactions
  std::vector<node_t> cnodes = std::vector<node_t> (); //all ordered circuit nodes
  for(std::vector<node_t>::iterator i=this->nodes.begin(); i!=this->nodes.end(); ++i){
    bool inserted = false;
    for(std::vector<node_t>::iterator j=cnodes.begin(); j!=cnodes.end(); ++j){
      if(this->edges[*i].size() > this->edges[*j].size()) {
        cnodes.insert(j, *i);
        inserted = true;
        break;
      }
      else if (this->edges[*i].size() == this->edges[*j].size() && this->redges[*i].size() > this->redges[*j].size()) {
        cnodes.insert(j, *i);
        inserted = true;
        break;
      }
    }
    if (!inserted) cnodes.push_back(*i); 
  }
  /*std::cout<<"Ordered Physical Qubits:"<<std::endl;
  for(std::vector<node_t>::iterator j=pnodes.begin(); j!=pnodes.end(); ++j){
    std::cout<<(*j)<<" ";
  }
  std::cout<<std::endl;
  std::cout<<"Ordered Logical Qubits:"<<std::endl;
  for(std::vector<node_t>::iterator j=cnodes.begin(); j!=cnodes.end(); ++j){
    std::cout<<(*j)<<" ";
  }
  std::cout<<std::endl;*/

  //std::map<node_t, node_t> phmap = std::map<node_t, node_t> ( ); //logical node to physical node mapping
  std::vector<node_t> stk = std::vector<node_t> ( );
  for(std::vector<node_t>::iterator i=cnodes.begin(); i!=cnodes.end();){  //pick the logical node
    if(phmap.end() == phmap.find(*i)){ //the logical node still not mapped                 
      for(std::vector<node_t>::iterator j=pnodes.begin(); j!=pnodes.end(); ++j){ //pick one physical node
        if (i != cnodes.end() && state[*j] == false) { //physical node not already mapped
          phmap[*i] = *j; //mapp the logical qubit *i to physical qubit *j
          state[*j] = true;
          //std::cout<<"Physical Mapping2:"<<*i<<" "<<*j<<std::endl;         
          //for(std::vector<node_t>::iterator k=this->edges[*i].begin(); k!=this->edges[*i].end(); ++k){
          //std::cout<<"Physical Mapping1:"<<*i<<" "<<*j<<std::endl;
          stk.push_back(*j); // push the mapped physical qubit into stack 
          i++; // get the next logical qubit 
          while (i != cnodes.end() && !stk.empty()){
            //node_t nd = stk.front();
            //stk.erase(stk.begin());
            node_t phn;
            int c_cost = 0;
            int p_cost = 0;
            bool unassigned = false;
            for(std::vector<node_t>::iterator t=stk.begin(); t!=stk.end(); ++t){ //Explore the neighbor of physical node *t
              //std::cout<<"Helloooo:"<<*t<<std::endl;
              node_t nd = *t;
              unassigned = false;
              //std::cout<<"Explore the neighbor of Node:"<<nd<<std::endl; 
              for(std::vector<node_t>::iterator l=pg.edges[nd].begin(); l!=pg.edges[nd].end(); ++l){ //Explore forward neighbor
                if (state[*l] == false) { // unmapped physical qubit *l
                  //std::cout<<"state of "<<*l<<" is "<< state[*l]<<std::endl;
                  unassigned = true;
                  c_cost = 0;
                  for(std::vector<node_t>::iterator k=cnodes.begin(); k!=cnodes.end(); ++k) // Get all mapped logical qubits
                    if(phmap.end() != phmap.find(*k)) c_cost += roundNum(pg.distance[*l][phmap[*k]]); // Compute the distance from *l
                  if (p_cost ==0 || p_cost > c_cost) {
                    phn = *l;
                    p_cost = c_cost;
                  }
                }
              }  
              for(std::vector<node_t>::iterator m=pg.redges[nd].begin(); m!=pg.redges[nd].end(); ++m){ //Explore backward neighbor
              //std::cout<<"I_node:"<<*i<<" P_node1:"<<*l<<" cost"<<cost<<" rcost"<<rcost<<std::endl;
                if (state[*m] == false) {  // unmapped physical qubit *m    
                  //std::cout<<"state of "<<*m<<" is "<< state[*m]<<std::endl;
                  unassigned = true;
                  c_cost = 0;
                  for(std::vector<node_t>::iterator k=cnodes.begin(); k!=cnodes.end(); ++k) // Get all mapped logical qubits
                    if(phmap.end() != phmap.find(*k)) c_cost += roundNum(pg.distance[*m][phmap[*k]]);  // Compute the distance from *l
                  if (p_cost ==0 || p_cost > c_cost) {
                    //std::cout<<"I_node:"<<*i<<" P_node1:"<<*l<<" P_node2:"<<*m<<" cost"<<cost<<" rcost"<<rcost<<std::endl;
                    phn = *m;
                    p_cost = c_cost;
                  }
                }
              }
              if (!unassigned) {
                unsigned pos = std::distance(stk.begin(), t);
                stk.erase(t);
                //std::cout<<"Hello:"<<*t<<std::endl;
                t = stk.begin() + pos - 1;
                //std::cout<<"Hello:"<<*t<<std::endl;
                //for (std::vector<node_t>::iterator t2=stk.begin(); t2!=stk.end(); ++t2)
                  //std::cout<<" "<<*t2;
                //std::cout<<std::endl;
              }
            }
            if(unassigned){
              phmap[*i] = phn;
              state[phn] = true; 
              //std::cout<<"state of "<<phn<<" is "<< state[phn]<<" cost "<<p_cost<<std::endl;
              stk.insert(stk.end(), phn);
              //for (std::vector<node_t>::iterator t2=t; t2!=stk.end(); ++t2)
                 //std::cout<<" "<<*t2;
              //std::cout<<std::endl;
              i++;
            } 
          } 
        } 
      }
    }
  }
}
  /*std::cout<<"Physical Mapping:"<<std::endl;
  for(std::vector<node_t>::iterator j=cnodes.begin(); j!=cnodes.end(); ++j){
    std::cout<<(*j)<<" "<<phmap[*j]<<std::endl;
  }
  std::cout<<std::endl;*/


