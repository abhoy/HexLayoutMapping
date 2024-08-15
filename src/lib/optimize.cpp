#include "optimize.hpp"


void getPath(grph::node_t p1, grph::node_t p2, std::vector<grph::node_t> & path, grph::Graph & pg){
  //std::cout<<"Hello"<<pg.getNodeInPath(p1, p2)<<std::endl;
  if ( pg.getNodeInPath(p1, p2) !=  p2){
    getPath(p1, pg.getNodeInPath(p1, p2), path, pg);
    getPath(pg.getNodeInPath(p1, p2), p2, path, pg);
  }
  else path.push_back(p2);
}
//ibm qx architecture
void getAllPaths(grph::node_t p1, grph::node_t p2, std::vector<std::vector<grph::node_t> > & paths, grph::Graph & pg, int i){
  grph::node_t px = pg.getNodeInPath(p1, p2); 
  if (px !=  p2){ //pg.getNodeInPath(p1, p2)
    //std::cout<<"1.["<<p1<<","<<p2<<"]="<<pg.getNodeInPath(p1, p2)<<std::endl;
    //std::cout<<i<<" Node Between "<<p1<<" "<<p2<<"-"<<pg.getNodeInPath(p1, p2)<<std::endl;
    paths[i].push_back(p1);
    getAllPaths(px, p2, paths, pg, i);

    //std::cout<<pg.getDistance(pg.getNodeInPath(p1, p2), p2)<<" : "<<pg.getDistance(pg.getNodeInPath2(p1, p2), p2)<<std::endl;
    grph::node_t py = pg.getNodeInPath2(p1, p2);
    bool isParent = false;
    for (int i = 0; i < paths.size(); i++){
      for (int j = 0; j < paths[i].size(); j++){
        if (paths[i][j] == py) isParent = true;
        break;
      }
    }
    if ( py !=  -1 && isParent == false ) { 
      //std::cout<<"2.["<<p1<<","<<p2<<"]="<<pg.getNodeInPath2(p1, p2)<<std::endl;
      paths.push_back(std::vector<grph::node_t> ());
      int k = paths.size()-1;
      //if (k > 0 && std::find(paths[k-1].begin(), paths[k-1].end(), p1) != paths[k-1].end()) {
      for(std::vector<grph::node_t>::iterator p = paths[i].begin(); p != paths[i].end(); ++p) {
        //std::cout<<" "<<*p;
        paths[k].push_back(*p);
        if ( *p == p1) break;
      }
      //std::cout<<paths.size()<<std::endl;
      //paths[k].push_back(p1);
      getAllPaths(py, p2, paths, pg, k);
    }
  }
  else {
    paths[i].push_back(p1);
    paths[i].push_back(p2);
  }
}
//hexagonal architecture
// path from p1 to p2
int getAllPaths(grph::node_t p1, grph::node_t p2, std::vector<std::vector<grph::node_t> > & paths, 
                                                                     std::pair<int, int> dim, int i){
  if (grph::getDistance(p1, p2, dim) == 0) { //std::cout<<std::endl; 
    paths[i].push_back(p2);
    return i + 1;
  }
  else {
    if (paths[0].size() == 0) paths[0].push_back(p1);
     // coordinate of node p1
    int y_p1 = p1 / dim.first;
    int x_p1 = (1 - (y_p1 % 2 )) + 2 * (p1 % dim.first);

    // coordinate of node p2
    int y_p2 = p2 / dim.first;
    int x_p2 = (1 - (y_p2 % 2 )) + 2 * (p2 % dim.first);

    int x_diff = abs(x_p1 - x_p2);
    int y_diff = abs(y_p1 - y_p2);
  
    int x_sign = x_p1 > x_p2 ? -1 : 1;
    int y_sign = y_p1 > y_p2 ? -1 : 1;  

    int x_max = 2 * dim.first - 1;
    int y_max = dim.second - 1;
    int x_p1_new = x_p1 + x_sign * 2; //moving in x direction
    grph::node_t p1_x_new =  (y_p1 * dim.first) + (x_p1_new / 2);
    if ( x_p1_new >= 0 && x_p1_new <= x_max){      //boundary checking of new node
      //if (x_p1_new <= x_p2 + x_sign * 1){
      if (grph::getDistance(p1_x_new, p2, dim) < grph::getDistance(p1, p2, dim)){
        //std::cout<<"X Next Node: "<<p1_x_new<<"("<<i<<")";//<<std::endl
        paths[i].push_back(p1_x_new);
        i = getAllPaths(p1_x_new, p2, paths, dim, i);        
      }
    }
    int y_p1_new = y_p1 + y_sign * 1; //moving in y direction
    int x_p1_new_y = x_p1 + x_sign * 1; // decreasing x distnce
    grph::node_t p1_y_new =  (y_p1_new * dim.first) + (x_p1_new_y / 2);    

    if ( x_p1_new_y >= 0 && x_p1_new_y <= x_max && y_p1_new >= 0 && y_p1_new <= y_max ){ //boundary checking of new node
      if (grph::getDistance(p1_y_new, p2, dim) < grph::getDistance(p1, p2, dim) ){
        //if it is an alternate path copy the previous path nodes
        if (grph::getDistance(p1_y_new, p2, dim) == grph::getDistance(p1_x_new, p2, dim)) {
          paths.push_back(std::vector<grph::node_t> ());
          i = paths.size() - 1;
          for (auto node : paths[i-1]) {
            if (node == p1_x_new) break;
            paths[i].push_back(node);
          }
        }
      
        //std::cout<<"Y d X Next Node: "<<p1_y_new<<"("<<i<<")";//<<std::endl
        paths[i].push_back(p1_y_new);
        i = getAllPaths(p1_y_new, p2, paths, dim, i);
      }
    }
    int x_p1_new_y2 = x_p1 + (-1) * x_sign * 1; //moving in y direction increasing x distnce
    grph::node_t p1_y2_new =  (y_p1_new * dim.first) + (x_p1_new_y2 / 2);
    if ( x_p1_new_y2 >= 0 && x_p1_new_y2 <= x_max && y_p1_new >= 0 && y_p1_new <= y_max ){ //boundary checking of new node       
      if (grph::getDistance(p1_y2_new, p2, dim) < grph::getDistance(p1, p2, dim) ){
        //if it is an alternate path copy the previous path nodes
        if (grph::getDistance(p1_y_new, p2, dim) == grph::getDistance(p1_y2_new, p2, dim)) {
          paths.push_back(std::vector<grph::node_t> ());
          i = paths.size() - 1;
          for (auto node : paths[i-1]) {
            if (node == p1_y_new) break;
            paths[i].push_back(node);
          }
        }      
        //std::cout<<"Y i X Next Node: "<<p1_y2_new<<"("<<i<<")";//<<std::endl
        paths[i].push_back(p1_y2_new);
        i = getAllPaths(p1_y2_new, p2, paths, dim, i);
      }
    }
    /*std::cout<<p1<<"("<<x_p1<<", "<<y_p1<<") "<<p2<<"("<<x_p2<<", "<<y_p2<<") "<<std::endl;
    std::cout<<p1_x_new<<"("<<x_p1_new<<", "<<y_p1<<") "<<std::endl;
    std::cout<<p1_y_new<<"("<<x_p1_new_y<<", "<<y_p1_new<<") ";//<<std::endl*/
  }
}

void removeDuplicatePath(std::vector<std::vector<grph::node_t> > & paths){
  //displayAllPaths(paths);
  for(std::vector<std::vector<grph::node_t> >::iterator p1 = paths.begin(); p1 != paths.end() - 1; ++p1){
    //std::cout<<"hello"<<paths.size()<<std::endl;
    unsigned pos1 = std::distance(paths.begin(),p1);
    for(std::vector<std::vector<grph::node_t> >::iterator p2 = p1 + 1; p2 != paths.end(); ++p2){
      //std::cout<<"world"<<pos1<<std::endl;
      bool identical = true;
      for(std::vector<grph::node_t>::iterator q1 = (*p1).begin(), q2 = (*p2).begin(); q1 != (*p1).end(); ++q1, ++q2){
        if (*q1 != *q2) identical = false;
      }
      if (identical) {
        unsigned pos2 = std::distance(paths.begin(),p2);
        paths.erase(p2);
        p2 = paths.begin() + pos2 - 1;
        if (p2 == paths.end()) break;
        //std::cout<<"world"<<pos2<<paths.size()<<std::endl;
      }
    }
    p1 = paths.begin() + pos1;
    if (p1 + 1 == paths.end()) break;
  }
  //displayAllPaths(paths);
}

void addInitNode(std::vector<std::vector<grph::node_t> > & paths, grph::node_t nd){
  for ( int i = 0; i < paths.size(); i++) paths[i].insert(paths[i].begin(), nd);
}

bool isInvalidPath(std::vector<grph::node_t> & path, std::map<rev::line_t, rev::line_t> & plmap){
  bool invalidPath = false;
  for (const auto & p : path){
    //std::cout<<*p<<" ";
    if (plmap.find(p) == plmap.end()) {
      invalidPath = true; 
      break;
    }
  }
  //std::cout<<std::endl;
  return invalidPath;
}

void displayAllPaths(std::vector<std::vector<grph::node_t> > & paths){
  for(std::vector<std::vector<grph::node_t> >::iterator p = paths.begin(); p != paths.end(); ++p){
    std::cout<<"Path: ";
    for(std::vector<grph::node_t>::iterator q = (*p).begin(); q != (*p).end(); ++q)
      std::cout<<*q<<" ";
    std::cout<<std::endl;
  }
}

void displayPath(std::vector<grph::node_t>& path){
  std::cout<<"Path: ";
  for(std::vector<grph::node_t>::iterator p = path.begin(); p != path.end(); ++p)
      std::cout<<*p<<" ";
  std::cout<<std::endl;
}


//IBM QX Mapping
mincost_t minCostPathIndex(std::vector<grph::node_t>& path, std::map<grph::node_t, grph::node_t> & temp_phmap, grph::Graph & pg, grph::Graph & temp_g){
  std::vector<grph::node_t> nodes = temp_g.getNodes();
  mincost_t mincost;
  bool empty = true;
  for(std::vector<grph::node_t>::iterator p = path.begin(); p != path.end()-1; ++p){
    grph::node_t control = *p;
    grph::node_t target = *(p+1);

    std::map<grph::node_t, grph::node_t> temp_phmap2 = std::map<grph::node_t, grph::node_t> ();
    std::map<grph::node_t, grph::node_t> temp_phmap3 = std::map<grph::node_t, grph::node_t> ();
    for(std::vector<grph::node_t>::iterator n = nodes.begin(); n != nodes.end(); ++n){
      temp_phmap2[*n] = temp_phmap[*n];
      temp_phmap3[temp_phmap[*n]] = *n;
    }
                 
    double cur_swap_cost = 0;
    //forward move in path
    for(std::vector<grph::node_t>::iterator pf = path.begin(); pf != p; ++pf){                  
      //perform swaping in temporary location
      cur_swap_cost += 1;
                   
      grph::node_t temp = temp_phmap2[temp_phmap3[*pf]];
      temp_phmap2[temp_phmap3[*pf]] = temp_phmap2[temp_phmap3[*(pf+1)]];
      temp_phmap2[temp_phmap3[*(pf+1)]] = temp;

      temp = temp_phmap3[*pf];
      temp_phmap3[*pf] = temp_phmap3[*(pf+1)];
      temp_phmap3[*(pf+1)] = temp;
    }
    //backward move in path
    for(std::vector<grph::node_t>::reverse_iterator pb = path.rbegin(),  end(p + 2); pb != end; ++pb){
      //perform swaping in temporary location
      cur_swap_cost += 1;
                   
      grph::node_t temp = temp_phmap2[temp_phmap3[*pb]];
      temp_phmap2[temp_phmap3[*pb]] = temp_phmap2[temp_phmap3[*(pb+1)]];
      temp_phmap2[temp_phmap3[*(pb+1)]] = temp;

      temp = temp_phmap3[*pb];
      temp_phmap3[*pb] = temp_phmap3[*(pb+1)];
      temp_phmap3[*(pb+1)] = temp;
    }
    
    //current control line is not valid for physical coupling map, i.e *(p+1)->*p  
    if (pg.getDistance(*p,*(p+1)) == RCOST){
      //perform swaping in temporary location
      cur_swap_cost += RCOST;

      grph::node_t temp = temp_phmap2[temp_phmap3[*p]];
      temp_phmap2[temp_phmap3[*p]] = temp_phmap2[temp_phmap3[*(p+1)]];
      temp_phmap2[temp_phmap3[*(p+1)]] = temp;

      temp = temp_phmap3[*p];
      temp_phmap3[*p] = temp_phmap3[*(p+1)];
      temp_phmap3[*(p+1)] = temp;
    }
    double coupling_cost = mappingCost(temp_phmap2, pg, temp_g);
    if (empty || coupling_cost < mincost.swap_cost){
      mincost.control = control;
      mincost.target = target;
      mincost.swaps = cur_swap_cost;
      mincost.swap_cost = coupling_cost;
      empty = false;
    }
  }
  return mincost;
}

//Hexagonal Mapping
mincost_t minCostPathIndex(std::vector<grph::node_t>& path, std::map<grph::node_t, grph::node_t> & phmap, std::pair<int, int> dim, grph::Graph & cg){
  std::vector<grph::node_t> nodes = cg.getNodes();
  mincost_t mincost;
  bool empty = true;
  //for(std::vector<grph::node_t>::iterator p = path.begin(); p != path.end()-1; ++p){
  for(int i = 0; i < path.size() - 1; i++){
    grph::node_t control = path[i];
    grph::node_t target = path[i + 1];

    std::map<grph::node_t, grph::node_t> temp_phmap = std::map<grph::node_t, grph::node_t> ();
    std::map<grph::node_t, grph::node_t> temp_phmap2 = std::map<grph::node_t, grph::node_t> ();
    for(const auto & n : nodes){
      temp_phmap[n] = phmap[n]; //logical to physical mapping
      temp_phmap2[phmap[n]] = n; //physical to logical mapping
    }
                 
    double cur_swap_cost = 0;
    //forward move in path
    for(int j = 0; j < i; j++){  //physical qubit from path                
      //perform swaping in temporary location
      cur_swap_cost += 1;
                   
      grph::node_t temp = temp_phmap[temp_phmap2[path[j]]]; //update mapping of logical to physical qubit
      temp_phmap[temp_phmap2[path[j]]] = temp_phmap[temp_phmap2[path[j + 1]]];
      temp_phmap[temp_phmap2[path[j + 1]]] = temp;

      temp = temp_phmap2[path[j]];      //update mapping of physical to logical qubit
      temp_phmap2[path[j]] = temp_phmap2[path[j + 1]];
      temp_phmap2[path[j + 1]] = temp;
    }
    //backward move in path
    for(int k = path.size() - 1; k < i + 1; k--){
      //perform swaping in temporary location
      cur_swap_cost += 1;
                   
      grph::node_t temp = temp_phmap[temp_phmap2[path[k]]]; //update mapping of logical to physical qubit
      temp_phmap[temp_phmap2[path[k - 1]]] = temp_phmap[temp_phmap2[path[k - 1]]];
      temp_phmap[temp_phmap2[path[k - 1]]] = temp;

      temp = temp_phmap2[path[k]];          //update mapping of physical to logical qubit
      temp_phmap2[path[k]] = temp_phmap2[path[k - 1]];
      temp_phmap2[path[k - 1]] = temp;
    }
    
    
    double coupling_cost = mappingCost(temp_phmap2, dim, cg);
    if (empty || coupling_cost < mincost.swap_cost){
      mincost.control = control;
      mincost.target = target;
      mincost.swaps = cur_swap_cost;
      mincost.swap_cost = coupling_cost;
      empty = false;
    }
  }
  return mincost;
}
//IBM QX
void updateMapping(std::vector<grph::node_t>& path, std::map<grph::node_t, grph::node_t> & temp_phmap, grph::Graph& pg, 
                             std::vector<grph::node_t> & nodes, grph::node_t control, grph::node_t target){

  std::map<grph::node_t, grph::node_t> temp_phmap2 = std::map<grph::node_t, grph::node_t> ();
  for(std::vector<grph::node_t>::iterator n = nodes.begin(); n != nodes.end(); ++n){
    temp_phmap2[temp_phmap[*n]] = *n;
  } 
  //std::cout<<"Hello12-6-1 "<<control<<" "<<target<<std::endl;
  //displayPath(path);             
  //forward move in path
  for(std::vector<grph::node_t>::iterator pf = path.begin(); *pf != control; ++pf){                  
    //std::cout<<"Hello12-6-1-1"<<std::endl;
    grph::node_t temp = temp_phmap[temp_phmap2[*pf]];
    temp_phmap[temp_phmap2[*pf]] = temp_phmap[temp_phmap2[*(pf+1)]];
    temp_phmap[temp_phmap2[*(pf+1)]] = temp;
    //std::cout<<"Hello12-6-1-2"<<std::endl;
    temp = temp_phmap2[*pf];
    temp_phmap2[*pf] = temp_phmap2[*(pf+1)];
    temp_phmap2[*(pf+1)] = temp;
    //std::cout<<"Hello12-6-1-3"<<std::endl;
  }
  //std::cout<<"Hello12-6-2"<<std::endl;
  //backward move in path
  for(std::vector<grph::node_t>::reverse_iterator pb = path.rbegin(); *pb != target; ++pb){
    
    grph::node_t temp = temp_phmap[temp_phmap2[*pb]];
    temp_phmap[temp_phmap2[*pb]] = temp_phmap[temp_phmap2[*(pb+1)]];
    temp_phmap[temp_phmap2[*(pb+1)]] = temp;

    temp = temp_phmap2[*pb];
    temp_phmap2[*pb] = temp_phmap2[*(pb+1)];
    temp_phmap2[*(pb+1)] = temp;
  }
  //std::cout<<"Hello12-6-3"<<std::endl;  
  //current control line is not valid for physical coupling map, i.e *(p+1)->*p  
  if (pg.getDistance(control,target) == RCOST){

    grph::node_t temp = temp_phmap[temp_phmap2[control]];
    temp_phmap[temp_phmap2[control]] = temp_phmap[temp_phmap2[target]];
    temp_phmap[temp_phmap2[target]] = temp;

    temp = temp_phmap2[control];
    temp_phmap2[control] = temp_phmap2[target];
    temp_phmap2[target] = temp;
  }
  //std::cout<<"Hello12-6-4"<<std::endl;
}

//Hexagonal Architecture
void updateMapping(std::vector<grph::node_t>& path, std::map<grph::node_t, grph::node_t> & phmap, std::pair<int, int> dim, 
                             std::vector<grph::node_t> & nodes, grph::node_t control, grph::node_t target){

  std::map<grph::node_t, grph::node_t> temp_phmap = std::map<grph::node_t, grph::node_t> ();
  for(const auto & n : nodes){
    temp_phmap[phmap[n]] = n; //current physical to logical mapping
  } 
  
  for(int i = 0; path[i] != control; i++){
    //std::cout<<"Hello12-6-1-1"<<std::endl;
    grph::node_t temp = phmap[temp_phmap[path[i]]];
    phmap[temp_phmap[path[i]]] = phmap[temp_phmap[path[i + 1]]];
    phmap[temp_phmap[path[i + 1]]] = temp;
    //std::cout<<"Hello12-6-1-2"<<std::endl;
    temp = temp_phmap[path[i]];
    temp_phmap[path[i]] = temp_phmap[path[i + 1]];
    temp_phmap[path[i + 1]] = temp;
    //std::cout<<"Hello12-6-1-3"<<std::endl;
  }
  //std::cout<<"Hello12-6-2"<<std::endl;
  //backward move in path
  for(int j = path.size() - 1; path[j] != target; j--){ 
    grph::node_t temp = phmap[temp_phmap[path[j]]];
    phmap[temp_phmap[path[j]]] = phmap[temp_phmap[path[j - 1]]];
    phmap[temp_phmap[path[j - 1]]] = temp;

    temp = temp_phmap[path[j]];
    temp_phmap[path[j]] = temp_phmap[path[j - 1]];
    temp_phmap[path[j - 1]] = temp;
  }  
}

//IBM QX
//computing coupling cost
double swapCost(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& cg, int base, unsigned depth){
  double swap_cost=0;
  grph::Graph temp_g = cg;
  std::vector<grph::node_t> nodes = temp_g.getNodes();
  std::map<grph::node_t, std::vector<grph::node_t> > edge_list = temp_g.getEdges();
  std::map<grph::node_t, bool> rev_map = std::map<grph::node_t, bool> ();
  std::map<grph::node_t, grph::node_t> temp_phmap = std::map<grph::node_t, grph::node_t> ();
  std::map<grph::node_t, grph::node_t> temp_phmap2 = std::map<grph::node_t, grph::node_t> ();
  
  for(std::vector<grph::node_t>::iterator n = nodes.begin(); n != nodes.end(); ++n){
    temp_phmap[*n] = phmap[*n];
    temp_phmap2[temp_phmap[*n]] = *n;
    rev_map[phmap[*n]] = false;
  }
  //std::cout<<"Hello12-2"<<std::endl;
  for (unsigned d = depth - 1; d >= 0; d--){
    double cur_weight = std::pow(base, d);
    bool edge_present = false;
    for(std::vector<grph::node_t>::iterator n = nodes.begin(); n != nodes.end(); ++n){
      for(std::vector<grph::node_t>::iterator m = edge_list[*n].begin(); m != edge_list[*n].end(); ++m){
        if (temp_g.getWeight(*n, *m) > 0) edge_present = true;
        //std::cout<<"Hello12-3"<<std::endl;
        if (temp_g.getWeight(*n, *m) >= cur_weight) {
          temp_g.updateEdgeWeight(*n, *m, cur_weight);
          //std::cout<<"Hello12-4"<<std::endl;
          if (pg.getDistance(temp_phmap[*n], temp_phmap[*m]) > 0) {
             //std::cout<<"Hello12-5x"<<std::endl;
             std::vector<std::vector<grph::node_t> > paths = std::vector<std::vector<grph::node_t> > ();
             paths.push_back(std::vector<grph::node_t> ());
             //std::cout<<"Hello12-5z: "<<temp_phmap[*n]<<" : "<<temp_phmap[*m]<<std::endl;
             getAllPaths(temp_phmap[*n],temp_phmap[*m], paths, pg, 0); //0, 3
             //std::cout<<"Hello12-5y"<<std::endl;
             removeDuplicatePath(paths);
             //addInitNode(paths, phmap[*n]);
             //displayAllPaths(paths);
             //exit(0);
             //std::cout<<"Hello12-5"<<std::endl;
             mincost_t mincost;
             bool empty = true;
             int index;
             for (int i=0; i < paths.size(); i++){
               if (isInvalidPath(paths[i], temp_phmap2)) continue;
               mincost_t tempcost = minCostPathIndex(paths[i], temp_phmap, pg, temp_g);
               if (empty || mincost.swap_cost > tempcost.swap_cost || (mincost.swap_cost == tempcost.swap_cost && mincost.swaps > tempcost.swaps)){
                 mincost.control = tempcost.control;
                 mincost.target = tempcost.target;
                 mincost.swaps = tempcost.swaps;
                 mincost.swap_cost = tempcost.swap_cost;
                 empty = false;
                 index = i;
               }            
             }
             //std::cout<<"Hello12-6"<<std::endl;
             updateMapping(paths[index], temp_phmap, pg, nodes, mincost.control, mincost.target);
             //std::cout<<"Hello12-7"<<std::endl;
             swap_cost += mincost.swaps;
             if (pg.getDistance(mincost.control,mincost.target) == RCOST){
               if (rev_map[mincost.control]) swap_cost -= RCOST/2;
               if (rev_map[mincost.target]) swap_cost -= RCOST/2; 
               rev_map[mincost.control] = true;  
               rev_map[mincost.target] = true;  
             }
             else{
               rev_map[mincost.control] = false;  
               rev_map[mincost.target] = false;  
             }
          }
        }
      }
    }
    if (!edge_present) break;
  } 
  return swap_cost;
}
//HEXAGONAL ARCHITECTURE mapping 
//computing coupling cost
double swapCost(std::map<grph::node_t, grph::node_t>& phmap, std::pair<int, int> dim, grph::Graph& cg, int base, unsigned depth){
  double swap_cost=0;
  grph::Graph temp_cg = cg;
  std::vector<grph::node_t> nodes = temp_cg.getNodes();
  std::map<grph::node_t, std::vector<grph::node_t> > edge_list = temp_cg.getEdges();
  //std::map<grph::node_t, bool> rev_map = std::map<grph::node_t, bool> ();
  std::map<grph::node_t, grph::node_t> temp_phmap = std::map<grph::node_t, grph::node_t> ();
  std::map<grph::node_t, grph::node_t> temp_phmap2 = std::map<grph::node_t, grph::node_t> ();
  
  for(const auto& n : nodes){ //logical qubit
    temp_phmap[n] = phmap[n]; //mapping of logical to physical qubit
    temp_phmap2[temp_phmap[n]] = n; //mapping of physical qubit to logical qubit
    //rev_map[phmap[*n]] = false;
  }
  //std::cout<<"Hello12-2"<<std::endl;
  for (unsigned d = depth - 1; d >= 0; d--){
    double cur_weight = std::pow(base, d);
    bool edge_present = false;
    for(const auto & n : nodes){  //logical qubit n
      for(const auto & m : edge_list[n]){ //logical qubit m involved in 2-qubit operation with qubit n
        if (temp_cg.getWeight(n, m) > 0) edge_present = true;
        //std::cout<<"Hello12-3"<<std::endl;
        if (temp_cg.getWeight(n, m) >= cur_weight) {
          temp_cg.updateEdgeWeight(n, m, cur_weight);
          //std::cout<<"Hello12-4"<<std::endl;
          if (grph::getDistance(temp_phmap[n], temp_phmap[m], dim) > 0) {
             //std::cout<<"Hello12-5x"<<std::endl;
             std::vector<std::vector<grph::node_t> > paths;
             paths.push_back(std::vector<grph::node_t> ());
             //std::cout<<"Hello12-5z: "<<temp_phmap[*n]<<" : "<<temp_phmap[*m]<<std::endl;
             getAllPaths(temp_phmap[n],temp_phmap[m], paths, dim, 0); //0, 3
             //std::cout<<"Hello12-5y"<<std::endl;
             //removeDuplicatePath(paths);
             //addInitNode(paths, phmap[*n]);
             //displayAllPaths(paths);
             //exit(0);
             //std::cout<<"Hello12-5"<<std::endl;
             mincost_t mincost;
             bool empty = true;
             int index;
             for (int i=0; i < paths.size(); i++){
               if (isInvalidPath(paths[i], temp_phmap2)) continue;
               mincost_t tempcost = minCostPathIndex(paths[i], temp_phmap, dim, temp_cg);
               if (empty || mincost.swap_cost > tempcost.swap_cost || (mincost.swap_cost == tempcost.swap_cost && mincost.swaps > tempcost.swaps)){
                 mincost.control = tempcost.control;
                 mincost.target = tempcost.target;
                 mincost.swaps = tempcost.swaps;
                 mincost.swap_cost = tempcost.swap_cost;
                 empty = false;
                 index = i;
               }            
             }
             //std::cout<<"Hello12-6"<<std::endl;
             updateMapping(paths[index], temp_phmap, dim, nodes, mincost.control, mincost.target);
             //std::cout<<"Hello12-7"<<std::endl;
             swap_cost += mincost.swaps;
             /*if (pg.getDistance(mincost.control,mincost.target) == RCOST){
               if (rev_map[mincost.control]) swap_cost -= RCOST/2;
               if (rev_map[mincost.target]) swap_cost -= RCOST/2; 
               rev_map[mincost.control] = true;  
               rev_map[mincost.target] = true;  
             }
             else{
               rev_map[mincost.control] = false;  
               rev_map[mincost.target] = false;  
             }*/
          }
        }
      }
    }
    if (!edge_present) break;
  }
  return swap_cost;
}

//IBM QX 
//estimation using coupling cost metric for global ordering
double mappingCost(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& cg){
  double cost=0;
  std::vector<grph::node_t> nodes = cg.getNodes();
  std::map<grph::node_t, std::vector<grph::node_t> > edge_list = cg.getEdges();
  for(const auto & n : nodes){
    for(const auto & m : edge_list[n]){
      cost += pg.getDistance(phmap[n], phmap[m]) * cg.getWeight(n, m);
      //std::cout<<"cost "<<pg.getDistance(phmap[*n], phmap[*m]) * cg.getWeight(*n, *m)<<"cg.getWeight(*n, *m) "<<cg.getWeight(*n, *m)<<std::endl;
    }  
  } 
  
  return cost;
}

//Hexagonal Architecture 
//estimation using coupling cost metric for global ordering
double mappingCost(std::map<grph::node_t, grph::node_t>& phmap, std::pair<int, int> dim, grph::Graph& cg){
  double cost=0;
  std::vector<grph::node_t> nodes = cg.getNodes();
  std::map<grph::node_t, std::vector<grph::node_t> > edge_list = cg.getEdges();
  for(const auto & n : nodes){
    for(const auto & m : edge_list[n]){
      cost += grph::getDistance(phmap[n], phmap[m], dim) * cg.getWeight(n, m);
      //std::cout<<"cost "<<pg.getDistance(phmap[*n], phmap[*m]) * cg.getWeight(*n, *m)<<"cg.getWeight(*n, *m) "<<cg.getWeight(*n, *m)<<std::endl;
    }  
  }   
  return cost;
}

//IBM QX 
//estimation using coupling cost metric for local ordering
double mappingCost(std::map<grph::node_t, grph::node_t>& phmap, std::map<rev::line_t, rev::line_t>& lgmap, 
                                                                                 grph::Graph& pg, grph::Graph& cg){
  double cost=0;
  std::vector<grph::node_t> nodes = cg.getNodes();
  std::map<grph::node_t, std::vector<grph::node_t> > edge_list = cg.getEdges();
  for(const auto & n : nodes){
    for(const auto & m : edge_list[n])
      cost += pg.getDistance(phmap[lgmap[n]], phmap[lgmap[m]]) * cg.getWeight(n, m);
  
  } 
  
  return cost;
}
//Hexagonal Architecture
//estimation using coupling cost metric for local ordering
double mappingCost(std::map<grph::node_t, grph::node_t>& phmap, std::map<rev::line_t, rev::line_t>& lgmap, 
                                                                                 std::pair<int, int> dim, grph::Graph& cg){
  double cost=0;
  std::vector<grph::node_t> nodes = cg.getNodes();
  std::map<grph::node_t, std::vector<grph::node_t> > edge_list = cg.getEdges();
  for(const auto & n : nodes){
    for(const auto & m : edge_list[n])              
      cost += grph::getDistance(phmap[lgmap[n]], phmap[lgmap[m]], dim) * cg.getWeight(n, m);
  
  } 
  
  return cost;
}


//IBM QX 
//estimation using remote CNOT cost metric for global ordering
double remoteCNOTCost(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg, grph::Graph& cg){
  double cost=0;
  std::vector<grph::node_t> nodes = cg.getNodes();
  std::map<grph::node_t, std::vector<grph::node_t> > edge_list = cg.getEdges();
  for(const auto & n : nodes){
    for(const auto & m : edge_list[n]){
      //std::cout<<"ph:"<<phmap[*n]<<" ph:"<<phmap[*m];
      //std::cout<<" lo:"<<*n<<" lo:"<<*m;
      //std::cout<<" cost:"<<4 * pg.getDistance(phmap[*n], phmap[*m]) * cg.getWeight(*n, *m)<<std::endl;
      cost += 4 * pg.getDistance(phmap[n], phmap[m]) * cg.getWeight(n, m);
      //std::cout<<"cost "<<pg.getDistance(phmap[*n], phmap[*m]) * cg.getWeight(*n, *m)<<"cg.getWeight(*n, *m) "<<cg.getWeight(*n, *m)<<std::endl;
    }  
  } 
  
  return cost;
}

//Hexagonal architecture
//estimation using remote CNOT cost metric for global ordering
double remoteCNOTCost(std::map<grph::node_t, grph::node_t>& phmap, std::pair<int, int> dim, grph::Graph& cg){
  double cost=0;
  std::vector<grph::node_t> nodes = cg.getNodes();
  std::map<grph::node_t, std::vector<grph::node_t> > edge_list = cg.getEdges();
  for(const auto & n : nodes){
    for(const auto & m : edge_list[n]){      
      cost += 4 * grph::getDistance(phmap[n], phmap[m], dim) * cg.getWeight(n, m);      
    }  
  }   
  return cost;
}

//IBM QX 
//estimation using remote CNOT cost metric for local ordering
double remoteCNOTCost(std::map<grph::node_t, grph::node_t>& phmap, std::map<rev::line_t, rev::line_t>& lgmap,
                                                                                      grph::Graph& pg, grph::Graph& cg){
  double cost=0;
  std::vector<grph::node_t> nodes = cg.getNodes();
  std::map<grph::node_t, std::vector<grph::node_t> > edge_list = cg.getEdges();
  for(const auto & n : nodes){
    for(const auto & m : edge_list[n]){      
      cost += 4 * pg.getDistance(phmap[lgmap[n]], phmap[lgmap[m]]) * cg.getWeight(n, m);
    }  
  } 
  
  return cost;
}

//Hexagonal architecture
//estimation using remote CNOT cost metric for local ordering
double remoteCNOTCost(std::map<grph::node_t, grph::node_t>& phmap, std::map<rev::line_t, rev::line_t>& lgmap, 
                                                                                 std::pair<int, int> dim, grph::Graph& cg){
  double cost=0;
  std::vector<grph::node_t> nodes = cg.getNodes();
  std::map<grph::node_t, std::vector<grph::node_t> > edge_list = cg.getEdges();
  for(const auto & n : nodes){
    for(const auto & m : edge_list[n])              
      cost += 4 * grph::getDistance(phmap[lgmap[n]], phmap[lgmap[m]], dim) * cg.getWeight(n, m);
  
  } 
  
  return cost;
}

double move_cost(pair<grph::node_t, grph::node_t>& swapped_nodes, grph::Graph& pg){
  //std::cout<<swapped_nodes.first<<", "<<swapped_nodes.second<<std::endl;
  return (pg.getDistance(swapped_nodes.first, swapped_nodes.second) + 1) * 3; // Number of CNOT operations for SWAPPing state
}

/*void optimize(rev::Circuit & ckt){
  std::vector<rev::Gate> gates = ckt.getGates();
  for(std::vector<rev::Gate>::iterator g=gates.begin(); g!=gates.end(); ++g){ 
    for(std::vector<rev::Gate>::reverse_iterator gr(g); gr != gates.rend(); ++gr){
      if (g->getType() != gr->getType() && g->getType() == rev::CX){
        if (g->getControls()[0] == gr->getTargets()[0] || g->getTargets()[0] == gr->getTargets()[0]){ 
          break; 
        }
      }
      else if (g->getType() != gr->getType() && gr->getType() == rev::CX){  
        if (gr->getControls()[0] == g->getTargets()[0] || gr->getTargets()[0] == g->getTargets()[0]){ 
          break; 
        }
      }
      else if (g->getType() == gr->getType() && g->getType() == rev::CX){
        if (g->getControls()[0] == gr->getControls()[0] && g->getTargets()[0] == gr->getTargets()[0]){
            //std::cout<<"Matched1"<<std::endl;
            int pos = std::distance(gates.begin(), gr.base()) - 1;
            gates.erase(g);
            //std::cout<<"Erase 2 gates"<<std::endl;
            gates.erase(--(gr.base()));
            g = gates.begin() + pos;
            if (g == gates.end()) --g;
            break;
        } 
        else if (g->getControls()[0] == gr->getTargets()[0] || g->getTargets()[0] == gr->getControls()[0]) break;
      }
      else if (g->getTargets()[0] == gr->getTargets()[0] && g->getControls().size() == 0 && gr->getControls().size() == 0
                                     && g->getTargets().size() == 1 && gr->getTargets().size() == 1){
            //std::cout<<"Matched2"<<std::endl;
            //g->display();
            //gr->display();
            gr->mergeGate(*g);
            //gr->display();
            gates.erase(g);
            //std::cout<<"Erase 1 gate"<<std::endl;
            int pos = std::distance(gates.begin(), gr.base()) - 1;            
            if (gr->isIdentity()) gates.erase(--(gr.base()));
            g = gates.begin() + pos;
            if (g == gates.end()) --g;
            break;
      }   
    } 
  } 
  ckt.setGates(gates);
}*/

double computeDepth(rev::Circuit & ckt){
  std::map<grph::node_t, double> depths = std::map<grph::node_t, double> ();
  std::vector<rev::Gate> gates = ckt.getGates();
  std::vector<rev::line_t> lines = ckt.getLines();
  //std::map<grph::node_t, bool > isFound = std::map<grph::node_t, bool > ();
  for(std::vector<rev::line_t>::iterator l=lines.begin(); l!=lines.end(); ++l){
    //isFound[*l] = false;
    depths[*l] = 0;
  }
  for(std::vector<rev::Gate>::iterator g=gates.begin(); g!=gates.end(); ++g){
    //g->display(std::cout);
    if (g->getType() == rev::CX){
      double cdepth = depths[g->getControls()[0]] > depths[g->getTargets()[0]] ? 
                                     depths[g->getControls()[0]] + 1 : depths[g->getTargets()[0]] + 1;
      depths[g->getControls()[0]] = cdepth;
      depths[g->getTargets()[0]] = cdepth;
    }
    else{
      depths[g->getTargets()[0]] = depths[g->getTargets()[0]] + 1;
    }      
  }
  double maxDepth = 0;
  //std::cout<<"\tDepth:";
  for(std::vector<rev::line_t>::iterator l=lines.begin(); l!=lines.end(); ++l){
    //isFound[*l] = false;
    //std::cout<<" "<<depths[*l];  
    if (maxDepth < depths[*l]) maxDepth = depths[*l];
  }
  //std::cout<<std::endl;  
  return maxDepth;
}

long computeCNOTGates(rev::Circuit & ckt){
  //std::map<grph::node_t, double> depths = std::map<grph::node_t, double> ();
  long CNOT_Gates = 0;
  std::vector<rev::Gate> gates = ckt.getGates();
  //std::vector<rev::line_t> lines = ckt.getLines();
  //std::map<grph::node_t, bool > isFound = std::map<grph::node_t, bool > ();
  /*for(std::vector<rev::line_t>::iterator l=lines.begin(); l!=lines.end(); ++l){
    //isFound[*l] = false;
    depths[*l] = 0;
  }*/
  for(auto & g : gates){
    //g->display(std::cout);
    if (g.getType() == rev::CX) CNOT_Gates += 1;      
  }  
  //std::cout<<std::endl;  
  return CNOT_Gates;
}

void generateRemoteCNOT(std::vector<rev::Gate>& gates, rev::Gate & g, std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& pg){
  rev::line_t control = phmap[g.getControls()[0]]; //mapped physical control qubit
  rev::line_t target = phmap[g.getTargets()[0]];  //mapped physical target qubit
  grph::distance_t d = pg.getDistance(control,target); 
  //std::cout<<" Distance: "<< d <<std::endl;
  //g.display(std::cout);
  while (d > 0 ){ //if distance between control and target qubit is greater than 0
     for (const auto & i : phmap){
       grph::distance_t d_new = pg.getDistance(i.second,target); 
       //find a mapped qubit adjacent to control qubit and closer to target qubit 
       if (control != i.second && target != i.second && (d - d_new) == 1) { 
         gates.push_back(rev::Gate(rev::CX, control, i.second));
         control = i.second;
         d = d_new;
         break;
       }       
     }
  }
  std::vector<rev::Gate> gates_seq = std::vector<rev::Gate> ();
  for ( auto g = gates.rbegin(); g !=gates.rend(); ++g) 
    gates_seq.push_back(rev::Gate(rev::CX, g->getControls()[0], g->getTargets()[0]));

  gates.push_back(rev::Gate(rev::CX, control, target)); 
  gates.insert(gates.end(), gates_seq.begin(), gates_seq.end());
  std::vector<rev::Gate> gates_seq2 = std::vector<rev::Gate> ();//(gates.begin() + 1, gates.end());
  for ( auto g = gates.begin() + 1; g + 1 !=gates.end(); ++g) 
    gates_seq2.push_back(rev::Gate(rev::CX, g->getControls()[0], g->getTargets()[0]));
  //std::cout<<"Hello2 at Generate Remote CNOT"<<std::endl;
  gates.insert(gates.end(), gates_seq2.begin(), gates_seq2.end()); 
  //std::cout<<"Hello3 at Generate Remote CNOT"<<gates.size()<<std::endl;
  
}
