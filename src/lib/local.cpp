#include "local.hpp"

//Hexagonal Approach 2
mincost_t minCostPathIndex(std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& lgmap,
                             std::map<rev::line_t, rev::line_t>& glmap, std::map<rev::line_t, rev::line_t>& plmap, 
                             std::map<grph::node_t, grph::node_t> & phmap, std::vector<rev::Gate>& gates,
                             std::vector<rev::line_t>& lines, std::pair<int, int> dim, grph::Graph & lg, 
                             unsigned pos, int base, Metric_Type metric){
  //double swap_cost = 0;
  //rev::line_t control;//, target;
  mincost_t mincost;
  mincost.swap_cost = 0;
  mincost.swaps = 0;
  //computing minimum cost insertion path index
  //for(std::vector<grph::node_t>::iterator p = path.begin(); p != path.end()-1; ++p){
  for(int i = 0; i < path.size() - 1; i++){
    rev::line_t cur_control = path[i];
    rev::line_t cur_target = path[i + 1];
    std::map<rev::line_t, rev::line_t> lgmap_tmp = std::map<rev::line_t, rev::line_t> (); //temporary local to global line map
    std::map<rev::line_t, rev::line_t> glmap_tmp = std::map<rev::line_t, rev::line_t> (); //temporary global to local line map
    for (const auto & l : lines) { //copy the current mapping to temporary location
      lgmap_tmp[l] = lgmap[l];	
      glmap_tmp[l] = glmap[l];	
    }
    double cur_swaps = 0;
    //forward move in path
    //for(std::vector<grph::node_t>::iterator pf = path.begin(); pf != p; ++pf){                  
    for(int j = 0; j < i; j++){  //physical qubit from path
      //perform swaping in temporary location
      cur_swaps += 1;
      //updating local to global map
      lgmap_tmp[glmap_tmp[plmap[path[j]]]] = plmap[path[j + 1]];
      lgmap_tmp[glmap_tmp[plmap[path[j + 1]]]] = plmap[path[j]];

      //updating global to local map
      rev::line_t l = glmap_tmp[plmap[path[j]]];
      glmap_tmp[plmap[path[j]]] = glmap_tmp[plmap[path[j + 1]]];
      glmap_tmp[plmap[path[j + 1]]] = l;    
    }
    //backward move in path
    //for(std::vector<grph::node_t>::reverse_iterator pb = path.rbegin(),  end(p + 2); pb != end; ++pb){
    for(int k = path.size() - 1; k > i + 1; k--){
      //perform swaping in temporary location
      cur_swaps += 1;
      //updating local to global map
      lgmap_tmp[glmap_tmp[plmap[path[k]]]] = plmap[path[k - 1]];
      lgmap_tmp[glmap_tmp[plmap[path[k - 1]]]] = plmap[path[k]];

      //updating global to local map
      rev::line_t l = glmap_tmp[plmap[path[k]]];
      glmap_tmp[plmap[path[k]]] = glmap_tmp[plmap[path[k - 1]]];
      glmap_tmp[plmap[path[k - 1]]] = l;    
    }
    
    //  grph::Graph lg = grph::Graph(lines, gates, pos + 1, grph::NEW, window, base);
    double cur_swap_cost;
    if (lg.isEmpty() == false) {
      if (metric == Metric_Type::COUPLING_COST)
        cur_swap_cost = mappingCost(phmap, lgmap_tmp, dim, lg);
      else
        cur_swap_cost = remoteCNOTCost(phmap, lgmap_tmp, dim, lg);
    }
    else 
      cur_swap_cost = cur_swaps;
    //std::cout<<"current control:"<<cur_control<<" current target:"<<cur_target<<" SWAP Cost: "<<cur_swap_cost<<std::endl;
    if (i == 0 || mincost.swap_cost > cur_swap_cost || (mincost.swap_cost == cur_swap_cost && mincost.swaps > cur_swaps)){
      mincost.swap_cost = cur_swap_cost;
      mincost.swaps = cur_swaps;
      mincost.control = cur_control;
      mincost.target = cur_target;
    }
  }
  //std::cout<<"current control:"<<mincost.control<<std::endl;
  //displayPath(path);
  return mincost;
}

//Hexagonal Approach 3
mincost_t minCostPathIndexv2(std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& lgmap,
                             std::map<rev::line_t, rev::line_t>& glmap, std::map<rev::line_t, rev::line_t>& plmap, 
                             std::map<grph::node_t, grph::node_t> & phmap, std::vector<rev::Gate>& gates,
                             std::vector<rev::line_t>& lines, std::pair<int, int> dim, grph::Graph & lg, 
                             unsigned pos, int base, Metric_Type metric){
  //double swap_cost = 0;
  //rev::line_t control;//, target;
  mincost_t mincost;
  mincost.swap_cost = 0;
  mincost.swaps = 0;
  //computing minimum cost insertion path index
  //for(std::vector<grph::node_t>::iterator p = path.begin(); p != path.end()-1; ++p){
  //lg.displayGraph();
  for(int i = 0; i < path.size() - 1; i++){
    rev::line_t cur_control = path[i];    
    std::map<rev::line_t, rev::line_t> lgmap_tmp = std::map<rev::line_t, rev::line_t> (); //temporary local to global line map
    std::map<rev::line_t, rev::line_t> glmap_tmp = std::map<rev::line_t, rev::line_t> (); //temporary global to local line map
    for (const auto & l : lines) { //copy the current mapping to temporary location
      lgmap_tmp[l] = lgmap[l];	
      glmap_tmp[l] = glmap[l];	
    }
    double cur_swaps = 0;
    //forward move in path
    //for(std::vector<grph::node_t>::iterator pf = path.begin(); pf != p; ++pf){                  
    for(int j = 0; j < i; j++){  //physical qubit from path
      //perform swaping in temporary location
      cur_swaps += 1;
      //updating local to global map
      lgmap_tmp[glmap_tmp[plmap[path[j]]]] = plmap[path[j + 1]];
      lgmap_tmp[glmap_tmp[plmap[path[j + 1]]]] = plmap[path[j]];

      //updating global to local map
      rev::line_t l = glmap_tmp[plmap[path[j]]];
      glmap_tmp[plmap[path[j]]] = glmap_tmp[plmap[path[j + 1]]];
      glmap_tmp[plmap[path[j + 1]]] = l;    
    }        
    //backward move in path
    //for(std::vector<grph::node_t>::reverse_iterator pb = path.rbegin(),  end(p + 2); pb != end; ++pb){
    for(int k = path.size() - 1; k > i ; k--){
      rev::line_t cur_target = path[k];

      if (k < path.size() - 1){
        //perform swaping in temporary location
        cur_swaps += 1;
        //updating local to global map
        lgmap_tmp[glmap_tmp[plmap[path[k]]]] = plmap[path[k + 1]];
        lgmap_tmp[glmap_tmp[plmap[path[k + 1]]]] = plmap[path[k]];

        //updating global to local map
        rev::line_t l = glmap_tmp[plmap[path[k]]];
        glmap_tmp[plmap[path[k]]] = glmap_tmp[plmap[path[k + 1]]];
        glmap_tmp[plmap[path[k + 1]]] = l;    
      }
      //  grph::Graph lg = grph::Graph(lines, gates, pos + 1, grph::NEW, window, base);
      double cur_swap_cost;                        
      if (lg.isEmpty() == false) {
        if (metric == Metric_Type::COUPLING_COST)
          cur_swap_cost = mappingCost(phmap, lgmap_tmp, dim, lg); 
        else
          cur_swap_cost = remoteCNOTCost(phmap, lgmap_tmp, dim, lg);
      }                                       //+ grph::getDistance(cur_control, cur_target, dim) * 4;
      else 
        cur_swap_cost = cur_swaps;// + grph::getDistance(cur_control, cur_target, dim) * 4;
      //std::cout<<"i:"<<i<<" k:"<<k<<std::endl;
      //std::cout<<"Control:"<<cur_control<<" target:"<<cur_target<<"\nswaps:"<<cur_swaps<<" C_SW_C:"<<cur_swap_cost<<std::endl;
      //if (lg.isEmpty() == false) 
        //std::cout<<"CC:"<<mappingCost(phmap, lgmap_tmp, dim, lg)<<" RC:"<<grph::getDistance(cur_control, cur_target, dim) * 4<<std::endl;

      //std::cout<<"current control:"<<cur_control<<" current target:"<<cur_target<<" SWAP Cost: "<<cur_swap_cost<<std::endl;
      if ((i == 0 && k == path.size() - 1) || mincost.swap_cost > cur_swap_cost 
                           || (mincost.swap_cost == cur_swap_cost && mincost.swaps > cur_swaps)){
        mincost.swap_cost = cur_swap_cost;
        mincost.swaps = cur_swaps;
        mincost.control = cur_control;
        mincost.target = cur_target;        
      }      
    }            
  } 
  //std::cout<<"Final:"<<std::endl;
  //std::cout<<"Control:"<<mincost.control<<" target:"<<mincost.target<<"\nswaps:"<<mincost.swaps<<" cc:"<<mincost.swap_cost<<std::endl;
  //std::cout<<"current control:"<<mincost.control<<std::endl;
  //displayPath(path);
  return mincost;
}


//Hexagonal Mapping Approach 1 : Metric = Coupling Cost, Global Ordering, Remote CNOT gate
void localorder(rev::Circuit & ckt, std::map<grph::node_t, grph::node_t> & phmap, std::pair<int, int> dim){
  std::vector<rev::Gate> gates = ckt.getGates();
  std::vector<rev::line_t> lines = ckt.getLines();
  //unsigned size = gates.size();
  
  std::map<rev::line_t, rev::line_t> plmap = std::map<rev::line_t, rev::line_t> (); //physical to logical line map
  std::vector<rev::line_t> new_lines = std::vector<rev::line_t> ();   
  
  for (const auto & l : lines){
    plmap[phmap[l]] = l;
    new_lines.push_back(phmap[l]);
  }
  /*std::cout<<"All Maps";
  for(std::map<rev::line_t, rev::line_t>::iterator l = plmap.begin(); l != plmap.end(); ++l)
        std::cout<<l->first<<":"<<l->second<<" ";
  std::cout<<std::endl;*/
  //unsigned gcount = 0;
  //unsigned tot_gate = gates.size();
  std::vector<rev::Gate> new_gates;
  for(auto g = gates.begin(); g!=gates.end(); ++g){ 
    if(1 == g->getControls().size()){
      //gcount++;  
      //std::cout<<"Hello1:"<<gcount<<std::endl;  
      if (grph::getDistance(phmap[g->getControls()[0]], phmap[g->getTargets()[0]], dim)==0){ // adjacent gate
        //std::cout<<"Hellox"<<std::endl;  
        g->updateLines(phmap[g->getControls()[0]],phmap[g->getTargets()[0]]);
        new_gates.push_back(*g);
        //g->display(std::cout);
      }
      else { //qubits are not adjacent
        //std::cout<<"Helloy"<<std::endl;        
        std::vector<std::vector<grph::node_t> > paths = std::vector<std::vector<grph::node_t> > ();
        paths.push_back(std::vector<grph::node_t> ());

        getAllPaths(phmap[g->getControls()[0]],phmap[g->getTargets()[0]], paths, dim, 0);
        //std::cout<<"Helloz"<<std::endl;        
        unsigned pos = std::distance(gates.begin(),g);
        
        //std::cout<<"Hello2"<<std::endl;
        int i = 0; //path index   
        for ( ; i < paths.size(); i++){
          //std::cout<<"Hello1"<<std::endl;
          if (!isInvalidPath(paths[i], plmap)) {
            //std::cout<<"Hello1"<<std::endl;
            break;
          }
          //std::cout<<"Hello2"<<paths.size()<<std::endl;
          /*std::cout<<"Path: "<<i<<":";
    	  for (auto node : paths[i]){
            /std::cout<<"Next Node: "<<node<<" ";
          }
          std::cout<<std::endl;*/          
        }
        //else { mincost.control = paths[k][0];
        //std::cout<<"selected widnow: "<<window<<" swap_cost: "<<mincost.swap_cost<<" swaps: "<<mincost.swaps<<std::endl;
        //std::cout<<"Hello3: "<<k<<" control: "<<mincost.control<<std::endl;
        //g->display();
        //displayPath(paths[k]);
        //std::cout<<"Hello 3: i="<<i<<std::endl;
        std::vector<rev::Gate> rcnot_gates = std::vector<rev::Gate> ();
        //Actual swap gate insertion in forward path
        //insertRCNOT(rcnot_gates, *g, paths[i], phmap[g->getControls()[0]], phmap[g->getTargets()[0]]);
        insertRCNOT(rcnot_gates, paths[i], phmap[g->getControls()[0]], phmap[g->getTargets()[0]]);
        //std::cout<<"Hello 4"<<std::endl;        
        //g->updateLines(mincost.control, mincost.target);
        //g->display(std::cout);
        //std::cout<<"Hello 5"<<std::endl;
        //gates.insert(g, rcnot_gates.begin(),rcnot_gates.end());        
        new_gates.insert(new_gates.end(), rcnot_gates.begin(),rcnot_gates.end());
        //std::cout<<"Hello 6"<<std::endl;
        //g = gates.begin() + pos + swap_gates.size();
        //gates.erase(g);
        //g = gates.begin() + pos + rcnot_gates.size();
        //std::cout<<"Hello 7"<<std::endl;
        //if (g == gates.end() ) break;
        //std::cout<<"Gates: "<<gates.size()<<" swaps: "<<swap_gates.size()<<" actual gates: "<<gates.size()-swap_gates.size()<<std::endl;
        //ckt.setGates(gates);
      }
    }
    else{
      g->updateLines(phmap[g->getTargets()[0]]);
      new_gates.push_back(*g);
      //g->display();
    //g->display(std::cout);
    }
    //std::cout<<"Hello6"<<std::endl;
  }
  ckt.setGates(new_gates);
  ckt.setLines(new_lines);
}

//Hexagonal Mapping Approach 2 : Metric =Coupling Cost, Global Ordering, Local Ordering, SWAP gate
void localorder(rev::Circuit & ckt, std::map<grph::node_t, grph::node_t> & phmap, std::pair<int, int> dim, 
                                                                    grph::Graph & lg, int base, Metric_Type metric){
  grph::Graph gr = lg;
  std::vector<rev::Gate> gates = ckt.getGates();
  std::vector<rev::line_t> lines = ckt.getLines();
  //unsigned size = gates.size();
  std::map<rev::line_t, rev::line_t> lgmap = std::map<rev::line_t, rev::line_t> (); //local to global line map
  std::map<rev::line_t, rev::line_t> glmap = std::map<rev::line_t, rev::line_t> (); //global to local line map
  std::map<rev::line_t, rev::line_t> plmap = std::map<rev::line_t, rev::line_t> (); //physical to logical line map
  std::vector<rev::line_t> new_lines = std::vector<rev::line_t> (); 
  
  for (const auto & l : lines){
    lgmap[l] = l;	
    glmap[l] = l;
    plmap[phmap[l]] = l;
    new_lines.push_back(phmap[l]);
  }
  /*std::cout<<"All Maps";
  for(std::map<rev::line_t, rev::line_t>::iterator l = plmap.begin(); l != plmap.end(); ++l)
        std::cout<<l->first<<":"<<l->second<<" ";
  std::cout<<std::endl;*/
  //std::cout<<"local ordering Approach 2, window"<<window<<std::endl;
  //unsigned gcount = 0;
  //unsigned tot_gate = gates.size();
  std::vector<rev::Gate> new_gates;
  for(auto g = gates.begin(); g!=gates.end(); ++g){ 
    if(1 == g->getControls().size()){
      //gcount++;  
      //std::cout<<"local ordering 21"<<std::endl;  
      if (grph::getDistance(phmap[lgmap[g->getControls()[0]]], phmap[lgmap[g->getTargets()[0]]], dim)==0){ // adjacent gate
        g->updateLines(phmap[lgmap[g->getControls()[0]]],phmap[lgmap[g->getTargets()[0]]]);
        new_gates.push_back(*g);
        //g->display(std::cout);
      }
      else { //qubits are not adjacent
              
        std::vector<std::vector<grph::node_t> > paths = std::vector<std::vector<grph::node_t> > ();
        paths.push_back(std::vector<grph::node_t> ());

        getAllPaths(phmap[lgmap[g->getControls()[0]]],phmap[lgmap[g->getTargets()[0]]], paths, dim, 0);

        //std::cout<<"Hello1"<<std::endl;
        //unsigned window;
        //if ( base == 2) window = 125;
        //else if ( base == 3) window = 80;
        //else if ( base == 4) window = 60;
        unsigned pos = std::distance(gates.begin(),g);        
        if (g + 1 != gates.end()){
          //g->display(std::cout);
          //std::cout<<"Hello2"<<window<<std::endl;
          if (metric == Metric_Type::COUPLING_COST)
            gr = grph::Graph(lines, gates, pos + 1, grph::NEW, window, base);
          else
            //gr.setWeight(g->getControls()[0], g->getTargets()[0], 
            //             lg.getWeight(g->getControls()[0], g->getTargets()[0]) - 1.0); //entire circuit
            
            gr = grph::Graph(lines, gates, pos + 1, window); // constant size window
            //gr = grph::Graph(lines, gates, pos + 1, 6, grph::window_type::variable); //variable size window, maxdeg for hexagonal architecture = 6
          //gr.displayGraph(); 
        }
        //std::cout<<"Hello2"<<std::endl;
        int k = 0; //path index
        mincost_t mincost;   
        bool empty = true;
        for (int i=0; i < paths.size(); i++){
          mincost_t tempcost;
          if (isInvalidPath(paths[i], plmap)) {
            //std::cout<<"Hello1"<<std::endl;
            continue;
          }
          //std::cout<<"Hello2A"<<paths.size()<<std::endl;
          /*std::cout<<"Path: "<<i<<":";
    	  for (auto node : paths[i]){
            /std::cout<<"Next Node: "<<node<<" ";
          }
          std::cout<<std::endl;*/
          tempcost = minCostPathIndex(paths[i], lgmap, glmap, plmap, phmap, gates, lines, dim, gr, pos, base, metric);
          //std::cout<<"current control:"<<tempcost.control<<" current target:"<<tempcost.target<<" SWAP Cost: "<<tempcost.swap_cost<<std::endl;
          //std::cout<<"Hello2B"<<paths.size()<<std::endl;
          if (empty || mincost.swap_cost > tempcost.swap_cost || (mincost.swap_cost == tempcost.swap_cost && mincost.swaps > tempcost.swaps)) {
            mincost.swap_cost = tempcost.swap_cost;
            mincost.swaps = tempcost.swaps;
            mincost.control = tempcost.control;
            mincost.target = tempcost.target;
            k = i;
            empty = false;
          }
          if (paths.size() * paths[0].size() > 20) break; //Minimize search when number of path and nodes in path increases
        }
        //else { mincost.control = paths[k][0];
        //std::cout<<"selected widnow: "<<window<<" swap_cost: "<<mincost.swap_cost<<" swaps: "<<mincost.swaps<<std::endl;
        //std::cout<<"Hello3: "<<k<<" control: "<<mincost.control<<std::endl;
        //g->display();
        //displayPath(paths[k]);
        //std::cout<<"Hello 3"<<std::endl;
        std::vector<rev::Gate> swap_gates = std::vector<rev::Gate> ();
        //Actual swap gate insertion in forward path
        insertSWAP(swap_gates, paths[k], lgmap, glmap, plmap, mincost.control, mincost.target);     
        //std::cout<<"Hello 4"<<std::endl;        
        g->updateLines(mincost.control, mincost.target);
        //g->display(std::cout);
        //std::cout<<"Hello 5"<<std::endl;
        //gates.insert(g,swap_gates.begin(),swap_gates.end());        
        new_gates.insert(new_gates.end(), swap_gates.begin(),swap_gates.end()); 
        new_gates.push_back(*g);
        //std::cout<<"Hello 6"<<std::endl;
        //g = gates.begin() + pos + swap_gates.size();
        //gates.erase(g);
        //g = gates.begin() + pos + swap_gates.size();
        //std::cout<<"Hello 7"<<std::endl;
        //if (g == gates.end() ) break;
        //std::cout<<"Gates: "<<gates.size()<<" swaps: "<<swap_gates.size()<<" actual gates: "<<gates.size()-swap_gates.size()<<std::endl;
        //ckt.setGates(gates);
      }
    }
    else{
      g->updateLines(phmap[lgmap[g->getTargets()[0]]]);
      new_gates.push_back(*g);
      //g->display();
    //g->display(std::cout);
    }
    //std::cout<<"Hello6"<<std::endl;
  }
  ckt.setGates(new_gates);
  ckt.setLines(new_lines);
}

//Hexagonal Mapping Approach 3 : Metric =Coupling Cost, Global Ordering, Local Ordering, SWAP + RCNOT gate
void localorderv2(rev::Circuit & ckt, std::map<grph::node_t, grph::node_t> & phmap, std::pair<int, int> dim, 
                                                                   grph::Graph & lg, int base, Metric_Type metric){
  grph::Graph gr = lg;
  std::vector<rev::Gate> gates = ckt.getGates();
  std::vector<rev::line_t> lines = ckt.getLines();
  //unsigned size = gates.size();
  std::map<rev::line_t, rev::line_t> lgmap = std::map<rev::line_t, rev::line_t> (); //local to global line map
  std::map<rev::line_t, rev::line_t> glmap = std::map<rev::line_t, rev::line_t> (); //global to local line map
  std::map<rev::line_t, rev::line_t> plmap = std::map<rev::line_t, rev::line_t> (); //physical to logical line map
  std::vector<rev::line_t> new_lines = std::vector<rev::line_t> (); 
  
  for (const auto & l : lines){
    lgmap[l] = l;	
    glmap[l] = l;
    plmap[phmap[l]] = l;
    new_lines.push_back(phmap[l]);
  }
  /*std::cout<<"All Maps";
  for(std::map<rev::line_t, rev::line_t>::iterator l = plmap.begin(); l != plmap.end(); ++l)
        std::cout<<l->first<<":"<<l->second<<" ";
  std::cout<<std::endl;*/
  //std::cout<<"local ordering Approach 3, window"<<window<<std::endl;
  //unsigned gcount = 0;
  //unsigned tot_gate = gates.size();
  std::vector<rev::Gate> new_gates;
  for(auto g = gates.begin(); g!=gates.end(); ++g){ 
    if(1 == g->getControls().size()){
      //gcount++;          
      if (grph::getDistance(phmap[lgmap[g->getControls()[0]]], phmap[lgmap[g->getTargets()[0]]], dim)==0){ // adjacent gate
        g->updateLines(phmap[lgmap[g->getControls()[0]]],phmap[lgmap[g->getTargets()[0]]]);
        new_gates.push_back(*g);
        //g->display(std::cout);
      }
      else { //qubits are not adjacent
              
        std::vector<std::vector<grph::node_t> > paths = std::vector<std::vector<grph::node_t> > ();
        paths.push_back(std::vector<grph::node_t> ());

        getAllPaths(phmap[lgmap[g->getControls()[0]]],phmap[lgmap[g->getTargets()[0]]], paths, dim, 0);

        //std::cout<<"Hello1"<<std::endl;
        //unsigned window;
        //if ( base == 2) window = 125;
        //else if ( base == 3) window = 80;
        //else if ( base == 4) window = 60;
        unsigned pos = std::distance(gates.begin(),g);
        if (g + 1 != gates.end()){
          //g->display(std::cout);
          //std::cout<<"Hello2"<<window<<std::endl;
          if (metric == Metric_Type::COUPLING_COST)
            gr = grph::Graph(lines, gates, pos + 1, grph::NEW, window, base);
          else            
            //gr.setWeight(g->getControls()[0], g->getTargets()[0], 
            //             lg.getWeight(g->getControls()[0], g->getTargets()[0]) - 1.0); //entire circuit
            
            gr = grph::Graph(lines, gates, pos + 1, window); // constant size window
            //gr = grph::Graph(lines, gates, pos + 1, 6, grph::window_type::variable); //variable size window, maxdeg for hexagonal architecture = 6
          //gr.displayGraph(); 
        }
        //std::cout<<"Hello2"<<std::endl;
        int k = 0; //path index
        mincost_t mincost;   
        bool empty = true;
        for (int i=0; i < paths.size(); i++){
          mincost_t tempcost;
          if (isInvalidPath(paths[i], plmap)) {
            //std::cout<<"Hello1"<<std::endl;
            continue;
          }
          //std::cout<<"Hello2"<<paths.size()<<std::endl;
          /*std::cout<<"Path: "<<i<<":";
    	  for (auto node : paths[i]){
            /std::cout<<"Next Node: "<<node<<" ";
          }
          std::cout<<std::endl;*/
          //std::cout<<"Hello2A"<<std::endl;
          tempcost = minCostPathIndexv2(paths[i], lgmap, glmap, plmap, phmap, gates, lines, dim, gr, pos, base, metric);
          //std::cout<<"Hello3"<<std::endl;
          //std::cout<<"current control:"<<tempcost.control<<" current target:"<<tempcost.target<<" SWAP Cost: "<<tempcost.swap_cost<<std::endl;
          if (empty || mincost.swap_cost > tempcost.swap_cost || (mincost.swap_cost == tempcost.swap_cost && mincost.swaps > tempcost.swaps)) {
            mincost.swap_cost = tempcost.swap_cost;
            mincost.swaps = tempcost.swaps;
            mincost.control = tempcost.control;
            mincost.target = tempcost.target;
            k = i;
            empty = false;
          }
          if (paths.size() * paths[0].size() > 20) break;
        }
        //else { mincost.control = paths[k][0];
        //std::cout<<"selected widnow: "<<window<<" swap_cost: "<<mincost.swap_cost<<" swaps: "<<mincost.swaps<<std::endl;
        //std::cout<<"Hello3: "<<k<<" control: "<<mincost.control<<std::endl;
        //g->display();
        //displayPath(paths[k]);
        //std::cout<<"Hello 3"<<std::endl;
        std::vector<rev::Gate> swap_gates = std::vector<rev::Gate> ();
        //Actual swap gate insertion in forward path
        //std::cout<<"Hello 4B"<<std::endl;
        insertSWAP(swap_gates, paths[k], lgmap, glmap, plmap, mincost.control, mincost.target);
        new_gates.insert(new_gates.end(), swap_gates.begin(),swap_gates.end());         
        //std::cout<<"Hello 4"<<swap_gates.size()<<std::endl;        
        if (grph::getDistance(mincost.control, mincost.target, dim) > 0){
          std::vector<rev::Gate> rcnot_gates = std::vector<rev::Gate> ();
          //Actual swap gate insertion in forward path
          //std::cout<<"Hello 5B"<<std::endl;
          //insertRCNOT(rcnot_gates, *g, paths[k], mincost.control, mincost.target);
          insertRCNOT(rcnot_gates, paths[k], mincost.control, mincost.target);
          new_gates.insert(new_gates.end(), rcnot_gates.begin(), rcnot_gates.end());
          //std::cout<<"Hello 5"<<rcnot_gates.size()<<std::endl; 
          //merging gates
          //swap_gates.insert(swap_gates.end(), rcnot_gates.begin(), rcnot_gates.end());
        }
        else{
          g->updateLines(mincost.control, mincost.target);
          new_gates.push_back(*g);
        }
        //g->display(std::cout);
        //std::cout<<"Hello 5"<<std::endl;
        //gates.insert(g,swap_gates.begin(),swap_gates.end());        
        //std::cout<<"Hello 6"<<std::endl;
        //g = gates.begin() + pos + swap_gates.size();
        //gates.erase(g);
        //g = gates.begin() + pos + swap_gates.size();
        //std::cout<<"Hello 7"<<std::endl;
        //if (g == gates.end() ) break;
        //std::cout<<"Gates: "<<gates.size()<<" swaps: "<<swap_gates.size()<<" actual gates: "<<gates.size()-swap_gates.size()<<std::endl;
        //ckt.setGates(gates);
      }
    }
    else{
      g->updateLines(phmap[lgmap[g->getTargets()[0]]]);
      new_gates.push_back(*g);
      //g->display();
    //g->display(std::cout);
    }
    //std::cout<<"Hello6"<<std::endl;
  }
  ckt.setGates(new_gates);
  ckt.setLines(new_lines);
}

//void insertRCNOT(std::vector<rev::Gate>& rcnot_gates, rev::Gate & g, std::vector<grph::node_t>& path, 
//                             rev::line_t control, rev::line_t target){
void insertRCNOT(std::vector<rev::Gate>& rcnot_gates, std::vector<grph::node_t>& path, 
                             rev::line_t control, rev::line_t target){
  bool flag = false;
  //std::cout<<"Hello 5x1"<<std::endl;
  for (int i = 0; i < path.size() - 1; i++){
    if (path[i] == control) flag = true;
    if (flag == true){
      rcnot_gates.push_back(rev::Gate(rev::CX, path[i], path[i + 1]));
    }
    if (path[i + 1] == target) break;
  }
  //std::cout<<"Hello 5x2"<<std::endl;      
  for (int i = rcnot_gates.size() - 2; i >= 0; i --){ 
    rcnot_gates.push_back(rev::Gate(rev::CX, rcnot_gates[i].getControls()[0], rcnot_gates[i].getTargets()[0]));
  }
  //std::cout<<"Hello 5x3"<<std::endl;
  for (int i = rcnot_gates.size() - 2; i >= 1; i --){ 
    rcnot_gates.push_back(rev::Gate(rev::CX, rcnot_gates[i].getControls()[0], rcnot_gates[i].getTargets()[0]));
  }
  //std::cout<<"Hello 5x4"<<std::endl;
  
  //std::cout<<rcnot_gates[1].getControls()[0]<<std::endl;
  //std::cout<<"Hello 5x5"<<std::endl;
  //std::cout<<rcnot_gates[1].getTargets()[0]<<std::endl;
  //std::cout<<"Hello 5x6"<<std::endl;
  //g.updateLines(rcnot_gates[1].getControls()[0], rcnot_gates[1].getTargets()[0]);
  //std::cout<<"Hello 5x7"<<std::endl;  
}

void insertSWAP(std::vector<rev::Gate>& swap_gates, std::vector<grph::node_t>& path, std::map<rev::line_t, rev::line_t>& lgmap,
                             std::map<rev::line_t, rev::line_t>& glmap, std::map<rev::line_t, rev::line_t>& plmap, 
                             rev::line_t control, rev::line_t target){
  for(int i = 0; path[i] != control; i++){                  
    //updating local to global map
    lgmap[glmap[plmap[path[i]]]] = plmap[path[i + 1]];
    lgmap[glmap[plmap[path[i + 1]]]] = plmap[path[i]];

    //updating global to local map
    rev::line_t l = glmap[plmap[path[i]]];
    glmap[plmap[path[i]]] = glmap[plmap[path[i + 1]]];
    glmap[plmap[path[i + 1]]] = l;              

           
    swap_gates.push_back(rev::Gate(rev::CX, path[i], path[i + 1]));           
    swap_gates.push_back(rev::Gate(rev::CX, path[i + 1], path[i]));
    swap_gates.push_back(rev::Gate(rev::CX, path[i], path[i + 1]));

  }
  //std::cout<<"Hello4"<<std::endl;
  //Actual swap gate insertion in backward path
  for(int j = path.size() - 1; path[j] != target; j--){
  //for(std::vector<grph::node_t>::reverse_iterator pb = paths[k].rbegin(); *pb != mincost.target; ++pb){
    //updating local to global map
    lgmap[glmap[plmap[path[j]]]] = plmap[path[j - 1]];
    lgmap[glmap[plmap[path[j - 1]]]] = plmap[path[j]];

    //updating global to local map
    rev::line_t l = glmap[plmap[path[j]]];
    glmap[plmap[path[j]]] = glmap[plmap[path[j - 1]]];
    glmap[plmap[path[j - 1]]] = l;    

    swap_gates.push_back(rev::Gate(rev::CX, path[j], path[j - 1]));
    swap_gates.push_back(rev::Gate(rev::CX, path[j - 1], path[j]));
    swap_gates.push_back(rev::Gate(rev::CX, path[j], path[j - 1]));
  }
}
