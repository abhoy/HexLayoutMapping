#include <cmath>
#include <cstdlib>
#include "global.hpp"


//Generate the initial population, each solution being of size equal to the size of any element from phmaps. 
void gen_initial_pop(std::vector<std::map<grph::node_t, grph::node_t> > & phmaps, std::map<grph::node_t, grph::node_t> & phmap, grph::Graph& lg){ //grph::Graph& pg, 
  //Random placement of vertices in 1D grid
  for(std::vector<std::map<grph::node_t, grph::node_t> >::iterator i=phmaps.begin(); i!=phmaps.end(); ++i){
    physicalMap(*i, phmap, lg);
  }
}

//Mapping logical qubits to physial qubits
void physicalMap(std::map<grph::node_t, grph::node_t> & phmap_new, std::map<grph::node_t, grph::node_t> & phmap, grph::Graph& lg){ //grph::Graph& pg, 
  int r;
  int x = lg.getNodes().size();
  //std::vector<grph::node_t> nodes = pg.getNodes();
  //for(std::vector<grph::node_t>::iterator n = nodes.begin(); n != nodes.end() && phmap_new.size() < x; ++n){
  for(std::map<grph::node_t, grph::node_t>::iterator n = phmap.begin(); n != phmap.end() ; ++n){
    do{
      r = rand() % x;  
      //std::cout<<r<<std::endl;
    }while(phmap_new.find(r) != phmap_new.end());
    phmap_new[r] = n->second;
  }
}

//Calculate the fitness values of the solutions in "phmaps". 
//IBM QX mapping
void calculate_fitness(std::vector<std::map<grph::node_t, grph::node_t> > & phmaps, std::vector<double>& fitness, 
                                                          grph::Graph& pg, grph::Graph& lg, int base, unsigned depth){
  fitness.clear();
  for(std::vector<std::map<grph::node_t, grph::node_t> >::iterator i=phmaps.begin(); i!=phmaps.end(); ++i){
    //std::cout<<"mappingCost(*i, pg, lg) "<<mappingCost(*i, pg, lg)<<std::endl;

    /*std::cout<<"Physical Mapping:"<<std::endl;
    for(std::vector<grph::node_t>::iterator j=lg.getNodes().begin(); j!=lg.getNodes().end(); ++j){
      std::cout<<setw(2)<<(*j)<<" ";
      //lmap[*j] = phmap_min[*j];
      //lmap2[*j] = *j;
    }
    std::cout<<std::endl;
    for(std::vector<grph::node_t>::iterator j=lg.getNodes().begin(); j!=lg.getNodes().end(); ++j){
      std::cout<<setw(2)<<(*i)[*j]<<" ";
    }
    std::cout<<std::endl<<"cost: "<<swapCost(*i, pg, lg, base, depth)<<std::endl;*/
    //std::cout<<"Hello12-1"<<std::endl;
    fitness.push_back(swapCost(*i, pg, lg, base, depth));//mappingCost(*i, pg, lg));
  }
}

//Calculate the fitness values of the solutions in "phmaps". 
//HEXAGONAL ARCHITECTURE mapping
void calculate_fitness(std::vector<std::map<grph::node_t, grph::node_t> > & phmaps, std::vector<double>& fitness, 
                                       std::pair<int, int> dim, grph::Graph& lg, Metric_Type metric){
  fitness.clear();
  for (auto & m : phmaps){
    if (metric == Metric_Type::COUPLING_COST)
      fitness.push_back(mappingCost(m, dim, lg));  // coupling cost
    else
      fitness.push_back(remoteCNOTCost(m, dim, lg)); //remote CNOT COST
  }
}

//Calculate percentage fitness of the solutions 
void  calculate_percent_fitness(std::vector<double>& fitness, std::vector<double>& percent_fitness){
  double total=0;
  for(std::vector<double>::iterator i=fitness.begin(); i!=fitness.end(); ++i){
    //std::cout<<"*i "<<(*i)<<" total "<<total<<std::endl;
    total += *i;
  }
  percent_fitness.clear();
  for(std::vector<double>::iterator i=fitness.begin(); i!=fitness.end(); ++i){
    if (total == 0.) percent_fitness.push_back((*i));
    else percent_fitness.push_back((*i)/total);
  }
}

//Sort the solutions in "phmaps" in ascending order of fitness values
void sort_gener_cur(std::vector<std::map<grph::node_t, grph::node_t> > & phmaps, std::vector<double>& fitness){
  for(std::vector<double>::iterator i=fitness.begin(); i!=fitness.end() - 1; ++i){
    for(std::vector<double>::iterator j= i + 1; j!=fitness.end(); ++j){
      if(*i > *j){
        double t = *i;
        *i = *j;
        *j = t;
        for(std::map<grph::node_t, grph::node_t>::iterator m=phmaps[std::distance(fitness.begin(),i)].begin(),
                                                           n=phmaps[std::distance(fitness.begin(),j)].begin(); 
                                                           m!=phmaps[std::distance(fitness.begin(),i)].end(); ++m, ++n){
          //grph::node_t tmp1 = m->first;
          grph::node_t tmp = m->second;
          //m->first = n->first;
          m->second = n->second;
          //n->first = tmp1;
          n->second = tmp;
        }
      }
    }
  }
}

//Select a solution from "phmaps" using roulette wheel method 
int  selection(std::vector<double>& percent_fitness){
  int  i;
  long double  r, cumul_sum;
      
  r = rand()/(double)RAND_MAX;
  cumul_sum = 0.0;

  for(std::vector<double>::iterator i=percent_fitness.begin(); i!=percent_fitness.end(); ++i){
    //std::cout<<"*i "<<(*i)<<" cumul_sum "<<cumul_sum<<r<<std::endl;
    cumul_sum += *i; 
    if  (r < cumul_sum)  return std::distance(percent_fitness.begin(),i);
  }
  //std::cout<<"Hello152 r "<<r<<" cumul_sum "<<cumul_sum<<std::endl;
  return 0; //Default
}

//Perform crossover of two solutions of "phmaps" at indices "pos1" and "pos2", and 
//store the new solutions in "phmaps_next" at the end.
void crossover_xy (std::vector<std::map<grph::node_t, grph::node_t> > & phmaps, 
                    std::vector<std::map<grph::node_t, grph::node_t> > & phmaps_next, int pos1, int pos2, int x){
  
  int crosspoint;
  do {
    crosspoint = rand() % (x-1);
  }while(x == crosspoint + 2); //Crossover requires at least two items
  //std::cout<<"Hello1541: Crosspoint:"<<crosspoint<<": x:"<<x<<std::endl;
  std::map<grph::node_t, grph::node_t> phmap1, phmap2;
  /*std::cout<<"Pos1: "<<pos1<<std::endl;;
  for (auto i : phmaps[pos1]){
    std::cout<<"\t"<<i.first<<":"<<i.second<<std::endl;
  }
  std::cout<<std::endl;
  std::cout<<"Pos2: "<<pos2<<std::endl;
  for (auto i : phmaps[pos2]){
    std::cout<<"\t"<<i.first<<":"<<i.second<<std::endl;
  }
  std::cout<<std::endl;*/

  for(auto i=phmaps[pos1].begin(), j=phmaps[pos2].begin(); i!=phmaps[pos1].end() && 
      std::distance(phmaps[pos1].begin(),i)<=crosspoint; ++i, ++j){  // Copy all elements upto the crossover point
    phmap1[i->first] = i->second;
    phmap2[j->first] = j->second;
    //std::cout<<"Copy:"<<i->first<<":"<<i->second<<std::endl;
  }
  //std::cout<<"Hello1542"<<std::endl;
  //std::cout<<"Crosspoint:"<<crosspoint<<std::endl;
  int index1, index2;
  std::map<int, bool> snodes1, snodes2; 
  std::map<grph::node_t, bool> crossover1, crossover2;
  int items = x - crosspoint - 1; //number of items to be cross over
  for (auto i = 0; i < items; i++)  //initially no crossover 
    snodes1[i] = false, snodes2[i] = false;  
  //std::cout<<"Hello1543"<<std::endl;
  while (items > 0){ //perform crossover
    do{
      index1 = rand() % (x-crosspoint-1); //select a mapping that is not crossover
    }while(snodes1[index1]);      
    //std::cout<<"Hello1544"<<std::endl;
    do{
      index2 = rand() % (x-crosspoint-1); //select another mapping that is not crossover
      //std::cout<<index2<<" "<<x-crosspoint-1<<std::endl;
    }while(snodes1[index2] || index1 == index2);  
    snodes1[index1] = true, snodes1[index2] = true;
    auto i = phmaps[pos1].begin(), j = phmaps[pos1].begin(); 
    std::advance(i, crosspoint + 1 + index1), std::advance(j, crosspoint + 1 + index2);   
    phmap1[i->first] = j->second, phmap1[j->first] = i->second;
    crossover1[i->second] = true, crossover1[j->second] = true;
    //std::cout<<"Index1: "<<index1<<"- Index2: "<<index2<<std::endl;
    //std::cout<<"Crossovered: "<<i->first<<"-"<<i->second<<" and "<<j->first<<"-"<<j->second<<std::endl;
    //std::cout<<"Hello1545"<<std::endl;
    do{
      index1 = rand() % (x-crosspoint-1); //select a mapping that is not crossover
    }while(snodes2[index1]);  
    //std::cout<<"Hello1546"<<std::endl;
    do{
      index2 = rand() % (x-crosspoint-1); //select another mapping that is not crossover
      //std::cout<<index2<<" "<<x-crosspoint-1<<std::endl;
    }while(snodes2[index2] || index1 == index2);  
    snodes2[index1] = true, snodes2[index2] = true;
    i = phmaps[pos2].begin(), j = phmaps[pos2].begin(); 
    std::advance(i, crosspoint + 1 + index1), std::advance(j, crosspoint + 1 + index2);   
    phmap2[i->first] = j->second, phmap2[j->first] = i->second;
    crossover2[i->second] = true, crossover2[j->second] = true;
    //std::cout<<"Index1: "<<index1<<"- Index2: "<<index2<<std::endl;
    //std::cout<<"Crossovered: "<<i->first<<"-"<<i->second<<" and "<<j->first<<"-"<<j->second<<std::endl;
    //std::cout<<"Hello1547"<<std::endl;
    items -=2; //two items crossovered
    if (items == 1) { //can not perform crossover
      i = phmaps[pos1].begin(), j = phmaps[pos2].begin();
      for (std::advance(i, crosspoint + 1); i != phmaps[pos1].end(); i++){ //copy items that are not crossover
        if (crossover1.find(i->second) == crossover1.end()) {
          phmap1[i->first] = i->second; //copy the item
          break;
        }
      }
      for (std::advance(j, crosspoint + 1); j != phmaps[pos2].end(); j++){
        if (crossover2.find(j->second) == crossover2.end()) {
          phmap2[j->first] = j->second; 
          break;
        }
      }
      break; 
    }
  }
      
  /*std::cout<<"Pos1: "<<pos1<<std::endl;;
  for (auto i : phmap1){
    std::cout<<"\t"<<i.first<<":"<<i.second<<std::endl;
  }
  std::cout<<std::endl;
  std::cout<<"Pos2: "<<pos2<<std::endl;
  for (auto i : phmap2){
    std::cout<<"\t"<<i.first<<":"<<i.second<<std::endl;
  }
  std::cout<<std::endl;*/
  //exit(0);
  phmaps_next.push_back(phmap1);
  phmaps_next.push_back(phmap2);
}


//Generate next generation "phmaps_next" from the current generation "phmaps"
//Copy the top NBEST solutions as it is; rest are generated through crossover using probability CROSSPROB 
void  generation_next(std::vector<std::map<grph::node_t, grph::node_t> > & phmaps, 
                    std::vector<std::map<grph::node_t, grph::node_t> > & phmaps_next, 
                    std::vector<double>& percent_fitness, int x){
    
  phmaps_next.clear();
  // Copy best "NBEST" solutions
  //std::cout<<"Hello151"<<std::endl;
  for(std::vector<std::map<grph::node_t, grph::node_t> >::iterator i = phmaps.begin(); i != phmaps.end() &&
           std::distance(phmaps.begin(),i)<NBEST; ++i){
    std::map<grph::node_t, grph::node_t> phmap;
    for(std::map<grph::node_t, grph::node_t>::iterator j = (*i).begin(); j != (*i).end(); ++j) phmap[j->first] = j->second;
    phmaps_next.push_back(phmap);    
  }
  //std::cout<<"Hello152"<<std::endl;
  int num_crossover = (POPSIZE - NBEST) / 2;  /* Number of crossovers to perform POPSIZE & NBIST should be even */
 
  for(int k = 0; k < num_crossover; k++){
    int pos1 = selection(percent_fitness);
    int pos2 = selection(percent_fitness);
    //if (pos1 == 0 && pos2 == 0) break;
    //std::cout<<"Hello153 pos1 "<<pos1<<" pos2 "<<pos2<<std::endl;
    if  ((rand()/(float)RAND_MAX) < CROSSPROB){
      //std::cout<<"Hello154"<<std::endl;
      crossover_xy (phmaps, phmaps_next, pos1, pos2, x);      /* Perform crossover with probability */
      //std::cout<<"Hello155"<<std::endl;
    }
    else{
      //std::cout<<"Hello156"<<std::endl;
      std::map<grph::node_t, grph::node_t> phmap1, phmap2;  
      for(std::map<grph::node_t, grph::node_t>::iterator i = phmaps[pos1].begin(), j = phmaps[pos2].begin(); i != phmaps[pos1].end(); ++i, ++j){
        phmap1[i->first] = i->second;
        phmap2[j->first] = j->second;
      }
      phmaps_next.push_back(phmap1);    /* Else copy as it is */
      phmaps_next.push_back(phmap2);    
      //std::cout<<"Hello157"<<std::endl;
    }    	
  }
}

//Perform mutation in the solution at "index" in "gener_1d_next"
void  mutation (std::vector<std::map<grph::node_t, grph::node_t> >& phmaps_next, int index, int x){
  grph::node_t first, second;
  int type = rand() % 6;  // a number between 0 and 5
  switch (type){
    case 0:  // swap two adjacent qubits along x direction
    case 1:{
      first = rand() % x; 
      second = first + 1;
      if (second < x)
        exchange(phmaps_next[index], first, second);
      }break; 

   /* case 2:  // swap an empty slot and a randomly chosen non-empty slot
      if (nvar == x) break;  // No empty slots; skip this
      do
        first = rand() % x; 
      while (phmaps_next[index].find(first) != phmaps_next[index].end());//[first] != EMPTY);
      do
        second = (rand()/(float)RAND_MAX) * x; 
      while (phmaps_next[index].find(second) != phmaps_next[index].end());//[second] == EMPTY);
      exchange(phmaps_next[index], first, second);
      break; */
    case 2: 
    case 3:  // rotate solution right
    case 4:{ 
      std::map<grph::node_t, grph::node_t>::reverse_iterator n = phmaps_next[index].rbegin();
      std::map<grph::node_t, grph::node_t>::reverse_iterator m = n;
      grph::node_t tmp = m->second;
      for(++m; m != phmaps_next[index].rend(); ++n, ++m) //for (int i=1; i<x; i++)
        n->second = m->second; 
      n->second = tmp;
      }break;

    case 5:  // swap two arbitrarily randomly chosen entries
      first = rand() % x;  
      second = rand() % x;  
      exchange(phmaps_next[index], first, second);
      break;
  }
}

//Exange two qubits in the solution of "phmap" at positions "a" and "b".
void exchange(std::map<grph::node_t, grph::node_t>& phmap, grph::node_t a, grph::node_t b){
  grph::node_t n = phmap[a];
  phmap[a] = phmap[b];
  phmap[b] = n;
}

//Perform mutation with probability MUTPROB, and copy back  "phmaps_next" into "phmaps" 
void  mutation_and_copy(std::vector<std::map<grph::node_t, grph::node_t> >& phmaps, std::vector<std::map<grph::node_t, grph::node_t> >& phmaps_next, int x){
  for(int i=NBEST; i<POPSIZE; i++){
    if((rand()/(float)RAND_MAX) < MUTPROB)
      mutation (phmaps_next, i, x);           // Mutate i-th solution in "gener_nxt"

    for(std::map<grph::node_t, grph::node_t>::iterator j = phmaps_next[i].begin(); j != phmaps_next[i].end(); ++j)
      phmaps[i][j->first] = j->second;
  }
}

//GA based ordering of the qubits so as to optimize the cost
//HEXAGONAL ARCHITECTURE mapping
void  ga_search(std::map<grph::node_t, grph::node_t>& phmap, std::pair<int, int> dim, grph::Graph& lg, Metric_Type metric){
  //N (POPSIZE) 1D-grid of size (x) initialized with garbage qubit line 
  std::vector<std::map<grph::node_t, grph::node_t> > phmaps(POPSIZE, std::map<grph::node_t, grph::node_t>( ));
  std::vector<std::map<grph::node_t, grph::node_t> > phmaps_next(POPSIZE, std::map<grph::node_t, grph::node_t>( ));
  std::vector<double> fitness(POPSIZE, 0);
  std::vector<double> percent_fitness(POPSIZE, 0);
  srand(time(NULL));
  //std::cout<<"Hello10"<<std::endl;
  gen_initial_pop(phmaps, phmap, lg);
  //std::cout<<"Before sorting"<<std::endl;
  //display_all_pop(phmaps, dim, lg, metric);
  //std::cout<<"Hello11"<<std::endl;
  //std::cout<<"Original Population"<<std::endl;
  //display_all_pop(phmaps, pg, lg);

  for(int numgen=1; numgen<=MAXGEN; numgen++){ //MAXGEN
    //std::cout<<"Hello12"<<std::endl;
    calculate_fitness(phmaps, fitness, dim, lg, metric);
    //std::cout<<"Hello13"<<std::endl;
    sort_gener_cur(phmaps, fitness);
    //std::cout<<"Before sorting"<<std::endl;
    //display_all_pop(phmaps, dim, lg, metric);
    if (fitness[0] == 0) break;  //minimal solution
    //std::cout<<"Hello14"<<std::endl;
    calculate_percent_fitness(fitness, percent_fitness);
    //std::cout<<"Hello15"<<std::endl;
    generation_next(phmaps, phmaps_next, percent_fitness, lg.getNodes().size());
    //std::cout<<"Hello16"<<std::endl;
    mutation_and_copy(phmaps, phmaps_next, lg.getNodes().size());
    //std::cout<<"Population after mutation and copy"<<std::endl;
    //display_all_pop(phmaps, pg, lg);*/
  }  
  //std::cout<<"Population after "<<MAXGEN<<" generation"<<std::endl;
  //display_all_pop(phmaps, dim, lg, metric);
  phmap.clear();
  if (fitness[0] != 0){
    calculate_fitness(phmaps, fitness, dim, lg, metric);
    sort_gener_cur(phmaps, fitness);
  }
  //std::cout<<"Fitness"<<fitness[0]<<std::endl;
  for(std::map<grph::node_t, grph::node_t>::iterator i = phmaps[0].begin(); i != phmaps[0].end(); ++i)
     phmap[i->first] = i->second;
}


//Display placement details in 1D grid.
void displayMapping(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& lg){
  std::cout<<"#Qubit mapping details: (logical) = Physical"<<std::endl;
  for(std::vector<grph::node_t>::iterator i = lg.getNodes().begin(); i != lg.getNodes().end(); ++i)
    std::cout<<(*i)<<"="<<phmap[*i]<<std::endl;
}

//Display the population and their corresponding cost.
//Hexagonal Architecture
void   display_all_pop(std::vector<std::map<grph::node_t, grph::node_t> > &phmaps, std::pair<int, int> dim, 
                                                                   grph::Graph& lg, Metric_Type metric){
  for(std::vector<std::map<grph::node_t, grph::node_t> >::iterator i = phmaps.begin(); i != phmaps.end(); ++i){
    displayMapping(*i, lg);
    if (metric == Metric_Type::COUPLING_COST)
      std::cout<<"Cost: "<<mappingCost(*i, dim, lg)<<std::endl;
    else
      std::cout<<"Cost: "<<remoteCNOTCost(*i, dim, lg)<<std::endl;    
  }
}

