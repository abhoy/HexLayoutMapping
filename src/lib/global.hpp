#include <iostream>
#include <vector>
#include <map>
#include "optimize.hpp"
#include "../core/graph.hpp"

#define EMPTY -1  // Default logical to physical mapping 
#define	POPSIZE	30 // Population size 30
#define	NBEST 4	// Number of best solutions to carry; must be even
#define CROSSPROB 0.70 // Probability of crossover
#define	MUTPROB	0.10	// Probability of mutation 
#define	MAXGEN 500	// Number of generations to run 500

#ifndef GLOBAL_HPP
#define GLOBAL_HPP


void gen_initial_pop(std::vector<std::vector<rev::line_t> > &phmaps, std::map<grph::node_t, grph::node_t> & phmap, grph::Graph& lg); //grph::Graph& pg,

void physicalMap(std::map<grph::node_t, grph::node_t> & phmap_new, std::map<grph::node_t, grph::node_t> & phmap, grph::Graph& lg); //grph::Graph& pg,


void calculate_fitness(std::vector<std::map<grph::node_t, grph::node_t> > & phmaps, std::vector<double>& fitness,  
                                                     std::pair<int, int> dim, grph::Graph& lg, Metric_Type metric); //hexagon

void calculate_percent_fitness(std::vector<double>& fitness, std::vector<double>& percent_fitness);

void sort_gener_cur(std::vector<std::map<grph::node_t, grph::node_t> > & phmaps, std::vector<double>& fitness);

int  selection(std::vector<double>& percent_fitness);

void crossover_xy (std::vector<std::map<grph::node_t, grph::node_t> > & phmaps, 
                    std::vector<std::map<grph::node_t, grph::node_t> > & phmaps_next, int pos1, int pos2, int x);

void generation_next(std::vector<std::map<grph::node_t, grph::node_t> > & phmaps, 
                    std::vector<std::map<grph::node_t, grph::node_t> > & phmaps_next, std::vector<double>& percent_fitness, int x);


void  mutation (std::vector<std::map<grph::node_t, grph::node_t> >& phmaps_next, int index, int x);

void exchange(std::map<grph::node_t, grph::node_t>& phmap, grph::node_t a, grph::node_t b);

void  mutation_and_copy(std::vector<std::map<grph::node_t, grph::node_t> >& phmaps, std::vector<std::map<grph::node_t, grph::node_t> >& phmaps_next, int x);

void  ga_search(std::map<grph::node_t, grph::node_t>& phmap, std::pair<int, int> dim, grph::Graph& lg, Metric_Type metric);//hexagon


void displayMapping(std::map<grph::node_t, grph::node_t>& phmap, grph::Graph& lg);

#endif /* GLOBAL_HPP */
