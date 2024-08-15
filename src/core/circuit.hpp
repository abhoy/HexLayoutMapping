/* ========================================================================== */
/*                                                                            */
/*   Circuit.hpp                                                              */
/*   (c) 2014 Author                                                          */
/*                                                                            */
/*   Reversible circuit class declaration                                     */
/*                                                                            */
/* ========================================================================== */
//#include<string>
//#include<vector>
#include "gate.hpp"
#ifndef CIRCUIT_HPP
#define CIRCUIT_HPP
using namespace std;
namespace rev{
  typedef unsigned window_t;
  //extern window_t window;
  int convertToInteger(string );
  void exchange(vector<line_t>& lines, int i, int j);
  class Circuit{
    private:
      //string version;
      int numvars;
      //vector<string> variables;
      vector<string> variables; //backward compitable
      vector<line_t> lines;
      vector<Gate> gates;
      map<line_t, string> var_names; 
      map<string, line_t> var_indices;
      map<line_t, double> CNOTDepth; 
    public:
      Circuit();
      Circuit(const Circuit & ckt);
      ~Circuit();
      void init();
      void readQASM(string);
      void readReal(string);
      void writeQASM(string, unsigned);
      void writeReal(string file);
      void displayQASM(ostream & os, unsigned size);
      void displayQisKit(ostream & os, unsigned size);
      vector<Gate> getGates();
      vector<line_t> getLines();
      void setGates(vector<Gate>& gates);
      void setLines(vector<line_t>& lines);
      unsigned twoQubitGateCount();
      unsigned getCNOTDepth();
      unsigned maxCNOTDepth();
      void reduceCNOTDepth(Gate& g);
      unsigned twoQubitGateCount(unsigned pos);
      void addQubit(string );
      void clear();
      void operator = (const Circuit & ckt); 
      friend void exchange(vector<line_t>& lines, int i, int j);
      friend int convertToInteger(string );
  };
}
#endif /* CIRCUIT_HPP */
