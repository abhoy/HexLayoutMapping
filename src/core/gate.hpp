/* ========================================================================== */
/*                                                                            */
/*   Gate.hpp                                                                 */
/*   (c) 2014 Author                                                          */
/*                                                                            */
/*   Reversible & quantum gate declaration                                    */
/*                                                                            */
/* ========================================================================== */
/*#ifndef std::string_HPP
#define std::string_HPP
#include<std::string>
#endif /* std::string_HPP */

/*#ifndef COLLECTION_HPP
#define COLLECTION_HPP
#include<vector>
#include<map>
#endif /* COLLECTION_HPP */
#include <vector>
#include <map>
#include <cmath>
#define	PI	3.142857143

#ifndef GATE_HPP
#define GATE_HPP
using namespace std;
namespace rev{
  //CCX => Toffoli gate and MCX => Multicontrol Toffoli gate
  enum {X, Y, Z, H, S, SDG, T, TDG, RX, RY, RZ, U1, U2, U3, CX, CCX, MCX, SWAP}; 
  
  typedef unsigned gate_t;
  typedef unsigned line_t;
  typedef string angle_t;
  class Gate;
  struct GateCompare{
   bool operator() (const Gate& lhs, const Gate& rhs);
  };
  class Gate{
    private:
      static int gate_count;
      int gate_id;
      vector<line_t> controls;
      vector<line_t> targets;
      gate_t type;
      angle_t theta;
      angle_t phi;
      angle_t lamda;
      
    public:
      Gate(gate_t type, line_t );
      //Gate(gate_t type, angle_t, line_t );
      //Gate(gate_t type, angle_t, angle_t, line_t );
      //Gate(gate_t type, angle_t, angle_t, angle_t, line_t );
      Gate(gate_t , line_t , line_t );
      Gate(std::string gate, std::map<string, line_t>& var_indices); //QASM type
      Gate(gate_t type, std::vector<line_t>& controls, line_t target); //Real type
      void display(ostream & os);
      void display2(ostream & os);
      gate_t getType();
      //angle_t getAngle();
      int getId() const;
      vector<line_t> getControls();
      vector<line_t> getTargets();
      void updateLines(line_t control, line_t target);
      void updateLines(line_t target);
      string removeZero(string str);
      double round_sig(double num, int sig);
      //void mergeGate(Gate g); //Need to redefine due angle type changed from double to string
      void makeEquivalent(); //Need to redefine due angle type changed from double to string
      //bool isIdentity();
  };
}
#endif /* GATE_HPP */
