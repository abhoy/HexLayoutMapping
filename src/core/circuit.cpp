/* ========================================================================== */
/*                                                                            */
/*   Circuit.cpp                                                              */
/*   (c) 2014 Author                                                          */
/*                                                                            */
/*   Reversible circuit class definition                                      */
/*                                                                            */
/* ========================================================================== */

#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <vector>
#include <map>
#include "circuit.hpp"
using namespace rev;

/* Interchange the position of i-th qubit with j-th qubit*/
void rev::exchange(vector<line_t>& lines, int i, int j){
  line_t line = lines[i];
  lines[i] = lines[j];
  lines[j] = line;
}
/* Display the Quantum circuit in QASM format*/
void Circuit::displayQASM(ostream & os, unsigned size){
  //os<<"#QASM Format: To run in IBM Qiskit environment uncomment the following two lines and the 4th from last line"<<std::endl;
  //os<<"#from qiskit import QuantumCircuit"<<std::endl;
  //os<<"#qasm_str = \"\"\""<<std::endl;
  os<<"OPENQASM 2.0;"<<std::endl;
  os<<"include \"qelib1.inc\";"<<std::endl;
  os<<"qreg q["<<size<<"];"<<std::endl;
  os<<"creg c["<<size<<"];"<<std::endl;
  for(vector<Gate>::iterator i=this->gates.begin(); i!=this->gates.end(); ++i)
    (*i).display(os);
  //os<<"#\"\"\""<<std::endl;
  //os<<"qc = QuantumCircuit.from_qasm_str(qasm_str)"<<std::endl;
  //os<<"qc.draw()"<<std::endl;
  //os<<"#End of File"<<std::endl;
}
/* Display the Quantum circuit in IBM QisKit format*/
void Circuit::displayQisKit(ostream & os, unsigned size){
  os<<"#QISKIT Format"<<std::endl;
  os<<"from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister, Aer"<<std::endl;
  os<<"from qiskit.compiler import transpile"<<std::endl;
  //os<<"q = QuantumRegister("<<size<<")"<<std::endl;
  //os<<"c = ClassicalRegister("<<size<<")"<<std::endl;
  os<<"qc = QuantumCircuit("<<size<<", "<<size<<")"<<std::endl;
  for(vector<Gate>::iterator i=this->gates.begin(); i!=this->gates.end(); ++i)
    (*i).display2(os);
  os<<"#qc.draw(output='mpl')"<<std::endl;
  os<<"backend = Aer.get_backend('qasm_simulator')"<<std::endl;
  os<<"trans_qc = transpile(qc, backend, basis_gates=['id','cx', 'u3'], optimization_level=0)"<<std::endl;
  os<<"#trans_qc.draw(output='mpl')"<<std::endl;
  os<<"trans_qc.qasm(formatted=True)"<<std::endl;
  os<<"#End of File"<<std::endl;
}
void Circuit::readQASM(string file){  
  //std::cout<<file.c_str()<<std::endl;
  string line;
  ifstream myfile(file.c_str());
  bool readGate = false;
  if(myfile.is_open()){
    while(getline(myfile,line)){ 
      if(0 == line.size() || 1 == line.size()) continue; //skip empty lines
      if('#' == line[0])  continue; //skip comments  
      istringstream iss(line);
      string command;
      iss >> command;
      //std::cout << command<<std::endl;
      if("OPENQASM" == command) continue; //Version
      else if("include" == command) continue; //header file
      else if("qreg" == command) continue; /*{//Variables q[16]
        iss >> command;
        this->numvars = convertToInteger(command.substr(2, command.length()-4));
        for (unsigned u = 0; u < this->numvars; ++u){
          ostringstream oss;
          this->lines.push_back(u);
          oss << "q["<<u<<"]";
	  this->var_names[u] = oss.str();
	  this->var_indices[oss.str()] = u;
        }
      }*/
      else if("creg" == command) continue; //skip classical registers
      else if ("barrier" == command) continue;      //skip
      else if ("measure" == command) continue;      //skip
      else{
        iss >> command;
        while(iss){
          //std::cout<<"command "<<command<<":"<<std::endl;
          //if(command.find(" ") != string::npos)
           // command.replace(command.find(" "), 0, "");
          //std::cout<<"command "<<command<<":"<<std::endl;
          auto pos1 = string::npos, pos2 = command.find(",");
          //int i =0;
          for ( ; pos2 != string::npos; ){       
            string qubit = pos1 != string::npos ? command.substr(pos1 + 1, pos2 - pos1 - 1) : command.substr(0, pos2); 
            addQubit(qubit);
            //std::cout<<"qubit1 :"<<qubit<<std::endl;
            //if(qubits.find(" ") != string::npos)
              //qubits.replace(qubits.find(" "), 0, "");
            //i++; if ( i == 4) break;
            pos1 = pos2, pos2 = command.find(",", pos1 + 1);  
            //std::cout<<"pos2"<<pos2<<"npos"<<string::npos<<std::endl;
            if (pos2 == string::npos){ 
              if ((pos2 = command.find(";")) != string::npos) break;            
              iss >> command;
              pos1 = string::npos, pos2 = command.find(",");
              if (pos2 == string::npos && (pos2 = command.find(";")) != string::npos) break;            
            }            
            //std::cout<<"pos2"<<pos2<<"npos"<<string::npos<<std::endl;
          }          
          if (pos2 == string::npos) pos2 = command.find(";");
          if (pos2 != string::npos) {
            string qubit = pos1 != string::npos ? command.substr(pos1 + 1, pos2 - pos1 - 1) : command.substr(0, pos2); 
            //std::cout<<"qubit2 :"<<qubit<<std::endl;
            addQubit(qubit);          
          }          
          iss >> command;
        }
        //std::cout<<"Hello"<<std::endl;
        this->gates.push_back(Gate(line, this->var_indices));
        //std::cout<<"Hello"<<std::endl;
        /*int epos = command.find(","); //multi qubit gate with qubits separated by comma (,)
        if (epos != std::string::npos){
	  int spos = 0; // first qubit position
          do{
          	string qubit = command.substr(spos,epos-spos); //command.find(",")
                std::cout<<"spos: "<<spos<<" epos:"<<epos<<std::endl;
                std::cout<<"qubit c:"<<qubit<<std::endl;
          	addQubit(qubit);
		spos = epos + 1;
                epos = command.find(",", spos);
          }while(epos != std::string::npos);
          std::cout<<"spos: "<<spos<<" epos:"<<epos<<std::endl;
          qubit = command.substr(spos, command.length()-spos-1); 
          std::cout<<"qubit t:"<<qubit<<std::endl;
	  addQubit(qubit);
        }
        else{
          string qubit = command.substr(0, command.length()-1); 
          std::cout<<"qubit t:"<<qubit<<std::endl;
          addQubit(qubit);
        }    
	*/
      }
    } 
  }
  else
    cout << "File can not be opened." <<endl;
}

void Circuit::readReal(string file){  
  //std::cout<<file.c_str()<<std::endl;
  string line;
  ifstream myfile(file.c_str());
  bool readGate = false;
  if(myfile.is_open()){
    int qubits = 0;
    while(getline(myfile,line)){ 
      if(0 == line.size() || 1 == line.size()) continue; //skip empty lines
      if('#' == line[0])  continue; //skip comments  
      istringstream iss(line);
      string command;
      iss >> command;
      //std::cout << command<<std::endl;
      if ( command == ".version" ) continue; //Version
      else if ( command == ".mode" ) continue;
      else if ( command == ".numvars" ) {
        iss >> command;
        qubits = std::stoi(command);
        continue;
      }
      else if ( command == ".variables" ) {
        while( iss >> command) addQubit(command);        
        continue;
      }
      else if ( command == ".inputs" ) continue;
      else if ( command == ".outputs" ) continue;
      else if ( command == ".constants" ) continue;
      else if ( command == ".garbage" ) continue;
      else if ( command == ".inputbus" ) continue;
      else if ( command == ".outputbus" ) continue;
      else if ( command == ".state" ) continue;
      else if ( command == ".module" ) continue;
      else if ( command == ".begin" ) continue;
      else if ( command == ".end" ) continue;
      else if ( command == ".define" ) {
        do{
          getline(myfile, line);
          istringstream iss(line);
          iss >> command;
        }while ( command != ".enddefine" );
        continue;
      }
      else {
        //std::cout <<"Toffoli"<<std::endl;            
        if (command[0] == 't') {   //Toffoli Gate
          int count = std::stoi(command.substr(1));
          if (count == 1 ) {
            iss >> command;
            this->gates.push_back(Gate(rev::X, this->var_indices[command]));
          }
          else if (count == 2) {	
            string control, target;
            iss >> control;
            iss >> target;
            this->gates.push_back(Gate(rev::CX, this->var_indices[control], this->var_indices[target]));
          }
          else {	            
            //std::cout <<"size:"<<count<<std::endl;
            std::vector<line_t> controls;
            string control, target;
            for (auto i = 0; i < count - 1; i++){
              iss >> control;
              controls.push_back(this->var_indices[control]);
            }
            iss >> target;
            rev::gate_t type;
            if (count == 3) type = rev::CCX;
            else type = rev::MCX;
            this->gates.push_back(Gate(type, controls, this->var_indices[target]));
          }    	
        }
        else if (command[0] == 'f') {  //Fredkin Gate                 
          int count = std::stoi(command.substr(1));
          if (count == 2){
            string target1, target2;
            iss >> target1;
            iss >> target2;
            this->gates.push_back(Gate(rev::SWAP, this->var_indices[target1], this->var_indices[target2]));
          }
          else {
            std::vector<line_t> controls;
            string control, target;
            for (auto i = 0; i < count - 1; i++){
              iss >> control;
              controls.push_back(this->var_indices[control]);
            }
            iss >> target;
            rev::gate_t type;
            if (count == 3) type = rev::CCX;
            else type = rev::MCX;
            this->gates.push_back(Gate(rev::CX, this->var_indices[target], controls[count - 2]));
            this->gates.push_back(Gate(type, controls, this->var_indices[target]));            
            this->gates.push_back(Gate(rev::CX, this->var_indices[target], controls[count - 2]));
          }
        }
        else if (command == "p") { //Peres Gate
          //int count = std::stoi(command.substr(1));
          //if (count == 3){
          std::vector<line_t> controls;
          string control, target;
          for (auto i = 0; i < 2; i++){
            iss >> control;
            controls.push_back(this->var_indices[control]);
          }
          iss >> target;
          this->gates.push_back(Gate(rev::CCX, controls, this->var_indices[target]));            
          this->gates.push_back(Gate(rev::CX, controls[0], controls[1]));                   
        }
        else if (command == "pi") { //Peres Gate
          //int count = std::stoi(command.substr(1));
          //if (count == 3){
          std::vector<line_t> controls;
          string control, target;
          for (auto i = 0; i < 2; i++){
            iss >> control;
            controls.push_back(this->var_indices[control]);
          }
          iss >> target;          
          this->gates.push_back(Gate(rev::CX, controls[0], controls[1]));                   
          this->gates.push_back(Gate(rev::CCX, controls, this->var_indices[target]));            
        }
        else if (command == "v") { //V Gate
          string control, target;
          iss >> control;
          iss >> target;
          this->gates.push_back(Gate(rev::H, this->var_indices[target]));            
          this->gates.push_back(Gate(rev::T, this->var_indices[control]));
          this->gates.push_back(Gate(rev::CX, this->var_indices[target], this->var_indices[control]));
          this->gates.push_back(Gate(rev::T, this->var_indices[target]));            
          this->gates.push_back(Gate(rev::TDG, this->var_indices[control]));
          this->gates.push_back(Gate(rev::CX, this->var_indices[target], this->var_indices[control]));
          this->gates.push_back(Gate(rev::H, this->var_indices[target]));          
        } 
        else if (command == "v+") { //V Gate
          string control, target;
          iss >> control;
          iss >> target;
          this->gates.push_back(Gate(rev::H, this->var_indices[target]));            
          this->gates.push_back(Gate(rev::TDG, this->var_indices[control]));
          this->gates.push_back(Gate(rev::CX, this->var_indices[target], this->var_indices[control]));
          this->gates.push_back(Gate(rev::TDG, this->var_indices[target]));            
          this->gates.push_back(Gate(rev::T, this->var_indices[control]));
          this->gates.push_back(Gate(rev::CX, this->var_indices[target], this->var_indices[control]));
          this->gates.push_back(Gate(rev::H, this->var_indices[target]));          
        }   
        else cout<<"Invalid gate: "<<command<<endl;              
      }
    } 
  }
  else
    cout << "File can not be opened." <<endl;
}

void Circuit::addQubit(string qubit){
  if (this->var_indices.find(qubit) == this->var_indices.end()){
    this->var_names[this->lines.size()] = qubit;
    this->var_indices[qubit] = this->lines.size();
    this->lines.push_back(this->lines.size());	
  }
  //std::cout<<"hello:"<<var_indices[qubit]<<" qubit"<<qubit<<std::endl;
}
unsigned Circuit::twoQubitGateCount(){
  unsigned count = 0;
  for(std::vector<rev::Gate>::iterator i=this->gates.begin(); i!=this->gates.end(); ++i)
    if(CX == i->getType()) count++;
  return count;
}
unsigned Circuit::getCNOTDepth(){
  this->CNOTDepth = std::map<rev::line_t, double> ();
  for(std::vector<rev::line_t>::iterator l=this->lines.begin(); l!=this->lines.end(); ++l){
    this->CNOTDepth[*l] = 0;
  }
  for(std::vector<rev::Gate>::iterator g=this->gates.begin(); g!=this->gates.end(); ++g){
    if(CX == g->getType()){
      double cdepth = this->CNOTDepth[g->getControls()[0]] > this->CNOTDepth[g->getTargets()[0]] ? 
                                     this->CNOTDepth[g->getControls()[0]] + 1 : this->CNOTDepth[g->getTargets()[0]] + 1;
      this->CNOTDepth[g->getControls()[0]] = cdepth;
      this->CNOTDepth[g->getTargets()[0]] = cdepth; 
                                   
    }
  }
  return maxCNOTDepth();
}

unsigned Circuit::maxCNOTDepth(){
  double maxDepth = 0;
  for(std::vector<rev::line_t>::iterator l=this->lines.begin(); l!=this->lines.end(); ++l){
    if (maxDepth < this->CNOTDepth[*l]) maxDepth = CNOTDepth[*l];
  }
  return maxDepth;
}
void Circuit::reduceCNOTDepth(Gate& g){
  //g.display(std::cout);
  double cdepth = this->CNOTDepth[g.getControls()[0]] < this->CNOTDepth[g.getTargets()[0]] ? 
                                     this->CNOTDepth[g.getControls()[0]] - 1 : this->CNOTDepth[g.getTargets()[0]] - 1;
  this->CNOTDepth[g.getControls()[0]] = cdepth;
  this->CNOTDepth[g.getTargets()[0]] = cdepth; 
  std::cout<<"Cnot Depth";
  for(std::vector<rev::line_t>::iterator l=this->lines.begin(); l!=this->lines.end(); ++l){
    std::cout<<this->CNOTDepth[*l]<<" ";
  }
  std::cout<<std::endl;
}
unsigned Circuit::twoQubitGateCount(unsigned pos){
  unsigned count = 0;
  for(std::vector<rev::Gate>::iterator i=this->gates.begin() + pos; i!=this->gates.end(); ++i)
    if(CX == i->getType()) count++;
  return count;
}

void Circuit::writeQASM(string file, unsigned size){
  char fileName[50];
  int i, ln = file.length();
  for(i=0; i < ln; i++)
    fileName[i] = file[i];
  fileName[i] = '\0';
  ofstream myfile(fileName);
  std::map<rev::line_t, string> vars;

  if(myfile.is_open()){
    myfile<<"OPENQASM 2.0;"<<std::endl;
    myfile<<"include \"qelib1.inc\";"<<std::endl;
    myfile<<"qreg q["<<size<<"];"<<std::endl;
    myfile<<"creg q["<<size<<"];"<<std::endl;
    for(vector<Gate>::iterator i=this->gates.begin(); i!=this->gates.end(); ++i)
      (*i).display(myfile);
    myfile.close();
  }
  else
    cout << "File can not be opened." <<endl;
}

void Circuit::writeReal(string file){
  char fileName[50];
  int i, ln = file.length();
  for(i=0; i < ln; i++)
    fileName[i] = file[i];
  fileName[i] = '\0';
  ofstream myfile(fileName);
  std::map<rev::line_t, string> vars;	
  variables.clear();
  int vcount = 0, gcount = 0;
  for(vector<line_t>::iterator i=lines.begin(); i!=lines.end(); ++i){
    std::ostringstream vrname;
    vrname<<"x"<<vcount++;
    vars[*i]=vrname.str();    
    //std::cout<<vrname.str()<<std::endl;
  }
  
  if(myfile.is_open()){
    //cout<<"hello"<<endl;
    myfile<<".version 1.0"<<endl;
    myfile<<".numvars "<<lines.size()<<endl;
    myfile<<".variables";
    //std::cout<<"Hello1"<<std::endl;
    for(vector<line_t>::iterator i=lines.begin(); i!=lines.end(); ++i)
        myfile<<" "<<vars[(*i)];
    myfile<<endl;
    myfile<<".begin"<<endl;
    for(vector<Gate>::iterator i=gates.begin(); i!=gates.end(); ++i){
      vector<line_t> controls = (*i).getControls();
      vector<line_t> targets = (*i).getTargets();
      //{X, Y, Z, H, S, SDG, T, TDG, RX, RY, RZ, U1, U2, U3, CX, CCX, MCX, SWAP}; 
      if(X == (*i).getType())
        myfile<<"x";
      else if(Y == (*i).getType())
        myfile<<"y"; 
      else if(Z == (*i).getType())
        myfile<<"z";
      else if(H == (*i).getType())
        myfile<<"h";
      else if(S == (*i).getType())
        myfile<<"s";
      else if(SDG == (*i).getType())
        myfile<<"s+";
      else if(T == (*i).getType())
        myfile<<"p";
      else if(TDG == (*i).getType())
        myfile<<"p+";
      else if(RX == (*i).getType())
        myfile<<"rx";
      else if(RY == (*i).getType())
        myfile<<"ry";
      else if(RZ == (*i).getType())
        myfile<<"rz";
      else if(U1 == (*i).getType())
        myfile<<"u1";
      else if(U2 == (*i).getType())
        myfile<<"u2";
      else if(U3 == (*i).getType())
        myfile<<"u3";
      else if(CX == (*i).getType() || CCX == (*i).getType() || MCX == (*i).getType())
        myfile<<"t";
      else if(SWAP == (*i).getType())
        myfile<<"f"; 
      else{
        myfile<<"unknown";
        return;
      }
      if(CX == (*i).getType() || CCX == (*i).getType() || MCX == (*i).getType() || SWAP == (*i).getType())
        myfile<<controls.size()+targets.size();
      for(vector<line_t>::iterator i=controls.begin(); i!=controls.end(); ++i)
        myfile<<" "<<vars[(*i)];
      for(vector<line_t>::iterator i=targets.begin(); i!=targets.end(); ++i)
        myfile<<" "<<vars[(*i)];
      myfile<<endl;
    }
    //std::cout<<"Hello3"<<std::endl;
    myfile<<".end"<<endl;
    myfile.close();
  }
  else
    cout << "File can not be opened." <<endl;
}

vector<Gate> Circuit::getGates(){
  return gates;
}
vector<line_t> Circuit::getLines(){
  return lines;
}

void Circuit::clear(){
  this->gates.clear();
  this->lines.clear();
}

void Circuit::setGates(vector<Gate>& ngates){
  this->gates.clear();
  for(vector<Gate>::iterator g=ngates.begin(); g != ngates.end(); ++g)
    this->gates.push_back(*g);
}
void Circuit::setLines(vector<line_t>& nlines){
  this->lines.clear();
  for(vector<line_t>::iterator l=nlines.begin(); l != nlines.end(); ++l)
    this->lines.push_back(*l);
  
}
void Circuit::init(){
  this->lines = vector<line_t> ();
  this->gates = vector<Gate> ();;
  this->var_names = map<line_t, string> (); 
  this->var_indices = map<string, line_t> (); 
}
Circuit::Circuit(){
 init();
}
Circuit::Circuit(const Circuit & ckt){
  init();
  for (vector<line_t>::const_iterator l = ckt.lines.begin(); l != ckt.lines.end(); ++l)
    this->lines.push_back(*l);
  for (vector<Gate>::const_iterator g = ckt.gates.begin(); g != ckt.gates.end(); ++g)
    this->gates.push_back(*g);
  for (map<line_t, string>::const_iterator vn = ckt.var_names.begin(); vn != ckt.var_names.end(); ++vn)
    this->var_names[vn->first] = vn->second;
  for (map<string, line_t>::const_iterator vi = ckt.var_indices.begin(); vi != ckt.var_indices.end(); ++vi)
    this->var_indices[vi->first] = vi->second;
  for (map<rev::line_t, double>::const_iterator cd = ckt.CNOTDepth.begin(); cd != ckt.CNOTDepth.end(); ++cd)
    this->CNOTDepth[cd->first] = cd->second;
}
void Circuit::operator = (const Circuit & ckt){
  init();
  for (vector<line_t>::const_iterator l = ckt.lines.begin(); l != ckt.lines.end(); ++l)
    this->lines.push_back(*l);
  for (vector<Gate>::const_iterator g = ckt.gates.begin(); g != ckt.gates.end(); ++g)
    this->gates.push_back(*g);
  for (map<line_t, string>::const_iterator vn = ckt.var_names.begin(); vn != ckt.var_names.end(); ++vn)
    this->var_names[vn->first] = vn->second;
  for (map<string, line_t>::const_iterator vi = ckt.var_indices.begin(); vi != ckt.var_indices.end(); ++vi)
    this->var_indices[vi->first] = vi->second;
  for (map<rev::line_t, double>::const_iterator cd = ckt.CNOTDepth.begin(); cd != ckt.CNOTDepth.end(); ++cd)
    this->CNOTDepth[cd->first] = cd->second;
}
Circuit::~Circuit(){
}
int rev::convertToInteger(string value){
  int i=0, result=0, len=value.length();
  while(' ' == value[i]) i++; //skipping initial spaces
  for(; i<len; i++){
    if(48 <= value[i] && 57 >= value[i])
      result = 10*result + value[i]-48;
  }
  return result;
}

