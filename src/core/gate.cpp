/* ========================================================================== */
/*                                                                            */
/*   Gate.cpp                                                                 */
/*   (c) 2014 Author                                                          */
/*                                                                            */
/*   Reversible & quantum gate definition                                     */
/*                                                                            */
/* ========================================================================== */
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include "gate.hpp"
using namespace rev;
int Gate::gate_count = 0;

string Gate::removeZero(string str){ 
    int i = str.length(); 
    while (str[i] == '0') 
       i--; 
  
    str.erase(i); 
    return str; 
}
 
double Gate::round_sig(double num, int sig){
  double result = static_cast<double>(static_cast<int>(num * pow(10, sig))) / pow(10, sig);
  return result;
}

bool GateCompare::operator() (const Gate& lhs, const Gate& rhs){
  return lhs.getId() < rhs.getId();
}

void Gate::updateLines(line_t control, line_t target){
  controls[0] = control;
  targets[0] = target;
}

void Gate::updateLines(line_t target){
  targets[0] = target;
}

vector<line_t> Gate::getControls(){
  return controls;
}

vector<line_t> Gate::getTargets(){
  return targets;
}

gate_t Gate::getType(){
  return type;
}

int Gate::getId() const{
  return gate_id;
}

Gate::Gate(gate_t type, line_t target){
  gate_id = ++gate_count;
  targets.push_back(target);
  this->type = type;
}
/*Gate::Gate(gate_t type, angle_t angle, line_t target){
 gate_id = ++gate_count;
 targets.push_back(target);
 this->type = type;
 this->angle1 = 0.;
 this->angle2 = 0.;
 this->angle3 = angle;
}
Gate::Gate(gate_t type, angle_t angle2, angle_t angle3, line_t target){
 gate_id = ++gate_count;
 targets.push_back(target);
 this->type = type;
 this->angle1 = PI/2;
 this->angle2 = angle2;
 this->angle3 = angle3;
}
Gate::Gate(gate_t type, angle_t angle1, angle_t angle2, angle_t angle3, line_t target){
 gate_id = ++gate_count;
 targets.push_back(target);
 this->type = type;
 this->angle1 = angle1;
 this->angle2 = angle2;
 this->angle3 = angle3;
}*/
void Gate::makeEquivalent(){
 switch(this->type){
   case U1:
   case U2:
   case U3:
   case RZ:
     break;
   case X:{
     this->type = U3;
     this->theta = "PI";
     this->phi = "0";
     this->lamda = "PI";
     break;
   }    
   case Y:{
     this->type = U3;
     this->theta = "PI";
     this->phi = "PI/2";
     this->lamda = "PI/2";
     break;
   }    
   case Z:{
     this->type = U1;
     this->theta = "0";
     this->phi = "0";
     this->lamda = "PI";
     break;
   }     
   case H:{
     this->type = U2;
     this->theta = "PI/2";
     this->phi = "0";
     this->lamda = "PI";
     break;
   }       
   case S:{
     this->type = U1;
     this->theta = "0";
     this->phi = "0";
     this->lamda = "PI/2";
     break;
   }
   case SDG:{
     this->type = U1;
     this->theta = "0";
     this->phi = "0";
     this->lamda = "-PI/2";
     break;
   }       
   case T:{
     this->type = U1;
     this->theta = "0";
     this->phi = "0";
     this->lamda = "PI/4";
     break;
   }
   case TDG:{
     this->type = U1;
     this->theta = "0";
     this->phi = "0";
     this->lamda = "-PI/4";
     break;
   }
   default:{
     std::cout<<"INVALID GATE TYPE"<<std::endl;
   }
 }
}
/*void Gate::mergeGate(Gate g){ // g x this
  this->makeEquivalent();
  g.makeEquivalent();
  if (this->type == U1 && g.type == U1) {
    this->lamda += g.lamda;
  }
  else if (this->type == U1 && (g.type == U2 || g.type == U3)) {
    this->type = g.type;
    this->theta = g.theta;
    this->phi = g.phi;
    this->lamda += g.lamda;
  }
  else if ((this->type == U2 || this->type == U3) && g.type == U1){
    this->phi += g.lamda;
  }
  else if ((this->type == U2 || this->type == U3) && (g.type == U2 || g.type == U3)){
    this->type = U3;
    if (round_sig(g.lamda + this->phi, 2) == round_sig(2 * PI, 3) || round_sig(g.lamda + this->phi, 2) == 0.){
      this->theta += g.theta;
      this->phi = g.phi;
    }
    else{
      this->theta -= g.theta;
      this->phi = g.phi + PI;
    }
  }
  if (round_sig(this->theta, 2) >= round_sig(2 * PI, 3)) this->theta -= (2 * PI);
  if (round_sig(this->phi, 2) >= round_sig(2 * PI, 3)) this->phi -= (2 * PI);
  if (round_sig(this->lamda, 2) >= round_sig(2 * PI, 3)) this->lamda -= (2 * PI);
}*/
/*bool Gate::isIdentity(){
  return (round_sig(theta, 2) == 0. || round_sig(theta, 2) == round_sig(2 * PI, 2)) && 
         (round_sig(phi + lamda, 2) == 0. || round_sig(phi + lamda, 2) == round_sig(2 * PI, 2));
}*/
Gate::Gate(gate_t type,line_t control , line_t target){
  gate_id = ++gate_count;
  if(rev::SWAP == type)
    targets.push_back(control);
  else
    controls.push_back(control);
  targets.push_back(target);
  this->type = type;
}

Gate::Gate(std::string gate, std::map<string, line_t>& var_indices){
  gate_id = ++gate_count;
  gate=gate.substr(gate.find_first_not_of(" \n\r\t\f\v"));
  
  stringstream check(gate);      
  string gate_type;
  check>>gate_type;
  //std::cout<<"Gate :"<<gate<<std::endl;
  //std::cout<<"gate_type :"<<gate_type<<std::endl;

  if (gate_type == "x")        type = X;	
  else if (gate_type == "y")   type = Y;
  else if (gate_type == "z")   type = Z;
  else if (gate_type == "h")   type = H;
  else if (gate_type == "s")   type = S;    
  else if (gate_type == "sdg") type = SDG;
  else if (gate_type == "t")   type = T;
  else if (gate_type == "tdg") type = TDG;
  else if (gate_type == "cx")  type = CX;
  else if (gate_type == "ccx")  type = CCX;
  else if (gate_type == "mcx")  type = MCX;
  else if (gate_type == "swap")  type = SWAP;
  else if (gate_type.substr(0,2) == "rz"){
    type = U1;//RZ;
    theta = "0.";
    phi = "0.";
    lamda = gate_type.substr(3,gate_type.length()-4);
    //lamda = std::atof(gate_type.substr(3,gate_type.length()-4).c_str());
    //cout<<gate_type.substr(3,gate_type.length()-4)<<" "<<angle3<<endl;
  }
  else if (gate_type.substr(0,2) == "u1"){
    type = U1;
    theta = "0.";
    phi = "0.";
    lamda = gate_type.substr(3,gate_type.length()-4);
    //lamda = std::atof(gate_type.substr(3,gate_type.length()-4).c_str());
    //cout<<gate_type.substr(3,gate_type.length()-4)<<" "<<angle3<<endl;
  }
  else if (gate_type.substr(0,2) == "u2"){
    type = U2;
    theta = "PI/2";
    string angles = gate_type.substr(3,gate_type.length()-4);
    size_t start, end = 0;
    start = angles.find_first_not_of(',', end);
    end = angles.find(',', start);
    phi = angles.substr(start, end - start);
    //phi = std::atof(angles.substr(start, end - start).c_str());
    start = angles.find_first_not_of(',', end);
    end = angles.find(',', start);
    lamda = angles.substr(start, end - start);    
    //lamda = std::atof(angles.substr(start, end - start).c_str());    
    //cout<<gate_type.substr(3,gate_type.length()-4)<<" "<<angle3<<endl;
  }
  else if (gate_type.substr(0,2) == "u3"){
    type = U3;
    string angles = gate_type.substr(3,gate_type.length()-4);
    size_t start, end = 0;
    start = angles.find_first_not_of(',', end);
    end = angles.find(',', start);
    theta = angles.substr(start, end - start);
    //theta = std::atof(angles.substr(start, end - start).c_str());
    start = angles.find_first_not_of(',', end);
    end = angles.find(',', start);
    phi = angles.substr(start, end - start);
    //phi = std::atof(angles.substr(start, end - start).c_str());
    start = angles.find_first_not_of(',', end);
    end = angles.find(',', start);
    lamda = angles.substr(start, end - start);
    //lamda = std::atof(angles.substr(start, end - start).c_str());
    //cout<<" "<<theta<<" "<<phi<<" "<<lamda<<endl;
    //cout<<gate_type.substr(3,gate_type.length()-4)<<" "<<angle3<<endl;
  }
  else cout<<"Invalid gate: "<<gate_type<<endl;
  

  string qubits;
  check >> qubits;
  //std::cout<<"Qubit :"<<qubits<<std::endl;  
  if(CX == type || CCX == type || MCX == type || SWAP == type){    
    for (auto pos1 = string::npos, pos2 = qubits.find(","); pos2 != string::npos;){       
      string qubit = pos1 != string::npos ? qubits.substr(pos1 + 1, pos2 - pos1 - 1) : qubits.substr(0, pos2); 
      if (SWAP == type) 
        targets.push_back(var_indices[qubit]);
      else 
        controls.push_back(var_indices[qubit]);
       //std::cout<<"qubit1 :"<<var_indices[qubit]<<std::endl;
      //if(qubits.find(" ") != string::npos)
        //qubits.replace(qubits.find(" "), 0, "");
      pos1 = pos2, pos2 = qubits.find(",", pos1 + 1);        
      if (pos2 == string::npos) {          
        if ( (pos2 = qubits.find(";")) == string::npos){            
          check >> qubits;
          pos1 = string::npos, pos2 = qubits.find(",");
          if (pos2 == string::npos) { 
            qubit = qubits.substr(0, qubits.find(";")); 
            //std::cout<<"qubit2 :"<<var_indices[qubit]<<std::endl;
            targets.push_back(var_indices[qubit]);
            break;
          }
        }
        else{
          qubit = qubits.substr(pos1 + 1, pos2 - pos1 - 1); 
          //std::cout<<"qubit2 :"<<var_indices[qubit]<<std::endl;
          targets.push_back(var_indices[qubit]);  
          break;
        }
      }
    }    
    //std::cout<<"C"<<controls[0]<<" T"<<targets[0]<<std::endl;
  }
  else { //arbitrary single qubit gate
    string qubit = qubits.substr(0, qubits.find(";"));
    //std::cout<<"qubit3 :"<<var_indices[qubit]<<std::endl;
    targets.push_back(var_indices[qubit]);  
  }   
}
//Real type
Gate::Gate(gate_t type, std::vector<line_t>& controls, line_t target){
  gate_id = ++gate_count;
  for (auto c : controls)
    this->controls.push_back(c);
  targets.push_back(target);
  this->type = type;
}
void Gate::display(ostream & os){
  //cout<<"Gate Type: ";
  switch(type){
    case X:
      os<<"x ";
      break;
    case CX:
      os<<"cx ";
      break;
    case CCX:
      os<<"ccx ";
      break;
    case MCX:
      os<<"mcx ";
      break;
    case Y:
      os<<"y ";
      break;
    case Z:
      os<<"z ";
      break;
    case H:
      os<<"h ";
      break;
    case S:
      os<<"s ";
      break;
    case SDG:
      os<<"sdg ";
      break;
    case T:
      os<<"t ";
      break;
    case TDG:
      os<<"tdg ";
      break;
    case RZ:{
      std::ostringstream strs;
      strs << lamda;
      std::string str = strs.str();
      os<<"rz("<<removeZero(str)<<") ";
      //cout<<angle3<<" "<<str<<" "<<removeZero(str)<<endl;
      break;
    }
    case U1:{
      std::ostringstream strs;
      strs << lamda;
      std::string str = strs.str();
      os<<"u1("<<removeZero(str)<<") ";
      //cout<<angle3<<" "<<str<<" "<<removeZero(str)<<endl;
      break;
    }
    case U2:{
      std::ostringstream strs2;
      strs2 << phi;
      std::string str2 = strs2.str();
      std::ostringstream strs3;
      strs3 << lamda;
      std::string str3 = strs3.str();
      os<<"u2("<<removeZero(str2)<<","<<removeZero(str3)<<") ";
      //cout<<angle3<<" "<<str<<" "<<removeZero(str)<<endl;
      break;
    }
    case U3:{
      std::ostringstream strs1;
      strs1 << theta;
      std::string str1 = strs1.str();
      std::ostringstream strs2;
      strs2 << phi;
      std::string str2 = strs2.str();
      std::ostringstream strs3;
      strs3 << lamda;
      std::string str3 = strs3.str();
      os<<"u3("<<removeZero(str1)<<","<<removeZero(str2)<<","<<removeZero(str3)<<") ";
      //cout<<angle3<<" "<<str<<" "<<removeZero(str)<<endl;
      break;
    }
    case SWAP:
      cout<<"swap ";
      break;
    default:
      os<<"Invalid choice"<<endl;
  }
  //cout<<"Controls:";
  for(vector<line_t>::iterator i=controls.begin(); i!=controls.end(); ++i)
    os<<"q["<<*i<<"],";
  //cout<<endl;
  //cout<<"Targets:";
  for(vector<line_t>::iterator i=targets.begin(); i!=targets.end(); ++i)
    os<<"q["<<*i<<"]"<< (targets.size() > 1 && (i + 1) != targets.end() ? "," : ";");
  os<<endl;
}  
void Gate::display2(ostream & os){
  //cout<<"Gate Type: ";
  switch(type){
    case X:
      os<<"qc.x(";
      break;
    case CX:
      os<<"qc.cx(";
      break;
    case CCX:
      os<<"qc.ccx(";
      break;
    case MCX:
      os<<"qc.mct([";
      break;
    case Y:
      os<<"qc.y(";
      break;
    case Z:
      os<<"qc.z(";
      break;
    case H:
      os<<"qc.h(";
      break;
    case S:
      os<<"qc.s(";
      break;
    case SDG:
      os<<"qc.sdg(";
      break;
    case T:
      os<<"qc.t(";
      break;
    case TDG:
      os<<"qc.tdg(";
      break;
    case RZ:{
      std::ostringstream strs;
      strs << lamda;
      std::string str = strs.str();
      os<<"qc.u1("<<removeZero(str)<<", ";
      //cout<<angle3<<" "<<str<<" "<<removeZero(str)<<endl;
      break;
    }
    case U1:{
      std::ostringstream strs;
      strs << lamda;
      std::string str = strs.str();
      os<<"qc.u1("<<removeZero(str)<<", ";
      //cout<<angle3<<" "<<str<<" "<<removeZero(str)<<endl;
      break;
    }
    case U2:{
      std::ostringstream strs2;
      strs2 << phi;
      std::string str2 = strs2.str();
      std::ostringstream strs3;
      strs3 << lamda;
      std::string str3 = strs3.str();
      os<<"qc.u2("<<removeZero(str2)<<","<<removeZero(str3)<<", ";
      //cout<<angle3<<" "<<str<<" "<<removeZero(str)<<endl;
      break;
    }
    case U3:{
      std::ostringstream strs1;
      strs1 << theta;
      std::string str1 = strs1.str();
      std::ostringstream strs2;
      strs2 << phi;
      std::string str2 = strs2.str();
      std::ostringstream strs3;
      strs3 << lamda;
      std::string str3 = strs3.str();
      os<<"qc.u3("<<removeZero(str1)<<","<<removeZero(str2)<<","<<removeZero(str3)<<", ";
      //cout<<angle3<<" "<<str<<" "<<removeZero(str)<<endl;
      break;
    }
    case SWAP:
      cout<<"qc.swap(";
      break;
    default:
      os<<"Invalid choice"<<endl;
  }
  //cout<<"Controls:";  
  if (controls.size() > 0) {
    for(vector<line_t>::iterator i=controls.begin(); i!=controls.end() - 1; ++i)
      os<<*i<<",";//os<<"q["<<*i<<"],";
    if (type == MCX) os<<controls[controls.size() - 1]<<"],";//os<<"q["<<controls[controls.size() - 1]<<"]],";
    else os<<controls[controls.size() - 1]<<",";
  }
  //cout<<endl;
  //cout<<"Targets:";
  for(vector<line_t>::iterator i=targets.begin(); i!=targets.end(); ++i)
    os<<*i<< (targets.size() > 1 && (i + 1) != targets.end() ? ", " : ")");
  //os<<"q["<<*i<<"]"<< (targets.size() > 1 && (i + 1) != targets.end() ? ", " : ");");
  os<<endl;
}
