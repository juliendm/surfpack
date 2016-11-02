/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "surfpack.h"
#include "AxesBounds.h"
#include "SurfpackParserArgs.h"
#include "SurfpackParser.h"
#include "SurfpackInterface.h"
#include "SurfpackInterpreter.h"
#include "SurfData.h"
#include "ModelFactory.h"
#include "SurfpackModel.h"

using std::cerr;
using std::cout;
using std::endl;
using std::ostream;
using std::runtime_error;
using std::string;
using std::vector;
using std::ostringstream;
///////////////////////////////////////////////////////////////////////////////
/////		SurfpackInterpreter namespace functions			  /////
///////////////////////////////////////////////////////////////////////////////

SurfpackInterpreter::SurfpackInterpreter() : parser(SurfpackParser::instance())
{

}

SurfpackInterpreter::~SurfpackInterpreter()
{

}

void SurfpackInterpreter::execute(const string* input_string, 
  const string* output_string)
{
  if (parser.yyparse(input_string, output_string) == 0) {
    commandLoop(cout , cerr);
  } else {
    cerr << "Command(s) not executed." << endl;
    cerr << SurfpackParser::cmdstream.str() << endl; 
  }
}

void SurfpackInterpreter::commandLoop(ostream& os, ostream& es)
{
  const vector<ParsedCommand>& fullCommands = parser.commandList();
  vector<Command>& commands = parser.comms;
  for (unsigned i = 0; i < commands.size(); i++) {
    try {
      //if (commands[i].isShellCommand()) {
      //  executeShellCommand(commands[i]);
      if (commands[i].first == "CreateSample") {
        execCreateSample(commands[i].second);
      } else if (commands[i].first== "CreateSurface") {
        execCreateSurface(commands[i].second);
      } else if (commands[i].first == "Evaluate") {
        execEvaluate(commands[i].second);
      } else if (commands[i].first == "Fitness") {
        execFitness(commands[i].second);
      } else if (commands[i].first == "Load") {
        execLoad(commands[i].second);
      } else if (commands[i].first == "Save") {
        execSave(commands[i].second);
      } else if (commands[i].first == "CreateAxes") {
        execCreateAxes(commands[i].second);
      } else {
        es << "Unrecognized command: " << commands[i].first << endl;
      }
    } catch (string& msg) {
      es << "Error: " << msg;
      es << "\n  Failed Instruction: " << fullCommands[i].cmdstring << endl;
    } catch (std::exception& e) {
      es << "Exception: " << e.what() 
	 << "\n  Failed Instruction: " << fullCommands[i].cmdstring << endl;
    } catch (...) {
      es << "Exception (unknown)\n  FailedInstruction: " 
	 << fullCommands[i].cmdstring << endl;
    }
  }
}
  
int getResponseIndex(const ArgList& arglist, const SurfData& sd)
{
  bool valid;
  bool is_response;
  unsigned int response_index;
  string response_name = SurfpackParser::parseIdentifier("response",
	arglist,false);
  if (response_name == "") {
      response_index = 
      SurfpackParser::parseInteger("response_index",arglist,valid,false);
      if (!valid) {
        // Neither response nor response_index specified
        // Default index is 0
        return 0;
      } else {
        return response_index;
      }
  } else {
    ostringstream os;
    valid = sd.varIndex(response_name, response_index, is_response);
    if (!valid) {
      os << "No response named '" << response_name << "' found." << endl;
      throw(os.str());
    } else if (!is_response) {
      os << "'" << response_name << "' is a predictor variable, but a"
	   << " response variable was requested" << endl;
      throw(os.str());
    } else {
      return response_index;
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
/////		SurfpackInterpreter namespace functions			  /////
///////////////////////////////////////////////////////////////////////////////

void SurfpackInterpreter::execLoad(ParamMap& args)
{
  bool valid;
  string filename = asStr(args["file"]);
  if (surfpack::hasExtension(filename,".sps") || 
      surfpack::hasExtension(filename,".bsps")) {
    execLoadSurface(args);
  } else if (surfpack::hasExtension(filename,".spd") ||
	     surfpack::hasExtension(filename,".bspd") ||
	     surfpack::hasExtension(filename,".dat")) {
    execLoadData(args);
  } else {
    throw string("Expected file extension: .sps/.bsps (surface) or " 
		 ".spd/.bspd/.dat (data)");
  }
}

void SurfpackInterpreter::execLoadData(ParamMap& args)
{
  try {
  SurfData* data = 0;
  string name = asStr(args["name"]); 
  string filename = asStr(args["file"]); 
  bool valid = false;
  int n_vars = asInt(args["n_predictors"],valid);
  if (valid) {
    // n_responses is required if n_vars is present
    int n_responses = asInt(args["n_responses"]); // fail on error
    int n_cols_to_skip = asInt(args["n_cols_to_skip"],valid);
    if (!valid) n_cols_to_skip = 0;
    data = SurfpackInterface::LoadData(filename,n_vars,n_responses,n_cols_to_skip);
  } else {
    data = SurfpackInterface::LoadData(filename);
  }
  assert(data);
  symbolTable.dataVars.insert(SurfDataSymbol(name,data));
  } catch (...) {
    cout << "Error caught in execLoadData: Did you forget to specify n_predictors and n_responses?" << endl;
  }
}

///\todo Add support for LoadSurface in interpreter
void SurfpackInterpreter::execLoadSurface(ParamMap& args)
{
  string name = asStr(args["name"]); 
  string filename = asStr(args["file"]); 
  SurfpackModel* model = SurfpackInterface::LoadModel(filename);
  symbolTable.modelVars.insert(SurfpackModelSymbol(name, model));
}

void SurfpackInterpreter::execSaveData(ParamMap& args)
{
  string data_name = asStr(args["data"]);
  string filename = asStr(args["file"]);
  SurfData* sd = symbolTable.lookupData(data_name);
  // Call SaveData 
  SurfpackInterface::Save(sd, filename);
}


void SurfpackInterpreter::execSaveSurface(const SurfpackModel* model, 
					  const string& filename)
{
  SurfpackInterface::Save(model, filename);
}


void SurfpackInterpreter::execSave(ParamMap& args)
{
  try {
    string filename = asStr(args["file"]); 
    // Don't automatically fail if either of these isn't defined
    bool valid_data;
    string data_name = asStr(args["data"],valid_data); 
    bool valid_model;
    string model_name = asStr(args["surface"],valid_model);
    if (!valid_data) {
      if (!valid_model) {
        // Do fail if both are missing
        throw string("Save command requires either 'surface' or 'data' argument");
      } else {
        SurfpackModel* model = symbolTable.lookupModel(model_name);
	execSaveSurface(model, filename);
      }
    } else {
      if (valid_model) {
        // Fail if both are specified 
        throw string("Save command may not have both 'surface' and 'data' arguments");
      } else {
        SurfData* sd = symbolTable.lookupData(data_name);
        sd->write(filename);
      }
    }
  } catch (string& e) {
    cout << "Error (Save): " << e << endl;
  } catch (std::exception& e) {
    cout << "Exception (Save): " << e.what() << endl;
  } catch (...) {
    cout << "Exception (Save, unknown)" << endl;
  }
}

//int getResponseIndex(const ArgList& arglist, const SurfData& sd)
//{
//  bool valid;
//  bool is_response;
//  unsigned int response_index;
//  string response_name = SurfpackParser::parseIdentifier("response",
//	arglist,false);
//  if (response_name == "") {
//      response_index = 
//      SurfpackParser::parseInteger("response_index",arglist,valid,false);
//      if (!valid) {
//        // Neither response nor response_index specified
//        // Default index is 0
//        return 0;
//      } else {
//        return response_index;
//      }
//  } else {
//    ostringstream os;
//    valid = sd.varIndex(response_name, response_index, is_response);
//    if (!valid) {
//      os << "No response named '" << response_name << "' found." << endl;
//      throw(os.str());
//    } else if (!is_response) {
//      os << "'" << response_name << "' is a predictor variable, but a"
//	   << " response variable was requested" << endl;
//      throw(os.str());
//    } else {
//      return response_index;
//    }
//  }
//}
//
void SurfpackInterpreter::execCreateSurface(ParamMap& args)
{
  // Extract the variable name for this SurfData object
  string name = asStr(args["name"]);
  string data = asStr(args["data"]);
  bool valid = false;
  SurfData* sd = symbolTable.lookupData(data);
  // Call CreateSurface
  SurfpackModelFactory* smf = ModelFactory::createModelFactory(args);
  SurfpackModel* model = smf->Build(*sd);
  delete smf;
  assert(model);
  symbolTable.modelVars.insert(SurfpackModelSymbol(name,model));
}

void SurfpackInterpreter::execCreateSample(ParamMap& args)
{
  // Extract the variable name for this SurfData object
  string name = asStr(args["name"]);
  string axes = asStr(args["axes"]);
  bool valid_grid_points;
  VecUns grid_points = asVecUns(args["grid_points"], valid_grid_points);
  bool valid_size = false;
  int n_points = asInt(args["size"],valid_size); 
  AxesBounds* ab = symbolTable.lookupAxes(axes);
  SurfData* sd = 0;
  if (valid_size) {
    if (valid_grid_points) { // both size and grid_points specified
      throw string("Cannot specify both size and grid_points");
    } else { // only size specified
	// MonteCarlo sample
      sd = SurfpackInterface::CreateSample(ab,n_points); 
    }
  } else {
    if (!grid_points.empty()) { // only grid_points specified
	// grid sample
      sd = SurfpackInterface::CreateSample(ab,grid_points); 
    } else { // neither specified
      throw string("Must specify either size or grid_points");
    }
  }
  bool valid_functions;
  VecStr test_functions = asVecStr(args["test_functions"], valid_functions);
  if (valid_functions) {
    SurfpackInterface::Evaluate(sd,test_functions);
  }
  symbolTable.dataVars.insert(SurfDataSymbol(name,sd));
}

void SurfpackInterpreter::execEvaluate(ParamMap& args)
{
  // Extract the variable name for this SurfData object
  string surf_name = asStr(args["surface"]);
  string data = asStr(args["data"]);
  SurfpackModel* model = symbolTable.lookupModel(surf_name);
  SurfData* sd = symbolTable.lookupData(data);
  // Call Evaluate
  VecDbl results = (*model)(*sd);
  sd->addResponse(results,args["label"]); 
}

void SurfpackInterpreter::execFitness(ParamMap& args)
{
  string surface = asStr(args["surface"]);
  bool valid_data;
  string data = asStr(args["data"],valid_data);
  SurfpackModel* model = symbolTable.lookupModel(surface);
  SurfData* sd = 0;
  if (valid_data) {
    sd = symbolTable.lookupData(data);
  }
  string metric = asStr(args["metric"]);
  // Extract the response index for the input data set
  bool valid_response_index;
  int response_index = asInt(args["response_index"],valid_response_index);
  bool valid_n;
  int n = asInt(args["n"],valid_n);
  if (!valid_response_index) { // No response_index was specified, use 0
    response_index = 0;
  }
  double fitness;
  if (valid_data) {
    fitness = SurfpackInterface::Fitness(model,sd,metric,response_index,n);
  } else {
    fitness = SurfpackInterface::Fitness(model,metric,response_index,n); 
  }
  cout << metric << " for " << surface;
  if (data != "") cout << " on " << data;
  cout << ": " << fitness << endl;
}
  
void SurfpackInterpreter::execCreateAxes(ParamMap& args)
{
  string name = asStr(args["name"]);
  string bounds = asStr(args["bounds"]);
  AxesBounds* ab = SurfpackInterface::CreateAxes(bounds);
  symbolTable.axesVars.insert(AxesBoundsSymbol(name,ab));
}
//
//void SurfpackInterpreter::
//  executeShellCommand(ParamMap& args)
//{
//  system(c.cmdstring.c_str());
//}

int SurfpackInterpreter::asInt(const string& arg)
{
  if (arg == "") throw string("Expected integer value");
  return atoi(arg.c_str());
}

std::string SurfpackInterpreter::asStr(const string& arg)
{
  if (arg == "") throw string("Expected string value");
  return arg;
}

double SurfpackInterpreter::asDbl(const string& arg)
{
  if (arg == "") throw string("Expected double value");
  return atof(arg.c_str());
}

VecUns SurfpackInterpreter::asVecUns(const string& arg)
{
  if (arg == "") throw string("Expected vector unsigned");
  return surfpack::toVec<unsigned>(arg); 
}

VecStr SurfpackInterpreter::asVecStr(const string& arg)
{
  if (arg == "") throw string("Expected vector string");
  return surfpack::toVec<std::string>(arg); 
}

int SurfpackInterpreter::asInt(const string& arg, bool& valid)
{
  if (arg == "") {
    valid = false;
    return 0;
  }
  valid = true;
  return atoi(arg.c_str());
}

std::string SurfpackInterpreter::asStr(const string& arg, bool& valid)
{
  if (arg == "") {
    valid = false;
    return "";
  }
  valid = true;
  return arg; 
}

double SurfpackInterpreter::asDbl(const string& arg, bool& valid)
{
  if (arg == "") {
    valid = false;
    return 0;
  }
  valid = true;
  return atof(arg.c_str()); 
}

VecUns SurfpackInterpreter::asVecUns(const string& arg, bool& valid)
{
  if (arg == "") {
    valid = false;
    return VecUns();
  }
  valid = true;
  return surfpack::toVec<unsigned>(arg); 
}

VecStr SurfpackInterpreter::asVecStr(const string& arg, bool& valid)
{
  if (arg == "") {
    valid = false;
    return VecStr();
  }
  valid = true;
  return surfpack::toVec<std::string>(arg); 
}
//********************************************************************
SurfpackInterpreter::SymbolTable::~SymbolTable() 
{ 
  for (SurfDataMap::iterator iter = dataVars.begin();
       iter != dataVars.end();
       ++iter) {
    delete iter->second; 
  }
  for (SurfpackModelMap::iterator iter = modelVars.begin();
       iter != modelVars.end();
       ++iter) {
    delete iter->second; 
  }
  for (AxesBoundsMap::iterator iter = axesVars.begin();
       iter != axesVars.end();
       ++iter) {
    delete iter->second; 
  }
}

SurfpackModel* 
SurfpackInterpreter::SymbolTable::lookupModel(const std::string name)
{
  SurfpackModel* result = modelVars[name];
  if (result == 0) {
    cout << "Bad lookup; table size:  " << modelVars.size() << endl;
    for (SurfpackModelMap::iterator itr = modelVars.begin();
        itr != modelVars.end(); ++itr) {
      cout << itr->first << " " << itr->second << endl;
    }
    
    string msg = "Model variable " + name + " not found in symbol table."; 
    throw msg;
  }
  return result;
}

SurfData* SurfpackInterpreter::SymbolTable::lookupData(string name)
{
  SurfDataMap::iterator iter = dataVars.find(name);
  if (iter == dataVars.end()) {
    string msg = "Data variable " + name + " not found in symbol table."; 
    throw msg;
  }
  assert(iter->second);
  return iter->second;
}

AxesBounds* SurfpackInterpreter::SymbolTable::lookupAxes(string name)
{
  AxesBoundsMap::iterator iter = axesVars.find(name);
  if (iter == axesVars.end()) {
    string msg = "Axes variable " + name + " not found in symbol table."; 
    throw msg;
  }
  assert(iter->second);
  return iter->second;
}
