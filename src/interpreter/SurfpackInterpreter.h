/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */


#include "surfpack_system_headers.h"
#include "AxesBounds.h"
#include "SurfpackParser.h"
class SurfData;
class ParsedCommand;
class SurfpackParser;
class SurfpackModel;

class SurfpackInterpreter
{
public:
  SurfpackInterpreter();
  ~SurfpackInterpreter();
  void execute(const std::string* input_string = 0, 
    const std::string* output_string = 0);
  void commandLoop(std::ostream& os = std::cout, std::ostream& es = std::cerr);
  void execCreateAxes(ParamMap& args);
  void execCreateSample(ParamMap& args);
  void execCreateSurface(ParamMap& args);
  void execEvaluate(ParamMap& args);
  void execFitness(ParamMap& args);
  void execLoad(ParamMap& args);
  void execLoadData(ParamMap& args);
  void execLoadSurface(ParamMap& args);
  void execSave(ParamMap& args);
  void execSaveData(ParamMap& args);
  void execSaveSurface(const SurfpackModel* model, const std::string& filename);
  void execShellCommand(ParamMap& args);

  ///\todo Collapse all of these conversion functions into one template function
  static int asInt(const std::string& arg);
  static std::string asStr(const std::string& arg);
  static double asDbl(const std::string& arg);
  static VecUns asVecUns(const std::string& arg);
  static VecStr asVecStr(const std::string& arg);
  static int asInt(const std::string& arg, bool& valid);
  static std::string asStr(const std::string& arg, bool& valid);
  static double asDbl(const std::string& arg, bool& valid);
  static VecUns asVecUns(const std::string& arg, bool& valid);
  static VecStr asVecStr(const std::string& arg, bool& valid);
protected:
  class command_error 
  {
  public:
    command_error(const std::string& msg_in = "", 
      const std::string& cmdstring_in = "") 
      : msg(msg_in), cmdstring(cmdstring_in) {}
    void print() { std::cerr << "Error in " << cmdstring << ":  " 
			     << msg << std::endl; }
  protected:
    std::string msg;
    std::string cmdstring;
  };
public:
  typedef std::pair<std::string, SurfData*> SurfDataSymbol;
  typedef std::map<std::string, SurfData*> SurfDataMap;
  typedef std::pair<std::string, AxesBounds*> AxesBoundsSymbol;
  typedef std::map<std::string, AxesBounds*> AxesBoundsMap;
  typedef std::pair<std::string, SurfpackModel*> SurfpackModelSymbol;
  typedef std::map<std::string, SurfpackModel*> SurfpackModelMap;
  
private:
  struct SymbolTable
  {
    SurfDataMap dataVars;
    SurfpackModelMap modelVars;
    AxesBoundsMap axesVars;
    ~SymbolTable();  
    SurfpackModel* lookupModel(const std::string);
    SurfData* lookupData(const std::string);
    AxesBounds* lookupAxes(const std::string);
  };

  struct SymbolTable symbolTable;
  SurfpackParser& parser;
};
    
