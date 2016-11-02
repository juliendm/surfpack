/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#ifndef MY_PARSING_STRUCTURES_H
#define MY_PARSING_STRUCTURES_H

#include "surfpack_system_headers.h"
#include "SurfpackParserArgs.h"

class FlexWrapper;

class ParsedCommand
{ 
public:
  ParsedCommand();
  ParsedCommand(bool shell_command);
  bool isShellCommand() const;
  bool shellCommand;
  std::string name;
  ArgList arglist;
  std::string cmdstring;
};

/// Singleton class.  Works in conjunction with flex and bison
/// to parse out Surfpack commands from an input file
class SurfpackParser 
{
public:
  // commands for use by the Bison parser
  void addCommandName();
  void printComms();
  void appendArg();
  void addArg();
  void addArgName();
  void addArgValIdent();
  void addArgValInt();
  void addArgValString();
  void addArgValReal();
  void addArgValTuple();
  void addArgValArgList();
  void addArgListToCommand();
  void addTupleVal();
  void newTuple();
  void pushNewArgList();
  void popArgList();
  void init();
  void storeCommandString();
  void shellCommand();
  /// returns the only instance of this class
  static SurfpackParser& instance();
  static std::ostringstream cmdstream;
  FlexWrapper& globalLexer();
  int yyparse(const std::string* input_string,const std::string* output_string);
  std::vector<ParsedCommand>& commandList();

  // commands for use by interpreter to parse out certain argument types
  static std::string parseIdentifier(const std::string& argname,
    const ArgList& arglist, bool throwExIfAbsent = true);
  static std::string parseStringLiteral(const std::string& argname,
    const ArgList& arglist, bool throwExIfAbsent = true);
  static int parseInteger(const std::string& argname,
    const ArgList& arglist, bool& valid, bool throwExIfAbsent = true);
  static std::vector<double> parseTuple(const std::string& argname,
    const ArgList& arglist, bool throwExIfAbsent = true);
  static std::vector<unsigned> parseUnsignedTuple(const std::string& argname,
    const ArgList& arglist, bool throwExIfAbsent = true);
  static std::vector<std::string> parseStringTuple(const std::string& argname,
    const ArgList& arglist, bool throwExIfAbsent = true);
  static std::vector<std::string> parseMultiString(const std::string& argname,
    const ArgList& arglist, bool throwExIfAbsent = true);
protected:
  
  class missing_argument 
  {
  public:
    missing_argument(const std::string& name_in = "", 
      const std::string& cmdstring_in = "") 
      : name(name_in), cmdstring(cmdstring_in) {}
    void print() { std::cerr << "Error in " << cmdstring << ":  " 
			     << ".  Required parameter \"" << name
			     << " is not specified." << std::endl; }
  protected:
    std::string name;
    std::string cmdstring;
  };
// Default constructor, copy constructor and assignment operator declared 
// protected and not implemented, in order to make sure that only one
// instance is created
  SurfpackParser();
  ~SurfpackParser();
  SurfpackParser(const SurfpackParser&);
  const SurfpackParser& operator=(const SurfpackParser&);

  // For keeping track of parsed commands
public:
  static std::string argname;
  static std::string argval;
  static ParamMap params;
  static std::vector<Command> comms;

private:
  std::vector<ParsedCommand> commands;
  ArgList* currentArgList;
  int currentArgIndex;
  int currentTupleIndex;
  FlexWrapper* global_lexer;
  Tuple* currentTuple;
  std::stack< ArgList > arglistStack;
};
#endif
