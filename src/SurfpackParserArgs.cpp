/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

#include "SurfpackParserArgs.h"

Triplet::Triplet() : min(0), max(0), numPts(0) {}

const std::string& Rval::argType() const
{
  const static std::string s("none");
  return s;
}

int Rval::getInteger() const
{
  static int dummy = 0;
  noSuchValue();
  return dummy;
}

double Rval::getReal() const
{
  static double dummy = 0.0;
  noSuchValue();
  return dummy;
}

const Tuple& Rval::getTuple() const
{
  static Tuple t;
  noSuchValue();
  return t;
}

const Triplet& Rval::getTriplet() const
{
  static Triplet t;
  noSuchValue();
  return t;
}

const std::string& Rval::getIdentifier() const
{
  static std::string s;
  noSuchValue();
  return s;
}

const std::string& Rval::getStringLiteral() const
{
  static std::string s;
  noSuchValue();
  return s;
}

const ArgList& Rval::getArgList() const
{
  static ArgList al;
  noSuchValue();
  return al;
}

void Rval::noSuchValue() const
{
  throw std::string("This Rval class does not have such a value");
}

Rval::~Rval()
{

}

RvalInteger::RvalInteger(int value_in) : value(value_in) 
{

}

int RvalInteger::getInteger() const
{
  return value;
}

Rval* RvalInteger::clone() const
{
  return new RvalInteger(value);
}

RvalReal::RvalReal(double value_in) : value(value_in) 
{

}

double RvalReal::getReal() const
{
  return value;
}

Rval* RvalReal::clone() const
{
  return new RvalReal(value);
}

RvalTuple::RvalTuple(const Tuple& value_in) : value(value_in) 
{

}

RvalTuple::RvalTuple(const std::vector<double>& value_in)
{
  value.resize(value_in.size());
  for (unsigned i = 0; i < value_in.size(); i++) {
    std::ostringstream os;
    os << value_in[i];
    value[i] = os.str();
  }
}

const Tuple& RvalTuple::getTuple() const
{
  return value;
}

Rval* RvalTuple::clone() const
{
  return new RvalTuple(value);
}

const std::string& RvalTuple::argType() const
{
  const static std::string s("tuple");
  return s;
}

const std::vector< double >& 
RvalTuple::asVectorDouble(std::vector< double >& result, const Tuple& tuple)
{
  result.resize(tuple.size());
  for (unsigned i = 0; i < tuple.size(); i++) {
    result[i] = std::atof(tuple[i].c_str());
  }
  return result;
}

const std::vector< std::string >& 
RvalTuple::asVectorString(std::vector< std::string >& result, 
  const Tuple& tuple)
{
  result.resize(tuple.size());
  for (unsigned i = 0; i < tuple.size(); i++) {
    result[i] = tuple[i];
  }
  return result;
}

RvalTriplet::RvalTriplet(const Triplet& value_in) : value(value_in) 
{

}

const Triplet& RvalTriplet::getTriplet() const
{
  return value;
}

Rval* RvalTriplet::clone() const
{
  return new RvalTriplet(value);
}

RvalIdentifier::RvalIdentifier(const std::string& value_in) : value(value_in) 
{

}

const std::string& RvalIdentifier::getIdentifier() const
{
  return value;
}

Rval* RvalIdentifier::clone() const
{
  return new RvalIdentifier(value);
}

RvalStringLiteral::RvalStringLiteral(const std::string& value_in) 
  : value(value_in) 
{

}

const std::string& RvalStringLiteral::getStringLiteral() const
{
  return value;
}

Rval* RvalStringLiteral::clone() const
{
  return new RvalStringLiteral(value);
}

RvalArgList::RvalArgList(const ArgList& value_in) : value(value_in)
{

}

const ArgList& RvalArgList::getArgList() const
{
  return value;
}

Rval* RvalArgList::clone() const
{
  return new RvalArgList(value);
}

const std::string& RvalArgList::argType() const
{
  const static std::string s("arglist");
  return s;
}

Arg::Arg(const Arg& other)
  : name(other.name), rval(0)
{
  if (other.rval) {
    this->rval = other.rval->clone();
  }  
}

const Arg& Arg::operator=(const Arg& other)
{
  this->name = other.name;
  delete this->rval;
  if (other.rval) {
    this->rval = other.rval->clone();
  } else {
    this->rval = 0;
  }
  return *this;
} 

Arg::Arg(const std::string& name_in, Rval* rval_in)
  : name(name_in), rval(rval_in)
{
}

Arg::Arg()
 : rval(0)
{

}

void Arg::setName(const std::string& name_in)
{
  name = name_in;
}

const Rval* Arg::getRVal() const
{
  return rval;
}

void Arg::setRVal(Rval* rval_in)
{
  delete rval;
  rval = rval_in;
}

Arg::~Arg()
{
  delete rval;
  rval = 0;
}

Arg Arg::makeArg(const std::string name, int rval)
{
  return Arg(name, new RvalInteger(rval));
}
