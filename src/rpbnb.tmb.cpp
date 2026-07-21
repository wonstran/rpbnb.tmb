// src/rpbnb.tmb.cpp — will be filled in Tasks 3–4
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () {
  DATA_SCALAR(placeholder);
  PARAMETER(dummy);
  Type nll = dummy * dummy;
  return nll;
}
