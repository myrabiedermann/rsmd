/************************************************
 *                                              *
 *                rs@md                         *
 *    (reactive steps @ molecular dynamics )    *
 *                                              *
 ************************************************/
/* 
 Copyright 2020 Myra Biedermann
 Licensed under the Apache License, Version 2.0 
*/
#pragma once

#include "definitions.hpp"
#include "parameters/parameters.hpp"

//
// a base class that implements
// the interface for all energy readers
//

class EnergyParserBase
{
  protected:
    EnergyParserBase() = default;


  public:
    virtual ~EnergyParserBase() = default;

    virtual REAL readPotentialEnergyDifference( const std::size_t& , const std::size_t& ) = 0;

    virtual void setup(const Parameters&) = 0;

};