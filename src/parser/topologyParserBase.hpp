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
#include "container/topology.hpp"
#include "parameters/parameters.hpp"

//
// a base class that implements
// the interface for all topology readers/writers
//

class TopologyParserBase
{
  protected:
    TopologyParserBase() = default;

  public:
    virtual void read( Topology&, const std::size_t&) = 0;
    virtual void readRelaxed( Topology&, const std::size_t&) = 0;
    virtual void write(Topology&, const std::size_t&) = 0;

    virtual ~TopologyParserBase() = default;
};