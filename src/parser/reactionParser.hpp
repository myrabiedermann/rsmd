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

#include "reaction/reactionBase.hpp"
#include "enhance/utility.hpp"

#include <fstream>

//
// a reaction input file reader (/ writer)
//

class ReactionParser
{

  public:
    ReactionBase read(const std::string&);
    std::string  writeExample();

};