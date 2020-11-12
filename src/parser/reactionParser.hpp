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