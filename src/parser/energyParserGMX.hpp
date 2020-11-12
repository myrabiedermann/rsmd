#pragma once


#include "parser/energyParserBase.hpp"

#include <sstream>
#include <fstream>

//
// energy parser class
// reads energy differences from 
// gromacs (GMX) files
//


class EnergyParserGMX: public EnergyParserBase
{
  private:
    bool computeLocalPotentialEnergy {false};
    bool computeSolvationPotentialEnergy {false};
    REAL potentialEnergyAverageTime {0.0};
    REAL readPotentialEnergy( const std::string& );
    REAL readSolvationEnergy( const std::string& );


  public:
    ~EnergyParserGMX() = default;
    EnergyParserGMX()  = default;

    REAL readPotentialEnergyDifference( const std::size_t&, const std::size_t& );        

    void setup(const Parameters&);
};