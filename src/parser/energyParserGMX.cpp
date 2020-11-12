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

#include "parser/energyParserGMX.hpp"

//
// setup
//
void EnergyParserGMX::setup(const Parameters& parameters)
{
    potentialEnergyAverageTime = parameters.getOption("reaction.averagePotentialEnergy").as<REAL>();
    computeLocalPotentialEnergy = parameters.getOption("reaction.computeLocalPotentialEnergy").as<bool>();
    computeSolvationPotentialEnergy = parameters.getOption("reaction.computeSolvationPotentialEnergy").as<bool>();
}


//
// read potential energies from last (couple of) steps before/after reactive step
// and return difference
//
REAL EnergyParserGMX::readPotentialEnergyDifference( const std::size_t& cycle, const std::size_t& lastReactiveCycle )
{
    std::stringstream filenameBefore, filenameAfter {};
    filenameBefore << lastReactiveCycle << "-md.xvg";
    filenameAfter << cycle << "-rs.xvg";

    REAL energyDifference = 0;
    energyDifference = readPotentialEnergy(filenameAfter.str()) - readPotentialEnergy(filenameBefore.str());

    if( computeSolvationPotentialEnergy )
    {
        energyDifference += (readSolvationEnergy("products_solvation.xvg") - readSolvationEnergy("reactants_solvation.xvg"));   
    }

    return energyDifference;
}  


//
// read potential energies from .xvg file
// average them if requested, else read only energy from last step
//
// IDEA: read in reverse order through file in order to save time ?!
//       (--> benchmark before implementing!)
REAL EnergyParserGMX::readPotentialEnergy( const std::string& filename )
{
    REAL time = 0;
    REAL potentialEnergy = 0;
    REAL tmp = 0;
    REAL timeMargin = 0;
    std::string line {};
    std::stringstream linestream {};
    std::size_t counter = 0;

    std::ifstream FILE( filename );
    if( ! FILE )
    {
        rsmdCRITICAL( "could not read file '" << filename << "', cannot extract potential energy");
    }

    // read last line from file, save time & potential energy,
    FILE.seekg(-2, std::ios::end);
    while( FILE.peek() != '\n' )
    {
        FILE.seekg(-1, std::ios::cur);
    }
    FILE >> time >> potentialEnergy;
    rsmdDEBUG( "reading potential energy " << potentialEnergy );
    

    // if no averaging requested --> return potentialEnergy
    // else set timeMargin correctly and read + average remaining values
    if( potentialEnergyAverageTime != 0 )
    {
        FILE.seekg(0, std::ios::beg);

        // set timeMargin accordingly
        timeMargin = time - potentialEnergyAverageTime;
        rsmdDEBUG( "potentialEnergyAverageTime = " << potentialEnergyAverageTime << " ps");
        rsmdDEBUG( "reading potential energies in [" << timeMargin << ", " << time << "] (ps)" );
        if( timeMargin < 0 )
        {
            rsmdWARNING( "potentialEnergyAverageTime is larger than total relaxation sequence time (" << time << " < " << potentialEnergyAverageTime << ")" );
            rsmdWARNING( " setting potentialEnergyAverageTime to " << time << " ps.")
            timeMargin = 0;
        }

        // read file line by line and extract potential energies
        potentialEnergy = 0;
        while( std::getline(FILE, line, '\n') )
        {
            if( line[0] == '#' || line[0] == '@' ) continue;

            linestream = std::stringstream(line);
            linestream >> time >> tmp;
            
            if( time >= timeMargin )
            {
                potentialEnergy += tmp;
                ++ counter;
            }
        }
        potentialEnergy /= counter;
        rsmdDEBUG( "potentialEnergy = " << potentialEnergy << " kJ/mol (averaged over " << counter << " data points)" );
    }

    FILE.close();
    return potentialEnergy;
}


//
// read interaction energies with solvent from .xvg file
// average them if requested, else read only energy from last step
//
// IDEA: read in reverse order through file in order to save time ?!
//       (--> benchmark before implementing!)
REAL EnergyParserGMX::readSolvationEnergy( const std::string& filename )
{
    REAL time = 0;
    REAL energy_lj = 0;
    REAL energy_coulomb = 0;
    REAL tmp1 = 0;
    REAL tmp2 = 0;
    REAL timeMargin = 0;
    std::string line {};
    std::stringstream linestream {};
    std::size_t counter = 0;

    std::ifstream FILE( filename );
    if( ! FILE )
    {
        rsmdCRITICAL( "could not read file '" << filename << "', cannot extract potential energy");
    }

    // read last line from file, save time & potential energy,
    FILE.seekg(-2, std::ios::end);
    while( FILE.peek() != '\n' )
    {
        FILE.seekg(-1, std::ios::cur);
    }
    FILE >> time >> energy_coulomb >> energy_lj;
    rsmdDEBUG( "reading energies lj = " << energy_lj << ", coulomb = " << energy_coulomb );

    // if no averaging requested --> return potentialEnergy
    // else set timeMargin correctly and read + average remaining values
    if( potentialEnergyAverageTime != 0 )
    {
        FILE.seekg(0, std::ios::beg);

        // set timeMargin accordingly
        timeMargin = time - potentialEnergyAverageTime;
        rsmdDEBUG( "potentialEnergyAverageTime = " << potentialEnergyAverageTime << " ps");
        rsmdDEBUG( "reading potential energies in [" << timeMargin << ", " << time << "] (ps)" );
    
        // read file line by line and extract potential energies
        energy_lj = 0;
        energy_coulomb = 0;
        while( std::getline(FILE, line, '\n') )
        {
            if( line[0] == '#' || line[0] == '@' ) continue;

            linestream = std::stringstream(line);
            linestream >> time >> tmp1 >> tmp2;
            
            if( time >= timeMargin )
            {
                energy_coulomb += tmp1;
                energy_lj += tmp2;
                ++ counter;
                rsmdDEBUG(tmp1 << " " << tmp2);
            }
        }
        energy_lj /= counter;
        energy_coulomb /= counter;
        rsmdDEBUG( "lj energy = " << energy_lj << ", coulomb energy = " << energy_coulomb << " kJ/mol (averaged over " << counter << " data points)" );
    }

    FILE.close();
    return (energy_lj + energy_coulomb);
}
