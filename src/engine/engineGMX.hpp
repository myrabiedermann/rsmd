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

#include "engine/engineBase.hpp"
#include "enhance/utility.hpp"

#include <thread>
#include <filesystem>

//
// a derived class that implements the gromacs engine
// derived from engineBase
//

class EngineGMX : public EngineBase
{
  private:
    // all sorts of necessary stuff like filenames, n of threads to use etc.
    std::string executablePath {};

    std::string mdp_file {};
    std::string mdp_file_relaxation {};
    std::string mdp_file_energy {};

    std::string nt_as_str {};
    std::string ntmpi_as_str {};
    std::string ntomp_as_str {};

    REAL  extensionTime {1};
    std::string  extensionTime_str {"1"};

    bool computeLocalPotentialEnergies {false};
    bool computeSolvationPotentialEnergies {false};
    bool averagePotentialEnergies {false};

    bool        saveRejectedFiles {false};
    std::vector<std::string>  rejectedFilekeys {};
    std::string backupPolicy {"-nobackup"};

    // helper functions
    void grompp( const std::string&, const std::string&, const std::string&, const std::string& );
    void grompp( const std::string&, const std::string&, const std::string&, const std::string&, const std::string& );
    void convert_tpr( const std::string&, const std::string&, const std::string& );
    void convert_tpr( const std::string&, const std::string&);
    void trjconv( const std::string&, const std::string&, const std::string&, const std::string& );
    void mdrun( const std::string& );
    void mdrun( const std::string&, const std::string&, const std::string& );
    void mdrunRerun( const std::string&, const std::string&, const std::string& );
    void energy( const std::string&, const std::string& );
    void energySolvation( const std::string&, const std::string& );
    void read_mdp( const std::string& );


  public:
    EngineGMX() = default;
    ~EngineGMX() = default;

    void setup(const Parameters&);
    void verifyExecutable();
    void runMD( const std::size_t& );
    void runMDInitial();
    void runMDAppending( const std::size_t&, const std::size_t& );
    bool runRelaxation( const std::size_t& );
    void runEnergyComputation( const std::size_t&, const std::size_t& );
    void cleanup( const std::size_t& );
};