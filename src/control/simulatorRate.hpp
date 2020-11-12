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

#include "control/simulatorBase.hpp"
#include <unordered_map>

//
// SimulatorRate class
// 
// inherits interface from SimulatorBase;
// implements reactiveStep() for a 
// hybrid MC/MD simulation with rate based acceptance criterion
//

class SimulatorRate : public SimulatorBase
{
  private:
    std::size_t nCyclesReaction {0};
    std::size_t nCyclesNoReaction {0};
    std::size_t nCyclesFailedFirstRelaxation {0};

    REAL rsFrequency {0};

    // some functions that need to be implemented in derived:
    void reactiveStep();
    bool acceptance(const ReactionCandidate&);

  public:
    SimulatorRate() = default;
    ~SimulatorRate() = default;

    // some functions that need to be implemented in derived:
    void finish();
    void setup(const Parameters&);

};