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
#include "unitSystem.hpp"
#include "parameters/parameters.hpp"
#include "container/universe.hpp"
#include "engine/engineGMX.hpp"
#include "parser/energyParserGMX.hpp"

//
// SimulatorBase class
// 
// a base class that implements the interface
// for all methods that perform simulations
//

class SimulatorBase
{
  protected:
    Universe                          universe {};
    std::unique_ptr<EngineBase>       mdEngine      {nullptr};
    std::unique_ptr<EnergyParserBase> energyParser  {nullptr};
    
    std::size_t currentCycle {1};
    std::size_t lastReactiveCycle {0};
    std::size_t nCycles {0};
    std::size_t nCyclesCompleted {0};

    bool          writeStatistics {false};
    std::ofstream STATISTICS_FILE {};

    std::unique_ptr<UnitSystem>  unitSystem {nullptr}; 

    // some generally usable functions:
    void mdSequence();

    // some functions that need to be implemented in derived:
    virtual void reactiveStep() = 0;
    virtual bool acceptance(const ReactionCandidate&) = 0;

    // make constructor protected to make the class purely virtual
    SimulatorBase() = default;

  public:
    virtual ~SimulatorBase() = default;

    // some generally usable functions:
    void run();
    void writeRestartFile(const Parameters&) const;

    // some functions that need to be implemented in derived:
    virtual void setup(const Parameters&);
    virtual void finish() = 0;

    auto getNCycles() const { return nCyclesCompleted; };

};

// need to include controller.hpp here!
// else Controller::signal() and Controller::SIGNAL cannot be used within SimulatorBase class
#include "control/controller.hpp"