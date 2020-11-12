#pragma once


#include "control/simulatorBase.hpp"

//
// SimulatorMetropolis class
// 
// inherits interface from SimulatorBase;
// implements reactiveStep() for a 
// hybrid MC/MD simulation with Metropolis acceptance criterion
//


class SimulatorMetropolis : public SimulatorBase
{
  private:
    std::size_t nCyclesAccepted {0};
    std::size_t nCyclesRejected {0};
    std::size_t nCyclesRejectedFailedRelaxation {0};

    std::map<std::string, std::size_t> nCyclesFailedRelaxation_reactions {};
    REAL temperature {0};

    // some functions that need to be implemented in derived:
    void reactiveStep();
    bool acceptance(const ReactionCandidate&);

  public:
    SimulatorMetropolis() = default;
    ~SimulatorMetropolis() = default;

    // some functions that need to be implemented in derived:
    void finish();
    void setup(const Parameters&);

};