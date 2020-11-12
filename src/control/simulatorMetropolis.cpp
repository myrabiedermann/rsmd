#include "control/simulatorMetropolis.hpp"


//
// setup stuff specific to hybrid MC/MD simulation with Metropolis 
// acceptance criterion
//
void SimulatorMetropolis::setup(const Parameters& parameters)
{
    rsmdLOG( "setting up the simulation world ..." );

    // setup general stuff
    SimulatorBase::setup(parameters);

    // setup specific stuff
    temperature = parameters.getOption("reaction.temperature").as<REAL>();

    // setup map for counting failed relaxations:
    for( const auto& reaction: universe.getReactionTemplates() )
    {
        nCyclesFailedRelaxation_reactions[reaction.getName()] = 0;
    }

    // check statistics file and write header
    STATISTICS_FILE << "#" << std::setw(9) << "cycle"
                    << std::setw(15) << "# candidates";
    STATISTICS_FILE << std::setw(30) << "chosen_reaction"
                    << std::setw(10) << "acc/rej" << '\n';

    rsmdLOG( "... setup done, time to start the simulation!" );
    rsmdLOG( std::flush << std::setprecision(3) );
}


//
// do reactive step
//
void SimulatorMetropolis::reactiveStep()
{
    // search for candidates
    universe.update(lastReactiveCycle);
    auto candidates = universe.searchReactionCandidates(); // returns shuffled vector of reaction candidates
    STATISTICS_FILE << std::setw(10) << currentCycle << std::setw(15) << candidates.size();
    if( candidates.size() > 0 )
    {
        // count candidates per reaction type
        std::map<std::string, int> counts {};
        for( const auto& rType: universe.getReactionTemplates() )
        {
            auto count = std::accumulate(candidates.begin(), candidates.end(), 0, [&rType](int a, auto b){ return a + (b.getName() == rType.getName() ? 1 : 0); });
            counts[rType.getName()] =  count;
        }
        // compute weights
        std::vector<REAL> weights {}; 
        std::transform(candidates.begin(), candidates.end(), std::back_inserter(weights),
                    [&](const auto& c) -> REAL { return std::exp(-1.0 * c.getActivationEnergy() / (temperature*unitSystem->getR())); });
        // pick a candidate at random (but weighted) and perform reaction
        auto& candidate = *enhance::random_weighted_choice(candidates.begin(), weights.begin(), weights.end());
        rsmdLOG( "testing reaction candidate ");
        rsmdLOG( candidate.shortInfo() );
        STATISTICS_FILE << std::setw(30) << candidate.getName();
        universe.react(candidate);

        // relaxation
        universe.write(currentCycle);
        if( mdEngine->runRelaxation(currentCycle) )
        {
            // check acceptance / reverse if rejected
            mdEngine->runEnergyComputation(currentCycle, lastReactiveCycle);
            if( acceptance(candidate) )
            {
                lastReactiveCycle = currentCycle;
                ++ nCyclesAccepted;
                STATISTICS_FILE << std::setw(10) << "acc";
                // read configuration after relaxation and check if sensible
                universe.readRelaxed(currentCycle);
                universe.checkMovement(candidate);
            }
            else
            {
                // read configuration after relaxation and check if sensible
                universe.readRelaxed(currentCycle);
                universe.checkMovement(candidate);
                mdEngine->cleanup(currentCycle);
                ++ nCyclesRejected;
                STATISTICS_FILE << std::setw(10) << "rej";
            }
        }
        else
        {
            rsmdLOG( "... reactive step rejected! (due to a failed relaxation)" );
            mdEngine->cleanup(currentCycle);
            ++ nCyclesRejectedFailedRelaxation;
            ++ nCyclesFailedRelaxation_reactions[candidate.getName()];
            STATISTICS_FILE << std::setw(10) << "rej_relax";
        }
    }
    else
    {
        rsmdLOG( "... no reaction candidates available.")
        STATISTICS_FILE << std::setw(30) << "none" << std::setw(10) << "none" << std::setw(10) << "none";
    }
    STATISTICS_FILE << '\n' << std::flush;
}


//
// check acceptance
//
bool SimulatorMetropolis::acceptance(const ReactionCandidate& candidate)
{
    REAL random = enhance::random(0.0, 1.0);

    REAL energyDifference = energyParser->readPotentialEnergyDifference(currentCycle, lastReactiveCycle);
    rsmdLOG( "... potential energy difference = " << energyDifference << " + " << candidate.getReactionEnergy() 
                                         << " = " << energyDifference + candidate.getReactionEnergy() << ' ' << unitSystem->energy );
    energyDifference += candidate.getReactionEnergy();
    
    REAL condition = std::exp( -1.0 * energyDifference / (unitSystem->getR() * temperature) );

    if( random < condition )
    {
        rsmdLOG( "... candidate accepted: " << random << " < " << condition );
        return true;
    }
    else
    {
        rsmdLOG( "... candidate rejected: " << random << " !< " << condition );
        return false;
    }
}



//
// finish & clean up
//
void SimulatorMetropolis::finish() 
{
    STATISTICS_FILE.close();

    rsmdLOG( "" );
    rsmdLOG( "finished rs@md simulation" );
    rsmdLOG( "total " << (nCyclesAccepted + nCyclesRejected + nCyclesRejectedFailedRelaxation) << " cycles have been performed:" );
    rsmdLOG( "      " << nCyclesAccepted << " accepted" );
    rsmdLOG( "      " << nCyclesRejected << " rejected" );
    rsmdLOG( "      " << nCyclesRejectedFailedRelaxation << " rejected due to a failed relaxation" );
    rsmdLOG( "failed relaxations happened for: ");
    for( const auto& element: nCyclesFailedRelaxation_reactions )
    {
        rsmdLOG( "      " << element.second << " " << element.first );
    }
    rsmdLOG( "" << std::flush );
}

