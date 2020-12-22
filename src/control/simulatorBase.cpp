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

#include "control/simulatorBase.hpp"

//
// a base class version of setup()
// --> implements all the general stuff
//
void SimulatorBase::setup(const Parameters& parameters)
{
    // ... of random engine
    const auto seed = parameters.getOption("rseed").as<std::size_t>();
    if( seed != 0 )
    {
        enhance::RandomEngine.setSeed( seed );
        rsmdLOG( "... setting random seed to " << seed);
    }
    else
    {
        rsmdLOG( "... using (true) random seed " << enhance::RandomEngine.getSeed() );
    }

    // ... of the mdEngine and energyParser
    switch( parameters.getEngineType() )
    {
        case ENGINE::GROMACS:   
            mdEngine = std::make_unique<EngineGMX>();
            energyParser = std::make_unique<EnergyParserGMX>();
            assert(mdEngine);
            assert(energyParser);
            mdEngine->setup(parameters);
            energyParser->setup(parameters);

            unitSystem = std::make_unique<UnitSystem>("nm", "ps", "kJ/mol", "K");
            assert(unitSystem);
            
            break;

        case ENGINE::NONE:
            rsmdCRITICAL( "md engine is set to none" );
            break;
    }

    // ... of the universe
    universe.setup(parameters);  

    // ... handle start/end cycle nr according to simulationMode
    nCycles = parameters.getOption("simulation.cycles").as<std::size_t>();
    switch( parameters.getSimulationMode() )
    {
        case SIMMODE::NEW:
            rsmdLOG( "... will start a new simulation from cycle = 0" );
            // ... statistics file
            if( parameters.getOption("statistics").as<std::string>() != "" )
            {
                STATISTICS_FILE.open( parameters.getOption("statistics").as<std::string>() );
                if( ! STATISTICS_FILE )
                {   // safety check
                    rsmdCRITICAL( "opening file " << parameters.getOption("statistics").as<std::string>() << " failed.")
                } 
            }
            break;
        
        case SIMMODE::RESTART:
            lastReactiveCycle = parameters.getOption("simulation.restartCycleFiles").as<std::size_t>();
            currentCycle = parameters.getOption("simulation.restartCycle").as<std::size_t>();
            rsmdLOG( "... will restart simulation from cycle = " << currentCycle );
            // ... statistics file
            STATISTICS_FILE.open( parameters.getOption("statistics").as<std::string>(), std::ostream::app );
            if( ! STATISTICS_FILE )
            {   // safety check
                rsmdCRITICAL( "opening file " << parameters.getOption("statistics").as<std::string>() << " failed.")
            } 
            break;
    }
}



//
// run simulation
//
void SimulatorBase::run()
{
    // initial md sequence
    if( currentCycle == 1 )
    {
        rsmdLOG("@ cycle 0 (initial md sequence)");
        mdEngine->runMDInitial();
    }

    while( currentCycle <= nCycles )
    {
        // check for signals
        if( Controller::SIGNAL.load() != 0 ) break;

        rsmdLOG("@ cycle " << currentCycle);
        rsmdDEBUG("@ cycle " << currentCycle);

        // reactive step
        reactiveStep();
        
        // check for signals
        if( Controller::SIGNAL.load() != 0 && ! Controller::CIVILISED_SHUTDOWN.load() ) break;

        // do md sequence
        mdSequence();

        ++ currentCycle;
        ++ nCyclesCompleted;
        
        rsmdLOG(std::flush);
    }
}



//
// do md sequence
//
void SimulatorBase::mdSequence()
{
    if( lastReactiveCycle == currentCycle )
    {
        mdEngine->runMD(currentCycle);
    }
    else
    {
        mdEngine->runMDAppending(currentCycle, lastReactiveCycle);
    }
}



//
// write a restart file
//
void SimulatorBase::writeRestartFile(const Parameters& parameters) const
{
    std::ofstream FILE( parameters.getOption("output").as<std::string>() );
    if( FILE.bad() ){
        rsmdCRITICAL("something went wrong with outstream to " << parameters.getOption("output").as<std::string>() );
    }
    else
    {
        rsmdLOG( "... writing program options for restarting to " << parameters.getOption("output").as<std::string>() );
    } 

    // [simulation]
    FILE << "[simulation]\n";
    FILE << "engine      = " << parameters.getOption("simulation.engine").as<std::string>() << '\n';
    FILE << "cycles      = " << parameters.getOption("simulation.cycles").as<std::size_t>() << '\n';
    FILE << "restart     = " << "on" << '\n';
    FILE << "restartCycle = " << currentCycle << '\n';
    FILE << "restartCycleFiles = " << lastReactiveCycle << '\n';
    FILE << '\n';

    // [reaction]
    FILE << "[reaction]\n";
    for( const auto& filename: parameters.getOption("reaction.file").as<std::vector<std::string>>() )
        FILE << "file        = " << filename << '\n';
    FILE << "mc          = " << (parameters.getOption("reaction.mc").as<bool>() ? "on" : "off") << '\n';
    FILE << "rate        = " << (parameters.getOption("reaction.rate").as<bool>() ? "on" : "off") << '\n';
    if( parameters.getOption("reaction.rate").as<bool>() )
    {
        FILE << "frequency   = " << parameters.getOption("reaction.frequency").as<REAL>() << '\n';
    }
    else
    {
        FILE << "temperature = " << parameters.getOption("reaction.temperature").as<REAL>() << '\n';
        FILE << "averagePotentialEnergy = " << parameters.getOption("reaction.averagePotentialEnergy").as<REAL>() << '\n';
        FILE << "computeLocalPotentialEnergy = " << (parameters.getOption("reaction.computeLocalPotentialEnergy").as<bool>() ? "on" : "off" ) << '\n';
        FILE << "computeSolvationPotentialEnergy = " << (parameters.getOption("reaction.computeSolvationPotentialEnergy").as<bool>() ? "on" : "off" ) << '\n';
    }
    FILE << "saveRejected = " << (parameters.getOption("reaction.saveRejected").as<bool>() ? "on" : "off") << '\n';
    FILE << '\n';

    // md engine related --> [gromacs], ...
    switch( parameters.getEngineType() )
    {
        case ENGINE::GROMACS:   
            FILE << "[gromacs]\n";
            FILE << "topology     = " << std::to_string(lastReactiveCycle) + ".top" << '\n'; 
            FILE << "coordinates  = " << std::to_string(lastReactiveCycle) + "-md.gro" << '\n';
            FILE << "mdp          = " << parameters.getOption("gromacs.mdp").as<std::string>() << '\n';
            FILE << "mdp.relaxation = " << parameters.getOption("gromacs.mdp.relaxation").as<std::string>() << '\n';
            if( parameters.getOption("reaction.computeSolvationPotentialEnergy").as<bool>() )
                FILE << "mdp.energy   = " << parameters.getOption("gromacs.mdp.energy").as<std::string>() << '\n';
            FILE << "backup       = " << (parameters.getOption("gromacs.backup").as<bool>() ? "on" : "off") << '\n';
            FILE << "nt           = " << parameters.getOption("gromacs.nt").as<int>() << '\n';
            FILE << "ntmpi        = " << parameters.getOption("gromacs.ntmpi").as<int>() << '\n';
            FILE << "ntomp        = " << parameters.getOption("gromacs.ntomp").as<int>() << '\n';
            break;

        case ENGINE::NONE:
            break;
    }

    FILE.close();
}
