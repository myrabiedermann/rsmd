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

#include "control/controller.hpp"

std::atomic<int>  Controller::SIGNAL = {0};
std::atomic<bool> Controller::CIVILISED_SHUTDOWN = {false};


//
// signal handling
//
void Controller::signal( int SIG )
{
    static std::size_t gotCalled = 0;
    ++ gotCalled;

    if( gotCalled == 1 )
    {
        rsmdLOG( "Received signal " << SIG << " ... attempting civilised shutdown ..." );
    }
    else if( gotCalled == 2 )
    {
        rsmdLOG( "Received signal " << SIG << " ... still attempting civilised shutdown ..." );
    }
    else if( gotCalled == 3 )
    {
        rsmdLOG( "Received signal " << SIG << " ... IMMEDIATE SHUTDOWN!" );
        std::exit( SIG );
    }

    if( SIG == SIGUSR1 ) CIVILISED_SHUTDOWN.store(true);
    SIGNAL.store( SIG );
}


//
// setup 
//
void Controller::setup( int argc, char* argv[] )
{
    // opening words && handling of program options
    start_time = std::chrono::system_clock::now();
    std::time_t t1 = std::chrono::system_clock::to_time_t(start_time);
    parameters = std::make_unique<Parameters>(argc, argv);
    
    // log program options
    std::cout << "  [LOG]  " << "entering program rs@md, " << std::put_time( std::localtime(&t1), "%F %T" ) << '\n';
    std::cout << "  [LOG]  " << "reading the following program options ... \n" << parameters->str();
    
    // setup simulator and pass program options along
    switch( parameters->getSimulationAlgorithm() )
    {
        case SIMALGORITHM::MC:
            simulator = std::make_unique<SimulatorMetropolis>();
            assert(simulator);
            break;
 
        case SIMALGORITHM::RATE:
            simulator = std::make_unique<SimulatorRate>();
            assert(simulator);
            break;
    }
    simulator->setup(*parameters);

}


//
// start execution
//
void Controller::start() 
{    
    if( SIGNAL.load() == 0 )
    {
        simulator->run();
    }
}


//
// stop execution
//
void Controller::stop()
{
    // cleanup if required, write restart files if execution has been interrupted in a civilised manner (via SIGUSR1) etc. 
    if( CIVILISED_SHUTDOWN.load() )
    {
        std::cout << "  [LOG]   " << "civilised shutdown, catched SIGUSR1.\n";
        simulator->writeRestartFile(*parameters);
    }
    else if( SIGNAL.load() != 0 )
    {
        std::cout << "[WARNING] " << "not a civilised shutdown!\n";
        std::cout << "          ... you might need to do some cleaning up before attempting a restart.\n";
    }
    
    // finish up
    simulator->finish();

    // compute total run time
    end_time = std::chrono::system_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::milliseconds>( end_time - start_time );
    auto runtime_hours = int(runtime.count() / 3600000);
    auto runtime_minutes = int( (runtime.count() % 3600000) / 60000);
    auto runtime_seconds = std::round( (runtime.count() % 60000) / 1000);
    std::cout << "  [LOG]   total run time: " << std::setfill('0')          // set field fill character to '0'
              << runtime_hours 
              << "::"
              << std::setw(2) << runtime_minutes  
              << "::"
              << std::setw(2) << runtime_seconds     
              << " (hh::mm::ss)" << '\n';
    
    // compute <time per cycle>
    auto timePerCycle = runtime.count() / simulator->getNCycles();
    auto timePerCycle_hours = int(timePerCycle / 3600000);
    auto timePerCycle_minutes = int( (timePerCycle % 3600000) / 60000);
    auto timePerCycle_seconds = std::round( (timePerCycle % 60000) / 1000);
    std::cout << "  [LOG]   time per cycle: " << std::setfill('0') 
              << timePerCycle_hours << "::"
              << std::setw(2) << timePerCycle_minutes << "::"
              << std::setw(2) << timePerCycle_seconds << " (hh::mm::ss)\n";

    // closing words
    std::time_t t2 = std::chrono::system_clock::to_time_t(end_time);
    std::cout << "\n  [LOG]   leaving program rs@md, " <<  std::put_time( std::localtime(&t2), "%F %T" ) << '\n';
}