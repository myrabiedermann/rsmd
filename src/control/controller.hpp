#pragma once


#include "definitions.hpp"
#include "control/simulatorMetropolis.hpp"
#include "control/simulatorRate.hpp"

#include <csignal>

#include <iostream>
#include <iomanip>

#include <ctime>
#include <chrono>

#include <atomic>


//
// controller class
// 
// performs signal handling
// holds the simulator class
//


class Controller
{
  private:
    std::unique_ptr<SimulatorBase> simulator  {nullptr}; 
    std::unique_ptr<Parameters>    parameters {nullptr};
    std::chrono::system_clock::time_point start_time {};
    std::chrono::system_clock::time_point end_time   {};

  public:
    //
    // store signal in atomic
    //
    static std::atomic<int>  SIGNAL;
    static std::atomic<bool> CIVILISED_SHUTDOWN;

    //
    // a static signal handling function
    //
    static void signal(int SIG);

    // 
    // controlling methods
    //
    void setup( int argc, char* argv[] );
    void start();
    void stop();
};