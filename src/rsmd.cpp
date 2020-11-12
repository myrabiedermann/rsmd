/************************************************
 *                                              *
 *                rs@md                         *
 *    (reactive steps @ molecular dynamics )    *
 *                                              *
 ************************************************/


#include "control/controller.hpp"


int main( int argc, char* argv[] )
{
    // register important signals in Controller class
    // allowing civilised shutdown              // num on linux (num on osx), description
    std::signal( SIGHUP,  Controller::signal ); // 1, terminal line hangup
    std::signal( SIGINT,  Controller::signal ); // 2, interrupt program
    std::signal( SIGQUIT, Controller::signal ); // 3, quit program
    std::signal( SIGILL,  Controller::signal ); // 4, illegal instruction
    std::signal( SIGTRAP, Controller::signal ); // 5, trace trap
    std::signal( SIGIOT,  Controller::signal ); // 6, IOT Trap, synonym for SIGABRT
    std::signal( SIGABRT, Controller::signal ); // 6, abort program
    std::signal( SIGBUS,  Controller::signal ); // 7 (10), BUS error (bad memory access)
    std::signal( SIGFPE,  Controller::signal ); // 8, floating point exception
    std::signal( SIGUSR1, Controller::signal ); // 10, usr1
    std::signal( SIGSEGV, Controller::signal ); // 11, invalid memory reference
    std::signal( SIGTERM, Controller::signal ); // 15, software termination signal

    // simulation setup, execution & finish
    Controller controller; 
    controller.setup( argc, argv );
    controller.start();
    controller.stop(); 

    // return with correct status
    if( Controller::SIGNAL.load() != 0 )
        return Controller::SIGNAL.load();
    else
        return EXIT_SUCCESS;
}

