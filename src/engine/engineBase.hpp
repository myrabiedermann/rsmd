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
#include "parameters/parameters.hpp"

#include <stdlib.h>
#include <sstream>
#include <unistd.h>
#include <sys/wait.h>
#include <csignal>
#include <cstring>

//
// a base class that implements
// the interface for all md engines
//

enum PIPE_FILE_DESCRIPTORS
{
    READ_FD = 0,
    WRITE_FD = 1
};

class EngineBase
{
  protected:
    EngineBase() = default;

    // execute a command in subprocess 
    template<typename... Args>
    void execute( const char*, Args&& ... args );

    // execute a command in subprocess with piped input
    template<typename... Args>
    void execute( std::string&, const char*, Args&& ... args );

  public:
    virtual ~EngineBase() = default;

    virtual void setup(const Parameters&) = 0;
    virtual void verifyExecutable() = 0;
    virtual void runMD( const std::size_t& ) = 0;
    virtual void runMDInitial( ) = 0;
    virtual void runMDAppending( const std::size_t&, const std::size_t& ) = 0;
    virtual bool runRelaxation( const std::size_t& ) = 0;
    virtual void runEnergyComputation( const std::size_t&, const std::size_t& ) = 0;
    virtual void cleanup( const std::size_t&) = 0;
};



//
// execute the command (+ cmdline options) given by args
//
template<typename... Args>
void EngineBase::execute( const char* cmd, Args&& ... args )
{
    std::string empty {};
    execute(empty, cmd, std::forward<Args>(args)...);
}



//
// execute the command (+ cmdline options) given by args with piped input (parent->child)
//
template<typename... Args>
void EngineBase::execute( std::string& pipeIn, const char* cmd, Args&& ... args )
{
    int status;
    int childIn[2], childOut[2];    // for piping input/output of the child process
    pid_t child_pid;
    char buffer[101] = {0};
    ssize_t readResult {0};
    std::string pipeOut {};

    // about the pipes:
    // the first integer in the respective fd (file descriptor) array (element 0) is set up and opened for reading, 
    // while the second integer (element 1) is set up and opened for writing. 
    // visually speaking, the output of fd[1] becomes the input for fd[0]. 
    // once again, all data traveling through the pipe moves through the kernel.

    // some verbosity:
    std::stringstream stream {};
    stream << "[EngineBase::execute()] running:";
    if( ! pipeIn.empty() )
    {
        std::string tmp {pipeIn}; 
        tmp.erase(std::remove(tmp.begin(), tmp.end(), '\n'), tmp.end());
        stream << " " << tmp << " |";
    }
    using expander = int[];
    (void) expander {0, (void(stream << ' ' << std::forward<const char*>(args)),0)...};
    rsmdDEBUG( stream.str() );
    
    // creating the pipes
    if( pipe(childIn) < 0)   // failure in creating a pipe
        rsmdCRITICAL( "failure in creating a pipe")
    if( pipe(childOut) < 0)  // failure in creating a pipe
        rsmdCRITICAL( "failure in creating a pipe" );

    // try forking a new process
    child_pid = fork();     
    // about fork():
    // returns -1 in case of failure,
    //          0 to the newly created child process and 
    //          child_pid (positive value) to the parent
    switch( child_pid )
    {
        case -1:    // fork failed
            rsmdCRITICAL( "fork failed" );
            break;
    
        case 0:     // fork successful, this section is only entered from within child
            // replace standard input with input part of pipe and close unused other half
            dup2(childIn[READ_FD], READ_FD); // dup2() includes closing of this file descriptor! 
            close(childIn[WRITE_FD]);

            // replace standard output with output part of pipe and close unused other half
            dup2(childOut[WRITE_FD], STDOUT_FILENO);
            dup2(childOut[WRITE_FD], STDERR_FILENO);            
            close(childOut[READ_FD]);

            // execute program:
            execlp( cmd, std::forward<const char*>(args)..., NULL );     // the program should receive its own command as argv[0]
            
            // should't return, so exit here
            std::exit(EXIT_FAILURE);

            break;
    
        default:    // fork successful, this section is only entered from within parent
            // write to childIn[1] (which will become stdin for the child) and close both sides of the pipe afterwards
            close( childIn[READ_FD] );
            write( childIn[WRITE_FD], pipeIn.c_str(), (strlen( pipeIn.c_str() )+1) );
            close( childIn[WRITE_FD]) ;

            // wait for state changes in the child (aka for it to return)
            // and handle exit status or any signals correctly
            do
            {
                wait( &status );
                if( WIFEXITED(status) )
                {
                    rsmdDEBUG( "[EngineBase::execute()] " << "exited: status = " << WEXITSTATUS(status) );
                }
                else if( WIFSIGNALED(status) )
                {
                    rsmdWARNING( "[EngineBase::execute()] " << "killed by signal " << WTERMSIG(status) );    
                    std::raise( WTERMSIG(status) ); 
                }
                else if( WIFSTOPPED(status) ) 
                {
                    rsmdWARNING( "[EngineBase::execute()] " << "stopped by signal " << WSTOPSIG(status) );
                    std::raise( WSTOPSIG(status) );
                }
                else if( WIFCONTINUED(status) )
                {
                    rsmdWARNING( "[EngineBase::execute()] " << "continued" );
                }
            } while( !WIFEXITED(status) && !WIFSIGNALED(status) );
            
            // close writing part of file descriptor
            close( childOut[WRITE_FD] );

            // get output from reading part of file descriptor
            while( (readResult = read(childOut[READ_FD], buffer, sizeof(buffer)-1)) > 0 )
            {
                buffer[readResult] = 0;
                pipeOut.append( buffer );
            }
        
            // finally also close reading part
            close(childOut[READ_FD]);

            // spill process output in case anything went wrong
            if( WEXITSTATUS(status) != 0 || WIFSIGNALED(status) || WIFSTOPPED(status) || WIFCONTINUED(status) )
            {
                rsmdWARNING( "process output was: \n" << pipeOut );
            }

            // throw exception in case exited with status != 0
            if( WEXITSTATUS(status) != 0 )
            {
                throw std::runtime_error("something went wrong in child process execution");
            }

            break;
    }   // end of switch()
}




