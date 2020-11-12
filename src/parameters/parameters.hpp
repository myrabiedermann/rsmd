#pragma once

#include "definitions.hpp"
#include "parser/reactionParser.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>

//
// parameter class
//
// ... can read program options from commandline or from file
// ... stores options in a boost::program_options::variables_map
// ... and implements getter functions for all options
//


enum ENGINE { NONE, GROMACS };
enum SIMMODE { NEW, RESTART };
enum SIMALGORITHM { RATE, MC };


class Parameters
{
  private:
    boost::program_options::variables_map parameterMap {};
    std::string     configFileName {};
    ENGINE          mdEngine {ENGINE::NONE};
    SIMMODE         simulationMode {SIMMODE::NEW};
    SIMALGORITHM    simulationAlgorithm {SIMALGORITHM::MC};

    std::string programName {"                     * * * * * * *\n                     *   rs@md   *\n                     * * * * * * *\n\n"};

    //
    // read program options from commandline / from file
    //
    void read(int, char* []);

    //
    // check for completeness / workability with given program options 
    //
    void check();

    //
    // formatted output
    //
    template<typename T>
    std::string formatted(const std::string&, const T&) const;

  public:
    Parameters(int, char* []);
    ~Parameters() = default;
    
    //
    // a universal getter for options in parameterMap
    //
    const auto& getOption(const std::string& s) const
    {
        if( parameterMap.count(s) )
            return parameterMap[s];  
        else
            throw std::logic_error("parameterMap does not contain " + s); 
    }

    //
    // some additional getters 
    //
    const auto& getEngineType() const { return mdEngine; }
    const auto& getSimulationMode() const { return simulationMode; }
    const auto& getSimulationAlgorithm() const {return simulationAlgorithm; }

    //
    // print options to string
    //
    std::string str() const;

};


template<typename T>
inline std::string Parameters::formatted(const std::string& name, const T& value) const
{
    std::stringstream stream {};
    stream << std::boolalpha << std::left
           << "    " << std::setw(30) << name << " = " << std::setw(30) << value;
    return stream.str();
}

template<>
inline std::string Parameters::formatted(const std::string& name, const std::vector<std::string>& values) const
{
    std::stringstream stream {};
    stream << std::boolalpha << std::left
           << "    " << std::setw(30) << name << " = ";
           
    std::string sep {""};
    for(const auto& val: values)
    {
        stream << sep << val;
        sep = ", ";
    } 

    return stream.str();
}