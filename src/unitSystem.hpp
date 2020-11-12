#pragma once

#include <string>
#include <sstream>


struct UnitSystem
{
    // units
    const std::string length;
    const std::string time;
    const std::string energy;
    const std::string temperature;

    // getter for natural constants
    double getR() const { return R; }

    // constructor
    UnitSystem(const std::string& l, const std::string& t, const std::string& e, const std::string& temp)
        : length(l)
        , time(t)
        , energy(e)
        , temperature(temp)
    {
        if( energy == "kJ/mol" || energy == "kJ" )
        {
            R = 0.00831446261815324;    // R in kJ/mol
        }
        else if( energy == "kCal/mol" || energy == "kCal" )
        {
            R = 0.00198720425864083;    // R in kCal/mol
        }
        else
        {
            std::stringstream message {};
            message << "could not understand given energy unit: " << energy;
            throw std::runtime_error(message.str());
        }
    }

  private: 
    // natural constants
    double R {0}; 
    
};