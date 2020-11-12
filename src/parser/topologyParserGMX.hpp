
#pragma once


#include "parser/topologyParserBase.hpp"
#include "enhance/utility.hpp"

#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <algorithm>
#include <filesystem>

//
// topology parser that reads/writes
// gromacs (GMX) topologies
// from .top & .gro files
//


class TopologyParserGMX : public TopologyParserBase
{
  private:
    std::string              systemName {};
    std::vector<std::string> topologyFileContent {};

    std::map<std::string, unsigned int> read_top( const std::string& );
    void read_gro( const std::string&, Topology&);
    void write_top(const std::string&, Topology&);
    void write_gro(const std::string&, Topology&);
    void write_index(const std::string&, const std::string&, Topology&);


  public:
    void read( Topology&, const std::size_t&);
    void readRelaxed( Topology&, const std::size_t&);
    void write(Topology&, const std::size_t&);

};
