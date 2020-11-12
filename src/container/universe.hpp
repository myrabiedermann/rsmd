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

#include "unitSystem.hpp"
#include "enhance/random.hpp"
#include "container/topology.hpp"
#include "reaction/reactionCandidate.hpp"
#include "parser/topologyParserGMX.hpp"
#include "parser/reactionParser.hpp"

//
// universe class
//
// holds:   - topologies (+parser) and
//          - reaction templates (+parser) 
//
// can do:  - search for reaction candidates in topology
//          - react a given candidate
//

class Universe
{
  private:  
    // topology related stuff
    Topology topologyOld {};
    Topology topologyNew {};
    Topology topologyRelaxed {};
    std::unique_ptr<TopologyParserBase> topologyParser {nullptr};

    // reaction related stuff
    std::vector<ReactionBase> reactionTemplates {};
    
    std::unique_ptr<UnitSystem> unitSystem {nullptr};

    //
    // repair a molecule in case it is broken across periodic boundaries
    //
    void makeMoleculeWhole(Molecule&, const REALVEC& dimensions);

  public:
    //
    // initial setup of the universe
    //
    void setup(const Parameters&); 

    //
    // update universe (at beginning of new cycle)
    //
    void update(const std::size_t&);

    //
    // write function
    //
    void write(const std::size_t&);

    //
    // read structure after relaxation
    //
    void readRelaxed(const std::size_t&);

    //
    // search for reaction candidates
    //
    std::vector<ReactionCandidate> searchReactionCandidates();

    //
    // check availability of given candidate
    //
    bool isAvailable(const ReactionCandidate&);

    //
    // react a given candidate
    //
    void react(ReactionCandidate&);

    //
    // check a given candidate for for 'physical meaningfulness'
    //
    void checkMovement(const ReactionCandidate&);
    
    //
    // some getters
    //
    const auto& getReactionTemplates() const { return reactionTemplates; }
    
};