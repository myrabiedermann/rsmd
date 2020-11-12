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

#include "container/containerBase.hpp"
#include "container/molecule.hpp"

#include <vector>
#include <algorithm>
#include <numeric>

//
// topology container
//
// derived from ContainerBase
// contains molecules and all kind of useful methods that work with/on these molecules
//          + box dimensions
//

class Topology
    : public ContainerBase< std::vector<Molecule> >
{
    REALVEC dimensions {0, 0, 0};
    std::vector<std::pair<std::size_t, std::size_t>> reactedMoleculeRecords {};
    std::vector<std::pair<std::size_t, std::size_t>> reactedAtomRecords {};

  public:
    //
    // getter/setter for dimensions
    //
    inline void        setDimensions(const REALVEC& d) { dimensions = d; }
    inline const auto& getDimensions()    const { return dimensions; }

    //
    // getter/setter for reaction records
    //
    inline void        addReactionRecord(const std::size_t& molid) 
    { 
        reactedMoleculeRecords.emplace_back(std::make_pair(molid, 0)); 
    }
    inline const auto& getReactionRecordsAtoms()     { return reactedAtomRecords; }
    inline const auto& getReactionRecordsMolecules() { return reactedMoleculeRecords; }
    const std::size_t& getReactionRecordMolecule(const std::size_t& oldmolid);

    //
    // add new molecules to this topology
    //
    inline auto addMolecule(Molecule m)    
    { 
        return data.emplace(end(), m); 
    }
    inline auto addMolecule(std::size_t id, std::string name) 
    { 
        auto it = data.emplace(end()); 
        it->setID(id); 
        it->setName(name); 
        return it; 
    }

    //
    // remove specific molecule(s)
    //
    void removeMolecule(Molecule&);
    void removeMolecule(std::size_t);

    //
    // check if specific molecule exists
    //
    bool containsMolecule(const Molecule&) const;
    bool containsMolecule(const std::size_t&) const;
    
    //
    // get specific molecules
    //
    const Molecule& getMolecule(std::size_t) const;
    std::vector<std::reference_wrapper<Molecule>> getMolecules(std::string);

    // 
    // get specific molecule, create it if not yet existing
    //
    Molecule& getAddMolecule(std::size_t, std::string);

    //
    // get moleculetypes
    //
    std::vector<std::string> getMoleculetypes() const;

    //
    // get # of atoms
    //
    inline const auto getNAtoms() const 
    { 
        return std::accumulate( begin(), end(), 0, [](int counter, const auto& m){ return counter + m.size(); } ); 
    }

    //
    // sort topology, i.e. rearrange and renumber everything
    //
    void sort();

    //
    // repair molecule that is broken across periodic boundaries
    //
    void repairMoleculePBC(Molecule&);


    //
    // check whether topology contains any molecules
    //
    inline bool empty() const 
    { 
        return ( data.size() == 0 ? true : false ); 
    }


    //
    // clear topology
    //
    inline void clear() 
    { 
        data.clear(); 
        dimensions.setZero(); 
        reactedAtomRecords.clear(); 
        reactedMoleculeRecords.clear();
    }
    inline void clearReactionRecords() 
    { 
        reactedMoleculeRecords.clear(); 
        reactedAtomRecords.clear(); 
    }

    //
    // befriend << operator
    //
    friend inline std::ostream& operator<<(std::ostream& os, const Topology& obj);
};


inline std::ostream& operator<<(std::ostream& os, const Topology& obj)
{
    os << "<Topology contains " << obj.size() << " molecules within box dimensions " << obj.dimensions << ">";

    return os;
}