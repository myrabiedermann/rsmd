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
#include "container/containerBase.hpp"
#include "container/atom.hpp"

#include <vector>
#include <algorithm>
#include <functional>

//
// molecule container
// 
// derived from ContainerBase
// contains atoms and all kind of useful methods that work with a molecule
//

class Molecule 
    : public ContainerBase<std::vector<Atom>>
{
    std::size_t molid    {0};
    std::string molname  {};

  public:
    //
    // getter/setter 
    //
    void        setID(std::size_t id)      { molid = id; }
    void        setName(std::string  name) { molname = name; }
    const auto& getID()      const { return molid; }
    const auto& getName()    const { return molname; }

    //
    // add new atoms to this molecule
    //
    inline auto addAtom(Atom a)                           { return data.emplace(end(), a); }
    inline auto addAtom(std::size_t id, std::string name) { auto it = data.emplace(end()); it->id = id; it->name = name; return it; }

    //
    // atom getters
    //
    const Atom& getAtom(std::size_t id) const;
   
    //
    // remove atoms from this molecule
    //
    void removeAtom(Atom& element);
    void removeAtom(std::size_t id);

    //
    // check if molecule contains a specific atom
    //
    bool containsAtom(Atom& element) const ;
    bool containsAtom(std::size_t id) const ;
    bool containsAtom(std::string name) const;

    //
    // check whether molecule contains any atoms
    //
    bool empty() const { return ( data.size() == 0 ? true : false ); }

    //
    // some useful operators
    //
    inline bool operator==(const Molecule& other) const { return this == std::addressof(other); }
	inline bool operator!=(const Molecule& other) const { return this != std::addressof(other); }
    inline bool operator<(const Molecule& other)  const { return this->getID() < other.getID(); }
    inline bool operator>(const Molecule& other)  const { return this->getID() > other.getID(); }

    friend inline std::ostream& operator << (std::ostream&, const Atom&);
};


inline std::ostream& operator<<(std::ostream& os, const Molecule& obj)
{
    os << "<Molecule: " << obj.getID() << ", " << obj.getName() << ", "
       << "contains " << obj.size() << " atoms>";
    
    return os;
}