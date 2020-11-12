
#include "container/molecule.hpp"


//
// atom getters
//
const Atom& Molecule::getAtom(std::size_t id) const
{ 
    rsmdDEBUG( "searching for atom " << id << " within molecule " << this->molname << ", " << this->molid );
    auto it = std::find_if( begin(), end(), [&id](auto& a){ return id == a.id; } );
    if( it == end() )   rsmdCRITICAL("couldn't find atom " << id << " in molecule " << this->molid);
    return std::cref(*it);
}


//
// remove atoms from this molecule
//
void Molecule::removeAtom(Atom& element)
{
    data.erase( std::remove_if( begin(), end(), [&](auto& a){ return element == a; } ), end() );
}

void Molecule::removeAtom(std::size_t id)
{
    // attention: removes all atoms that have identifier id
    data.erase( std::remove_if(begin(), end(), [&](auto& a){ return id == a.id; }), end() );
}


//
// check if molecule contains a specific atom
//
bool Molecule::containsAtom(Atom& element) const
{
    auto it = std::find_if( begin(), end(), [&](auto& a){ return element == a; } );
    return (it == end() ? false : true);
}

bool Molecule::containsAtom(std::size_t id) const
{
    auto it = std::find_if( begin(), end(), [&](auto& a){ return id == a.id; } );
    return (it == end() ? false : true);
}

bool Molecule::containsAtom(std::string name) const
{
    auto it = std::find_if( begin(), end(), [&](auto& a){ return name == a.name; } );
    return (it == end() ? false : true);
}
