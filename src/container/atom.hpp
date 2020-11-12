
#pragma once

#include "definitions.hpp"

#include <iostream>
#include <iomanip>
#include <string>


//
// atom container 
//

struct Atom
{
    std::size_t id          {0};
    std::string name        {};
    REALVEC     position    {0, 0, 0};
    REALVEC     velocity    {0, 0, 0};

	inline bool operator==(const Atom& other) const { return this == std::addressof(other); }
	inline bool operator!=(const Atom& other) const { return this != std::addressof(other); }
    inline bool operator<(const Atom& other)  const { return this->id < other.id; }
    inline bool operator>(const Atom& other)  const { return this->id > other.id; }

    friend inline std::ostream& operator << (std::ostream&, const Atom&);
};



inline std::ostream& operator<<(std::ostream& os, const Atom& obj)
{
    os << "<Atom: " << obj.id << ", " << obj.name << ", "
       << "positions: " << obj.position << ", "
       << "velocities: " << obj.velocity
       << ">";
    
    return os;
}