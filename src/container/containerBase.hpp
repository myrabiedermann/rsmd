
#pragma once

#include <iterator>


//
// a container base class to derive from
//
// implements all iterator-related functions like begin() / end() etc.
// for use in loops
// as well as container-related operators like [] and ()
//


template<typename T>
struct ContainerBase
{
    T data {};

    inline auto& operator()(std::size_t i)                 { return data[i]; };
    inline constexpr auto& operator()(std::size_t i) const { return data[i]; };
   
    inline auto& operator[](std::size_t i)                 { return data[i]; };
    inline constexpr auto& operator[](std::size_t i) const { return data[i]; };
 
    inline auto begin()         { return std::begin(data); };
    inline auto end()           { return std::end(data); };
 
    inline auto begin()   const { return std::begin(data); };
    inline auto end()     const { return std::end(data); };
 
    inline auto cbegin()  const { return std::cbegin(data); };
    inline auto cend()    const { return std::cend(data); };
 
    inline auto rbegin()        { return std::rbegin(data); };
    inline auto rend()          { return std::rend(data); };
 
    inline auto rbegin()  const { return std::rbegin(data); };
    inline auto rend()    const { return std::rend(data); };
 
    inline auto crbegin() const { return std::crbegin(data); };
    inline auto crend()   const { return std::crend(data); };

    inline auto size()    const { return data.size(); };

    inline auto&       front()         { return *this->begin(); }
    inline const auto& front()   const { return *this->begin(); }
    inline auto&       back()          { auto tmp = this->end(); --tmp; return *tmp; }
    inline const auto& back()    const { auto tmp = this->end(); --tmp; return *tmp; }

    virtual ~ContainerBase() = default;

  protected:
    ContainerBase() = default;

};