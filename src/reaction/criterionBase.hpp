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

//
// a base class for reaction criterions
// like distances, angles etc
//
// -> derived from containerBase
// -> holds pairs of indices for each atom (molix, atomix)
//    between which the criterion is evaluated
// -> also holds minValue, maxValue and latestValue for the criterion
// -> has the (virtual) interface for a validity-check function, which will be 
//    implemented in the derived classes

class CriterionBase
    : public ContainerBase<std::vector<std::pair<std::size_t, std::size_t>>>
{
  protected:
    //
    // this is an abstract base class
    // and should only be derived from
    //
    CriterionBase() = default;

    REAL minValue {0};
    REAL maxValue {0};
    REAL latestValue{0};

    //
    // clone_impl(): the actual implementation of the clone() functionality
    // --> supposed to be implemented in the derived classes in order to 
    // guarantee correct deep-copying / cloning of objects 
    //
    virtual CriterionBase* clone_impl() const = 0;

  public:
    virtual ~CriterionBase() = default;

    // 
    // get type of criterion
    //
    virtual std::string getType() const = 0;

    //
    // get/set threshold values for the criterion
    //
    void setThresholds(const REAL& min, const REAL& max)  { minValue = min; maxValue = max; }
    void setThresholds(const std::pair<REAL, REAL>& values)  { minValue = values.first; maxValue = values.second; }
    void setMin(const REAL& value) { minValue = value; }
    void setMax(const REAL& value) { maxValue = value; }
    const auto& getMin() const { return minValue; }
    const auto& getMax() const { return maxValue; }

    //
    // get latest value of criterion
    //
    const auto&  getLatest() const { return latestValue; }

    //
    // setter for atom indices
    //
    void addAtomIndices(const std::size_t& molix, const std::size_t& atomix)
    {
        data.emplace_back( std::make_pair(molix, atomix) );
    }
    void addAtomIndices(const std::pair<std::size_t, std::size_t>& indices)
    {
        data.push_back( indices );
    }

    // 
    // check validity of criterion
    // (only to derive from)
    //
    virtual bool valid(const std::vector<Molecule>&, const REALVEC&) = 0;

    //
    // clone(): implements interface to clone an object of type CriterionBase
    //
    auto clone() const { return std::unique_ptr<CriterionBase>(clone_impl()); }
    // auto clone() const { return std::make_unique<CriterionBase>(clone_impl()); }


    friend inline std::ostream& operator << (std::ostream&, const CriterionBase&);
};


inline std::ostream& operator<<(std::ostream& os, const CriterionBase& obj)
{
    os << "<Criterion involving";
    
    std::string separator {" "};
    for( auto ixpair: obj.data )
    {
        os << separator << "(" << ixpair.first + 1 << ", " << ixpair.second + 1 << ")";
        separator = ", ";
    } 
    os << ">";
    
    return os;
}

