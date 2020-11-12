
#pragma once

#include "reaction/criterionBase.hpp"
#include "enhance/math_utility.hpp"

#include <cassert>



// 
// distance criterion
//
class CriterionDistance
    : public CriterionBase
{
  protected: 
    //
    // clone_impl(): the actual implementation of the clone() functionality
    // --> supposed to be implemented in the derived classes in order to 
    // guarantee correct deep-copying / cloning of objects 
    //
    virtual CriterionDistance* clone_impl() const override
    {
        return new CriterionDistance(*this); 
    }

  public:
    virtual std::string getType() const override { return "distance"; }

    bool valid(const std::vector<Molecule>& reactants, const REALVEC& boxDimensions)
    {
        assert( data.size() == 2 );

        latestValue = enhance::distance(reactants[data[0].first](data[0].second),
                                        reactants[data[1].first](data[1].second), boxDimensions);
        return ( (latestValue >= minValue && latestValue <= maxValue) ? true : false );
    }
};



// 
// angle criterion
//
class CriterionAngle
    : public CriterionBase
{
  protected: 
    //
    // clone_impl(): the actual implementation of the clone() functionality
    // --> supposed to be implemented in the derived classes in order to 
    // guarantee correct deep-copying / cloning of objects 
    //
    virtual CriterionAngle* clone_impl() const override
    {
        return new CriterionAngle(*this); 
    }

  public:
    virtual std::string getType() const override { return "angle"; }

    bool valid(const std::vector<Molecule>& reactants, const REALVEC& boxDimensions)
    {
        assert( data.size() == 3 );

        latestValue = enhance::angle(reactants[data[0].first][data[0].second],
                                     reactants[data[1].first][data[1].second],
                                     reactants[data[2].first][data[2].second], boxDimensions);
        return ( (latestValue >= minValue && latestValue <= maxValue) ? true : false);
    }
};



// 
// dihedral criterion
//
class CriterionDihedral
    : public CriterionBase
{
  protected: 
    //
    // clone_impl(): the actual implementation of the clone() functionality
    // --> supposed to be implemented in the derived classes in order to 
    // guarantee correct deep-copying / cloning of objects 
    //
    virtual CriterionDihedral* clone_impl() const override
    {
        return new CriterionDihedral(*this); 
    }

  public:
    virtual std::string getType() const override { return "dihedral"; }

    bool valid(const std::vector<Molecule>& reactants, const REALVEC& boxDimensions)
    {
        assert( data.size() == 4 );

        latestValue = enhance::dihedral(reactants[data[0].first][data[0].second], 
                                        reactants[data[1].first][data[1].second], 
                                        reactants[data[2].first][data[2].second], 
                                        reactants[data[3].first][data[3].second], boxDimensions);
        return ( (latestValue >= minValue && latestValue <= maxValue) ? true : false);
    }
};