
#pragma once

#include "reaction/reactionBase.hpp"
#include "reaction/criterionDerived.hpp"


//
// a derived class to store a specific reaction candidate 
//
// derived from reactionBase and
// is initialised with a reactionBase object
// via copyConstructor 
//
// implements methods to update reactant molecules (set positions, velocities, update indices etc), 
// apply translation tables to create updated product molecules
//

class ReactionCandidate
    : public ReactionBase
{
  private:
    using ReactionBase::ReactionBase;

    //
    // write to stream
    //
    std::string str() const override;

  public:
    //
    // construct via copy constructor from ReactionBase object
    // of via copy construction from ReactionCandidate object
    // or via move constructor/move assignment
    //
    ReactionCandidate(const ReactionBase&);                                // copy constructor (from base class instance)
    ReactionCandidate(const ReactionCandidate&)             = default;     // copy constructor
    ReactionCandidate()                                     = delete;      // constructor
    ReactionCandidate& operator=(const ReactionCandidate&)  = delete;      // copy assignment
    ReactionCandidate(ReactionCandidate&&)                  = default;     // move constructor, required for std::vector::emplace_back()
    ReactionCandidate& operator=(ReactionCandidate&&)       = default;     // move assignment, required to shuffle std::vector<ReactionCandidate>
    ~ReactionCandidate()                                    = default;     // destructor
   
    //
    // get current reaction rate value / distance value 
    // for the first distance criterion
    //
    REAL getCurrentReactionRateValue() const;
    REAL getCurrentDistanceValue() const;

    // 
    // update reactant molecules (set positions / velocities
    // and update ID's)
    //
    void updateReactant(const std::size_t, const Molecule&);

    // 
    // apply transitionTables
    //
    void applyTransitions();
    
    // 
    // apply translationTables
    //
    void applyTranslations();

    //
    // check validity of all criterions
    //
    bool valid(const REALVEC&);

    //
    // write to stream - short version
    //
    std::string shortInfo() const;

};

