
#pragma once

#include "definitions.hpp"
#include "container/molecule.hpp"
#include "reaction/criterionDerived.hpp"

#include <string>
#include <vector>
#include <algorithm>
#include <iterator>


struct TransitionTable
{
    // while the rest of the program works with ID's as identifiers for molecules/atoms, 
    // here, real indices are used in order to access the correct atoms/molecules
    // in the reaction template
    std::size_t oldMolix {0};
    std::size_t oldix    {0};
    std::size_t newMolix {0};
    std::size_t newix    {0};

    TransitionTable(const std::size_t& ix1, const std::size_t& ix2, const std::size_t& ix3, const std::size_t& ix4)
     : oldMolix( ix1 )
     , oldix( ix2 )
     , newMolix( ix3 )
     , newix( ix4 ) 
    {}
};

struct TranslationTable
{
    // while the rest of the program works with ID's as identifiers for molecules/atoms, 
    // here, real indices are used in order to access the correct atoms/molecules
    // in the reaction template
    std::pair<std::size_t, std::size_t> indices1 {};
    std::pair<std::size_t, std::size_t> indices2 {};
    REAL        value  {0};

    TranslationTable(const std::pair<std::size_t, std::size_t>& ix1, const std::pair<std::size_t, std::size_t>& ix2, const REAL& val)
      : indices1(ix1)
      , indices2(ix2)
      , value(val)
    {}
};

//
// a base class to store a reaction template in
//
// holds: reactants (vector<Molecule>), products (<vector<Molecule>),
// transitionTables (one for each reactant),
// energy, reactionRateValues, criterions (vector<CriterionBase>)
//

class ReactionBase
{
  protected:
    std::string              name {};
    std::vector<Molecule>    reactants {};
    std::vector<Molecule>    products {};
    std::vector<TransitionTable> transitionTables {};
    std::vector<TranslationTable>   translationTables {};
    REAL                     reactionEnergy {0};
    REAL                     activationEnergy {0};
    std::vector<std::pair<REAL, REAL>> reactionRate {};
    std::vector<std::unique_ptr<CriterionBase>> criterions {};

    //
    // write info to a string
    // needs to be virtual because it is overwritten in derived class (ReactionCandidate)
    //
    virtual std::string str() const;
    
  public:
    ReactionBase() = default;      // constructor
    //
    // rule of five:
    // need to implement desctructor, copy assignment/constructor and move assignment/constructor
    // if required, due to criterions beeing unique_ptr
    //
    virtual ~ReactionBase() = default;      // destructor
    ReactionBase(const ReactionBase&);      // copy constructor
    ReactionBase(ReactionBase&&) = default; // move constructor
    ReactionBase& operator=(const ReactionBase&) = delete;  // copy assignment  
    ReactionBase& operator=(ReactionBase&&)      = default; // move assignment (required for std::swap(ReactionBase&, ReactionBase&))

    //
    // getters/setters
    //
    inline void         setName(const std::string& n) { name = n; }
    inline const auto&  getName()               const { return name; } 

    inline void         setReactionEnergy(const REAL& e) { reactionEnergy = e; }
    inline const auto&  getReactionEnergy()        const { return reactionEnergy; }
    
    inline void         setActivationEnergy(const REAL& e) { activationEnergy = e; }
    inline const auto&  getActivationEnergy()        const { return activationEnergy; }

    inline void         setRate( const std::vector<std::pair<REAL, REAL>> r ) { reactionRate = r; }
    inline const auto&  getRate()                                       const { return reactionRate; }

    const auto          getReactant(const std::size_t&) const;
    const auto&         getReactants()      const { return reactants; }
    auto&               getReactants()            { return reactants; }

    const auto          getProduct(const std::size_t&) const;
    const auto&         getProducts()       const { return products; }
    auto&               getProducts()             { return products; }


    // 
    // add stuff
    //
    Molecule& getAddReactant(const std::size_t&);
    Molecule& getAddProduct(const std::size_t&);
    
    void addTransition(const std::size_t&, const std::size_t&, const std::size_t&, const std::size_t&);

    void addCriterion(const std::vector<std::pair<std::size_t, std::size_t>>&, const std::pair<REAL, REAL>&);

    void addTranslation(const std::vector<std::pair<std::size_t, std::size_t>>&, const REAL&);

    void consistencyCheck() const;

    //
    // write to stream
    //
    friend inline std::ostream& operator << (std::ostream& stream, const ReactionBase& reaction)
    { 
        stream << reaction.str();
        return stream;

    }
    
};
