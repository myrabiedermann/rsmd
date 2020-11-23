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

#include "reaction/reactionCandidate.hpp"

//
// construct via copy constructor from ReactionBase object
// 
ReactionCandidate::ReactionCandidate(const ReactionBase& other)
    : ReactionBase(other)
{}



//
// get current reaction rate value
//
REAL ReactionCandidate::getCurrentReactionRateValue() const
{
    // assumes that the first criterion is the correct distance criterion
    // the fact that it is a "distance" criterion is never checked, 
    // but how to assure that it is the 'correct' one?
    // --> add notes for users ...
    REAL currentDistanceValue = criterions[0]->getLatest(); 
    REAL currentRateValue = reactionRate[0].second;    
    for( const auto& pair: reactionRate )
    {
        if( pair.first > currentDistanceValue )    break;
        currentRateValue = pair.second;
    }
    return currentRateValue;
}


//
// get current distance value
// (for first distance criterion)
//
REAL ReactionCandidate::getCurrentDistanceValue() const
{
    return criterions[0]->getLatest();
}


// 
// update reactant molecules (set positions / velocities
// and update ID's)
//
// assumes that atoms in reactant and 'real' molecule are listed in exactly the same order!
void ReactionCandidate::updateReactant(const std::size_t reactantix, const Molecule& molecule)
{   
    Molecule& reactant = reactants[reactantix];

    for(auto& reactantAtom: reactant)
    {
        const Atom& atom = molecule[reactantAtom.id - 1];
        reactantAtom.id = atom.id;
        reactantAtom.position = atom.position;
        reactantAtom.velocity = atom.velocity;
    }
    reactant.setID( molecule.getID() );  

    rsmdDEBUG("updated reactant molecule " << reactantix + 1);
    rsmdDEBUG(reactant);
    #ifndef NDEBUG
        for(const auto& a: reactant) rsmdDEBUG(a);
    #endif
}



// 
// apply transitionTables
//
void ReactionCandidate::applyTransitions()
{
    for(auto& tt: transitionTables)
    {
        products[tt.newMolix](tt.newix).id = reactants[tt.oldMolix](tt.oldix).id;
        products[tt.newMolix](tt.newix).position = reactants[tt.oldMolix](tt.oldix).position;
        products[tt.newMolix](tt.newix).velocity = reactants[tt.oldMolix](tt.oldix).velocity;
    }

    #ifndef NDEBUG
    rsmdDEBUG("performed transition from reactants -> products");
    for(const auto& product: products)
    {
        rsmdDEBUG(product);
        for(const auto& atom: product) rsmdDEBUG(atom);
    }
    #endif
}



// 
// apply translationTables
//
void ReactionCandidate::applyTranslations()
{
    for(const auto& tt: translationTables)
    {
        rsmdLOG("... performing translation for product atom: " 
                << products[tt.indices1.first](tt.indices1.second).name
                << " towards/away from "
                << products[tt.indices2.first](tt.indices2.second).name );
        // first: compute connection vector between atoms
        REALVEC vector = products[tt.indices2.first](tt.indices2.second).position 
                       - products[tt.indices1.first](tt.indices1.second).position;
		REAL distance = vector.norm();
        vector /= vector.norm();
        rsmdDEBUG( "    position before: " << products[tt.indices1.first](tt.indices1.second).position );
        rsmdLOG( "    distance before: " << distance );
		// second: translate atom along the connection vector
        products[tt.indices1.first](tt.indices1.second).position += tt.value * vector;
        vector = products[tt.indices2.first](tt.indices2.second).position
			   - products[tt.indices1.first](tt.indices1.second).position;
		distance = vector.norm();
		rsmdDEBUG( "    position after: " << products[tt.indices1.first](tt.indices1.second).position );
		rsmdLOG( "    distance after: " << distance );
    }
}



//
// check validity of all criterions
//
bool ReactionCandidate::valid(const REALVEC& boxDimensions)
{
    rsmdDEBUG("checking validity of all criterions ...");
    for( const auto& criterion: criterions )
    {
        rsmdDEBUG(*criterion);
        if( ! criterion->valid(reactants, boxDimensions) )
        {
            rsmdDEBUG( "... INVALID: " << criterion->getLatest() << " not in [" << criterion->getMin() << ", " << criterion->getMax() << "]" );
            rsmdDEBUG( "... skipping any further criterions" );
            rsmdDEBUG(" ");
            return false;
        } 
        rsmdDEBUG( "... VALID: " << criterion->getLatest() << " is in [" << criterion->getMin() << ", " << criterion->getMax() << "]" )
            
    }
    rsmdDEBUG( "... all criterions are valid!" );
    rsmdDEBUG(" ");
    return true;
}



//
// write to string - short version
//
std::string ReactionCandidate::shortInfo() const
{
    std::stringstream stream;
    stream << "<Reaction " << name << ", \n";

    stream << rsmdALL_formatting << "  reactants: ";
    std::string separator {""};
    for(const auto& mol: reactants)
    {
        stream << separator << mol.getID() << " " << mol.getName();
        separator = ", ";
    }
    stream << '\n';

    stream << rsmdALL_formatting << "  products: ";
    separator = "";
    for(const auto& mol: products)
    {
        stream << separator << mol.getID() << " " << mol.getName();
        separator = ", ";
    }
    stream << ">";

    return stream.str();
}


//
// write info in a string
//
std::string ReactionCandidate::str() const 
{ 
    std::stringstream stream;
    stream << "<Reaction '" << name << "', \n";
    
    stream << rsmdALL_formatting << "  reactants: ";
    for(const auto& mol: reactants)
        stream << mol.getID() << " " << mol.getName() << ", ";
    stream << '\n';

    stream << rsmdALL_formatting << "  products: ";
    for(const auto& mol: products)
        stream << mol.getID() << " " << mol.getName() << ", ";
    stream << '\n';

    stream << rsmdALL_formatting << "  criterions: ";
    for(const auto& c: criterions)
    {
        std::string separator = "\n              ";
        for(const auto& ixPair: *c)
        {
            stream << separator << "reactant " << ixPair.first + 1 << " atom " << ixPair.second + 1;
            separator = ", ";
        }
        stream << '\n' << rsmdALL_formatting << "    |-> " << (c->getLatest() <= c->getMax() && c->getLatest() >= c->getMin() ? "valid" : "not valid") << "  (value: " << c->getLatest() << ", thresholds " << c->getMin() << " - " << c->getMax() << ')';
    }
    stream << '\n';

    stream << rsmdALL_formatting << ">";

    return stream.str();
}


