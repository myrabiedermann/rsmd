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

#include "reactionBase.hpp"

//
// copy constructor
// 
ReactionBase::ReactionBase(const ReactionBase& other)
    : name(other.name)
    , reactants(other.reactants)
    , products(other.products)
    , transitionTables(other.transitionTables)
    , translationTables(other.translationTables)
    , reactionEnergy(other.reactionEnergy)
    , activationEnergy(other.activationEnergy)
    , reactionRate(other.reactionRate)
{
    for( auto& c: other.criterions )
    {
        criterions.push_back(c->clone());
    }
}



const auto ReactionBase::getReactant(const std::size_t& molid) const
{
    // attention: returns first molecule that matches molid (assumes that molid is unique)
    auto it = std::find_if( std::begin(reactants), std::end(reactants), [&molid](auto& m){ return molid == m.getID(); } );
    if( it == std::end(reactants) )   rsmdCRITICAL("couldn't find molecule in reactants: " << molid);
    return std::cref(*it);
}

const auto ReactionBase::getProduct(const std::size_t& molid) const
{
    // attention: returns first molecule that matches molid (assumes that molid is unique)
    auto it = std::find_if( std::begin(products), std::end(products), [&molid](auto& m){ return molid == m.getID(); } );
    if( it == std::end(products) )   rsmdCRITICAL("couldn't find molecule in products: " << molid);
    return std::cref(*it);
}


Molecule& ReactionBase::getAddReactant(const std::size_t& molid) 
{
    auto it = std::find_if( std::begin(reactants), std::end(reactants), [&molid](auto& m){ return molid == m.getID(); });
    if( it == std::end(reactants) ) 
    {
        it = reactants.emplace(std::end(reactants));
        it->setID(molid);
    }
    return std::ref(*it);
}

Molecule& ReactionBase::getAddProduct(const std::size_t& molid) 
{
    auto it = std::find_if( std::begin(products), std::end(products), [&molid](auto& m){ return molid == m.getID(); });
    if( it == std::end(products) ) 
    {
        it = products.emplace(std::end(products));
        it->setID(molid);
    }
    return std::ref(*it);
}



void ReactionBase::addTransition(const std::size_t& oldMolix, const std::size_t& oldix, const std::size_t& newMolix, const std::size_t& newix) 
{ 
    transitionTables.emplace_back(oldMolix, oldix, newMolix, newix);
}



void ReactionBase::addCriterion(const std::vector<std::pair<std::size_t, std::size_t>>& ixList, const std::pair<REAL, REAL>& thresholds) 
{ 
    auto it = criterions.end();
    if( ixList.size() == 2 )
    {
        it = criterions.emplace( std::end(criterions), std::make_unique<CriterionDistance>() );
    }
    else if( ixList.size() == 3 )
    {
        it = criterions.emplace( std::end(criterions),std::make_unique<CriterionAngle>() );
    }
    else if( ixList.size() == 4 )
    {
        it = criterions.emplace( std::end(criterions), std::make_unique<CriterionDihedral>() );
    }
    else
    {
        rsmdCRITICAL("no criterion involving more than 4 atoms has been implemented yet");
    }

    for(auto ix: ixList) it->get()->addAtomIndices(ix);
    it->get()->setThresholds(thresholds );
}


void ReactionBase::addTranslation(const std::vector<std::pair<std::size_t, std::size_t>>& indices, const REAL& value)
{
    assert(indices.size() == 2);
    translationTables.emplace_back(indices[0], indices[1], value);
}


//
// consistency check:
// - check that at least one reactant molecule is listed
// - check that at least one product molecule is listed
// - check that at least one criterion of type 'distance' is listed 
//   (and is the first criterion, i.e. index 0)
// - check if all atoms named in transitionTables or
//   within criterions actually exist within reactants/products
//
void ReactionBase::consistencyCheck() const
{
    // check reactants / products
    if( reactants.size() == 0 )
    {
        rsmdEXIT( "error in input: no reactant molecule was found" );
    }
    if( products.size() == 0 )
    {
        rsmdEXIT( "error in input: no product molecule was found" );
    }
    // check for distance criterion
    if( criterions[0]->getType() != "distance" )
    {
        rsmdEXIT( "error in input: the first listed criterion needs to be a distance" );
    }
    // check for consistency within reactants/products/transitionTables
    for( const auto& tt: transitionTables )
    {
        if( tt.oldMolix >= reactants.size() )
        {
            rsmdEXIT( "error in input directive [products]: given atom (" << tt.oldMolix + 1 << ", " << tt.oldix + 1 << ") doesn't exist in reactants");
        }
        else if( tt.oldix >= reactants[tt.oldMolix].size() )
        {
            rsmdEXIT( "error in input directive [products]: given atom (" << tt.oldMolix + 1 << ", " << tt.oldix + 1 << ") doesn't exist in reactants");
        }
        
        if( tt.newMolix >= products.size() )
        {
            rsmdEXIT( "error in input directive [products]: given atom (" << tt.newMolix + 1 << ", " << tt.newix + 1 << ") doesn't exist in products");
        }
        else if( tt.newix >= products[tt.newMolix].size() )
        {
            rsmdEXIT( "error in input directive [products]: given atom (" << tt.newMolix + 1 << ", " << tt.newix + 1 << ") doesn't exist in products");
        }
    }
    // check for consistency within reactants/products/movementTables
    for( const auto& tt: translationTables )
    {
        if( tt.indices1.first >= products.size() )
        {
            rsmdEXIT( "error in input directive [translations]: given atom (" << tt.indices1.first + 1 << ", " << tt.indices1.second + 1 << ") doesn't exist in products");
        }
        else if( tt.indices1.second >= products[tt.indices1.first].size() )
        {
            rsmdEXIT( "error in input directive [translations]: given atom (" << tt.indices1.first + 1 << ", " << tt.indices1.second + 1 << ") doesn't exist in products");
        }
        if( tt.indices2.first >= products.size() )
        {
            rsmdEXIT( "error in input directive [translations]: given atom (" << tt.indices2.first + 1 << ", " << tt.indices2.second + 1 << ") doesn't exist in products");
        }
        else if( tt.indices2.second >= products[tt.indices2.first].size() )
        {
            rsmdEXIT( "error in input directive [translations]: given atom (" << tt.indices2.first + 1 << ", " << tt.indices2.second + 1 << ") doesn't exist in products");
        }
    }
    // check for consistency within reactants/products/criterions
    for( const auto& criterion: criterions )
    {
        for( const auto& ixs: *criterion )
        {
            if( ixs.first >= reactants.size() )
            {
                rsmdEXIT( "error in input directive [criterions]: given atom (" << ixs.first + 1 << ", " << ixs.second + 1 << ") doesn't exist in reactants" );
            }
            else if ( ixs.second >= reactants[ixs.first].size() )
            {
                rsmdEXIT( "error in input directive [criterions]: given atom (" << ixs.first + 1 << ", " << ixs.second + 1 << ") doesn't exist in reactants" );
            }
        }
        
        if( criterion->getMin() >= criterion->getMax() )
        {
            rsmdEXIT( "error in input directive [criterions]: it seems that you have interchanged minimum and maximum value" );
        }
    }
}


//
// write info in a string
//
std::string ReactionBase::str() const
{
    std::stringstream stream;
    stream << "<Reaction '" << name << "', \n";
    
    stream << rsmdALL_formatting << "  reactants: ";
    for(const auto& mol: reactants)
    {
        stream << mol.getID() << " " << mol.getName() << ", ";
    }
    stream << '\n';

    stream << rsmdALL_formatting << "  products: ";
    for(const auto& mol: products)
        stream << mol.getID() << " " << mol.getName() << ", ";
    stream << '\n';

    // add +1 for indices in translation tables, movement tables and criterion templates 
    // to match numbering in input files
    stream << rsmdALL_formatting << "  transitions reactant -> product: ";
    for(const auto& tt: transitionTables)
        stream << "\n              (" << tt.oldMolix + 1 << ", " << tt.oldix + 1 << ") -> (" << tt.newMolix + 1 << ", " << tt.newix + 1 << ") ";
    stream << '\n';

    stream << rsmdALL_formatting << "  translational movements: ";
    for(const auto& mt: translationTables)
    {
        stream << "\n              (" << mt.indices1.first + 1 << ", " << mt.indices1.second + 1 << ") (" << mt.indices2.first + 1 << ", " << mt.indices2.second + 1 << ") " << mt.value;
    }
    stream << '\n';

    stream << rsmdALL_formatting << "  criterions: ";
    for(const auto& c: criterions)
    {
        stream << "\n              ";
        for(const auto& pair: *c)
            stream << "(" << pair.first + 1 << ", " << pair.second + 1 << ")   ";
        stream << "[" << c->getMin() << ", " << c->getMax() << "]";
    }
    stream << '\n';

    stream << rsmdALL_formatting << "  reaction energy: " << reactionEnergy << '\n';
    stream << rsmdALL_formatting << "  activation energy: " << activationEnergy << '\n';

    stream << rsmdALL_formatting << "  rate: ";
    for(const auto& r: reactionRate)
        stream << "\n              " << r.first << "  " << r.second;
    stream << '\n';

    stream << rsmdALL_formatting << ">";

    return stream.str();
}