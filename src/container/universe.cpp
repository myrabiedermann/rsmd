
#include "container/universe.hpp"

//
// initial setup of the universe
//
void Universe::setup(const Parameters& parameters)     
{
    switch( parameters.getEngineType() )
    {
        case ENGINE::GROMACS:   
            topologyParser = std::make_unique<TopologyParserGMX>();
            assert(topologyParser);

            unitSystem = std::make_unique<UnitSystem>("nm", "ps", "kJ/mol", "K");
            assert(unitSystem);

            break;

        case ENGINE::NONE:
            rsmdCRITICAL( "md engine is set to none" );
            break;
    }

    // read reaction templates from files
    auto reactionFiles = parameters.getOption("reaction.file").as<std::vector<std::string>>();
    rsmdLOG("... reading reaction templates ... ");
    ReactionParser reactionParser {};
    for( const auto& file: reactionFiles )
    {
        auto reaction = reactionParser.read( file );
        
        // some verbose printing: 
        rsmdLOG( "... from file '" << file << "': ")
        rsmdLOG( reaction );

        rsmdLOG( "... checking for consistency in provided input for reaction '" << reaction.getName() << "' ...");
        // check that reaction template contains required input 
        // for the chosen simulation algorithm
        switch( parameters.getSimulationAlgorithm() )
        {
            case SIMALGORITHM::MC:
                // --> check for energy != 0
                if( reaction.getReactionEnergy() == 0 )
                    rsmdWARNING( "    reaction energy == 0, are you sure that is correct?" );
                break;

            case SIMALGORITHM::RATE:
                // --> check for reaction rate
                if( reaction.getRate().size() == 0 )
                    rsmdWARNING( "    no reaction rate input, are you sure that is correct?" );
                break;
        }
        // check for consistency within reactants/products/criterions
        reaction.consistencyCheck();
        rsmdLOG( "... consistency check done. everything seems fine.");
        
        reactionTemplates.emplace_back(reaction);
    }
}


// 
// update topologies
//
void Universe::update(const std::size_t& cycle) 
{
    topologyOld.clear();
    topologyNew.clear();
    topologyRelaxed.clear();

    topologyParser->read(topologyOld, cycle);
    topologyOld.clearReactionRecords();
    topologyNew = topologyOld;
}


//
// write (new) topology to file 
//
void Universe::write(const std::size_t& cycle)
{
    topologyNew.sort();
    topologyParser->write(topologyNew, cycle);
}


//
// read relaxed configuration from file
//
void Universe::readRelaxed(const std::size_t& cycle)
{
    topologyRelaxed.clear();
    topologyParser->readRelaxed(topologyRelaxed, cycle);
}


//
// check given reaction candidate (in topologyRelaxed) for 
// 'physical meaningfulness' after relaxation, 
// i.e. how much the corresponding atoms moved
//
void Universe::checkMovement(const ReactionCandidate& candidate)
{
    // first: compute typical length in system against which to check
    REAL volume = topologyNew.getDimensions()[0] * topologyNew.getDimensions()[1] * topologyNew.getDimensions()[2];
    // REAL typicalDistance = std::sqrt( (3.0 * volume) / (4.0 * M_PI * topologyNew.getNAtoms()) );
    REAL typicalDistance = std::cbrt( (3.0 * volume) / (4.0 * M_PI * topologyNew.getNAtoms()) );

    for( auto& molecule: candidate.getProducts() )
    {
        // get same molecule in topologyRelaxed
        std::size_t newMolID = topologyNew.getReactionRecordMolecule(molecule.getID());
        auto& newMolecule = topologyRelaxed.getMolecule(newMolID);

        // go through molecule and compute movement of atoms
        auto atomBefore = molecule.begin();
        auto atomAfter  = newMolecule.begin();
        while( atomBefore != molecule.end() || atomAfter != newMolecule.end() )
        {
            auto distance = enhance::distance(*atomBefore, *atomAfter, topologyNew.getDimensions());

            if( distance > 3 * typicalDistance )
            {
                rsmdWARNING( std::setprecision(3) << "... atom " << atomAfter->name << " " << atomAfter->id << " of molecule " << newMolecule.getName() << " " << newMolecule.getID() 
                        << " moved more than three times the typical distance: " << distance << ' ' << unitSystem->length << " ( > 3 * " << typicalDistance << ' ' << unitSystem->length << ")");
            }
            else if( distance > 2 * typicalDistance )
            {
                rsmdWARNING( std::setprecision(3) << "... atom " << atomAfter->name << " " << atomAfter->id << " of molecule " << newMolecule.getName() << " " << newMolecule.getID() 
                        << " moved more than twice the typical distance: " << distance << ' ' << unitSystem->length << " ( > 2 * " << typicalDistance << ' ' << unitSystem->length << ")");
            }
            else
            {
                rsmdDEBUG( "... atom " << atomAfter->name << " " << atomAfter->id << " of molecule " << newMolecule.getName() << " " << newMolecule.getID() 
                        << " moved: " << distance << ' ' << unitSystem->length);
            }
            ++ atomBefore;
            ++ atomAfter;
        }
    }
}


//
// check if a candidate is still available
//
bool Universe::isAvailable( const ReactionCandidate& candidate )
{
    bool reactantsAreAvailable = true;
    for( const auto& reactant: candidate.getReactants() )
    {
        if( ! topologyNew.containsMolecule(reactant) )
        {
            rsmdDEBUG( "couldn't find molecule " << reactant.getName() << " " << reactant.getID() << " in topology" );
            reactantsAreAvailable = false;
            break;
        }
    }
    return reactantsAreAvailable;
}

void Universe::makeMoleculeWhole(Molecule& molecule, const REALVEC& dimensions)
{
    rsmdLOG( "... repairing molecule in case it is broken across periodic boundaries: " << molecule );
    Atom& referenceAtom = molecule.front();
    for(auto& atom: molecule)
    {   
        REALVEC before = atom.position;
        REALVEC distance = (atom.position - referenceAtom.position);
        bool moved = false;
        for( std::size_t i=0; i<3; ++i )
        {
            atom.position[i] -= static_cast<int>( distance[i] / (0.5 * dimensions[i]) ) * dimensions[i];
            if( static_cast<int>(distance[i] / (0.5*dimensions[i]) ) != 0 )    moved = true;
        }
        REALVEC after = atom.position;
        if( moved )
        {
            rsmdLOG( "    before: " << before );
            rsmdLOG( "    after: " << after );
        }
    }
}

//
// react a given candidate
// (checks for whether the molecules are still available need to happen before!)
//
void Universe::react(ReactionCandidate& candidate)
{
    rsmdDEBUG( "performing reaction for candidate " << candidate.shortInfo() );
   
    // reactant --> product translation 
    candidate.applyTransitions();
    // make products whole
    for(auto& product: candidate.getProducts())
    {
        makeMoleculeWhole(product, topologyNew.getDimensions());
    }
    // apply translational movements of product atoms
    candidate.applyTranslations();

    // apply changes to topology
    auto highestMolID = std::max_element( std::begin(topologyNew), std::end(topologyNew), [](const auto& mol1, const auto& mol2){ return mol1.getID() < mol2.getID(); } )->getID(); 
    for( const auto& reactant: candidate.getReactants() )
    {
        topologyNew.removeMolecule( reactant.getID() );    
    }
    for( auto& product: candidate.getProducts() )
    {
        product.setID( ++highestMolID );
        auto molecule __attribute__((unused)) = topologyNew.addMolecule( product );
        topologyNew.addReactionRecord( highestMolID );
        // topologyNew.repairMoleculePBC( *molecule );
        rsmdDEBUG( "new molecule " << molecule->getName() << " got ID " << molecule->getID() );
    }
}


//
// search for reaction candidates
//
std::vector<ReactionCandidate> Universe::searchReactionCandidates()
{
    // search for possible reaction candidates and return them if they match all criteria
    std::vector<ReactionCandidate> reactionCandidates {};

    for( auto& reactionTemplate: reactionTemplates )
    {
        if( reactionTemplate.getReactants().size() == 1 )
        {
            for( auto& reactant: topologyOld.getMolecules( reactionTemplate.getReactants()[0].getName() ) )
            {
                rsmdDEBUG( "checking reaction candidate: " << reactant.get().getName() << ", " << reactant.get().getID() );
                reactionCandidates.push_back( reactionTemplate );
                reactionCandidates.back().updateReactant( 0, reactant.get() );
                if( ! reactionCandidates.back().valid(topologyOld.getDimensions()) ) reactionCandidates.pop_back();
            }
        }
        else if( reactionTemplate.getReactants().size() == 2 )
        {
            for( auto& reactant1: topologyOld.getMolecules( reactionTemplate.getReactants()[0].getName() ) )
            {
                for( auto& reactant2: topologyOld.getMolecules( reactionTemplate.getReactants()[1].getName() ) )
                {
                    if( reactant1.get() == reactant2.get() ) continue;
                    if( reactant1.get().getName() == reactant2.get().getName() && reactant1.get().getID() > reactant2.get().getID() ) continue;
                    rsmdDEBUG( "checking reaction candidate: " << reactant1.get().getName() << ", " << reactant1.get().getID() << " + " << reactant2.get().getName() << ", " << reactant2.get().getID() );
                    reactionCandidates.push_back( reactionTemplate );
                    reactionCandidates.back().updateReactant( 0, reactant1.get() );
                    reactionCandidates.back().updateReactant( 1, reactant2.get() );
                    if( ! reactionCandidates.back().valid(topologyOld.getDimensions()) )   reactionCandidates.pop_back();

                }
            }
        }
        else if( reactionTemplate.getReactants().size() == 3 )
        {
            for( auto& reactant1: topologyOld.getMolecules( reactionTemplate.getReactants()[0].getName() ) )
            {
                for( auto& reactant2: topologyOld.getMolecules( reactionTemplate.getReactants()[1].getName() ) )
                {
                    if( reactant1.get() == reactant2.get() ) continue;
                    if( reactant1.get().getName() == reactant2.get().getName() && reactant1.get().getID() > reactant2.get().getID() ) continue; 
                    for( auto& reactant3: topologyOld.getMolecules( reactionTemplate.getReactants()[2].getName() ) )
                    {
                        if( reactant1.get() == reactant3.get() || reactant2.get() == reactant3.get() )  continue;
                        if( reactant2.get().getName() == reactant3.get().getName() && reactant2.get().getID() > reactant3.get().getID() ) continue;
                        rsmdDEBUG( "checking reaction candidate: " << reactant1.get().getName() << ", " << reactant1.get().getID() 
                                                          << " + " << reactant2.get().getName() << ", " << reactant2.get().getID() 
                                                          << " + " << reactant3.get().getName() << ", " << reactant3.get().getID() );
                        reactionCandidates.push_back( reactionTemplate );
                        reactionCandidates.back().updateReactant( 0, reactant1.get() );
                        reactionCandidates.back().updateReactant( 1, reactant2.get() );
                        reactionCandidates.back().updateReactant( 2, reactant3.get() );
                        if( ! reactionCandidates.back().valid(topologyOld.getDimensions()) )   reactionCandidates.pop_back();
                    }
                }
            }
        }
        else
        {
            rsmdCRITICAL("attention: more than 3 reactants per reaction is currently not implemented!");
        }
        
    }

    // shuffle candidates
    enhance::shuffle(reactionCandidates.begin(), reactionCandidates.end());

    return reactionCandidates;
}

