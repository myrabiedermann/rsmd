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

#include "parser/reactionParser.hpp"

//
// read reaction from file
//
ReactionBase ReactionParser::read(const std::string& reactionFile)
{
    ReactionBase reaction {};
    std::vector<std::pair<REAL, REAL>> reactionRate {};

    std::ifstream FILE (reactionFile);
	if( ! FILE )
    {   // safety check
		rsmdCRITICAL( reactionFile << " doesn't exist, cannot read reaction")
	} 
    else 
    {
        std::string currentDirective {};
        while( FILE.good() )
        {
            // get next line
            std::string line;
            std::getline(FILE, line, '\n');
            std::stringstream linestream(line);

            // interpret line and save accordingly
            if( line.empty() )
            {
                // empty line --> skip
                continue;
            }
			else if( line.find_first_not_of(" \t\n\v\f\r") == std::string::npos )
			{
				// line contains only spaces --> skip
				continue;
			}
            else if( line[line.find_first_not_of(' ')] == '#' )
            {
                // comment line --> skip
                continue;
            }
            else if( line.find('[') != std::string::npos )
            {
                // begin of new directive
                auto pos = line.find("[");
                auto pos2 = line.find("]");
                if(pos2 == std::string::npos) rsmdCRITICAL("something is wrong with file formatting, couldn't find enclosing ]");
                currentDirective = line.substr(pos+1, (pos2-pos-1));
                currentDirective = enhance::trimString(currentDirective);
            } 
            else 
            {
                // first: check if part of line is a comment
                auto pos = line.find("#");
                if( pos != std::string::npos )
                {
                    line.erase(pos, std::string::npos);
                }
                // second: parse content of line
                if( currentDirective == "name" )
                {
                    // remove white spaces before/after name
                    line = enhance::trimString(line);
                    reaction.setName(line);
                }
                else if( currentDirective == "reactants" )
                {
                    // # molNr   molName     atomName  atomNr
                    int molid = 0, atomid = 0;
                    std::string molname {},  atomname {};
                    linestream >> molid >> molname >> atomname >> atomid;

                    Molecule& mol = reaction.getAddReactant(molid);
                    mol.setName(molname);
                    mol.addAtom(atomid, atomname);
                }
                else if( currentDirective == "products" )
                {
                    // # molNr      molName     atomName    atomNr  origin->molNr   origin->atomNr
                    int molNr = 0, atomNr = 0, oldmolNr= 0, oldatomNr = 0;
                    std::string molname {}, atomname {};
                    linestream >> molNr >> molname >> atomname >> atomNr >> oldmolNr >> oldatomNr;
                    
                    Molecule& mol = reaction.getAddProduct(molNr);
                    mol.setName(molname);
                    mol.addAtom(atomNr, atomname);
                    // remove -1 in order to transform atomNr/molNr in indices
                    reaction.addTransition(oldmolNr-1, oldatomNr-1, molNr-1, atomNr-1);
                }
                else if( currentDirective == "criteria" )
                {
                    std::string type {};
                    std::vector<std::pair<std::size_t, std::size_t>> atomIDs {};
                    std::pair<std::size_t, std::size_t> thisID {};
                    REAL minValue = 0;
                    REAL maxValue = 0;
                    linestream >> type;
                    
                    if( type == "DIST" or type == "dist" or type == "Dist" )
                    {
                        for(int i=0; i<2; ++i)
                        {
                            linestream >> thisID.first >> thisID.second;
                            atomIDs.push_back(thisID);
                        }
                    }
                    else if( type == "ANG" or type == "ang" or type == "Ang" )
                    {
                        for(int i=0; i<3; ++i)
                        {
                            linestream >> thisID.first >> thisID.second;
                            atomIDs.push_back(thisID);
                        }
                    }
                    else if( type == "DIH" or type == "dih" or type == "Dih" )
                    {
                        for(int i=0; i<4; ++i)
                        {
                            linestream >> thisID.first >> thisID.second;
                            atomIDs.push_back(thisID);
                        }
                    }
                    linestream >> minValue >> maxValue;
                    
                    // transform molID/atomID to an index by removing -1
                    for(auto& ix: atomIDs) 
                    {
                        ix.first -= 1; 
                        ix.second -= 1;
                    }

                    reaction.addCriterion(atomIDs, std::make_pair(minValue, maxValue));
                }
                else if( currentDirective == "translations" )
                {
                    std::vector<std::pair<std::size_t, std::size_t>> atomIDs {};
                    REAL value = 0;
                    
                    for(int i=0; i<2; ++i)
                    {
                        std::pair<std::size_t, std::size_t> thisID {};
                        linestream >> thisID.first >> thisID.second;
                        atomIDs.push_back(thisID);
                    }
                    linestream >> value;
                   
                    // transform molID/atomID to an index by removing -1
                    for(auto& ix: atomIDs) 
                    {
                        ix.first -= 1; 
                        ix.second -= 1;
                    }

                    reaction.addTranslation(atomIDs, value);
                }
                else if( currentDirective == "energy" )
                {
                    REAL value = 0;
                    linestream >> value;
                    
                    reaction.setReactionEnergy(value);
                }
                else if( currentDirective == "activation" )
                {
                    REAL value = 0;
                    linestream >> value;
                    
                    reaction.setActivationEnergy(value);
                }
                else if( currentDirective == "rate" )
                {
                    // # distance	rate 
                    std::pair<REAL, REAL> values {};
                    linestream >> values.first >> values.second;
                    
                    reactionRate.emplace_back(values);
                }
            }

        } // end of while loop through file
        // finally: sort and add rate values:
        std::stable_sort( std::begin(reactionRate), std::end(reactionRate), [](const auto& lhs, const auto& rhs){ return lhs.first < rhs.first; });
        reaction.setRate(reactionRate);
    }
    return reaction;
}



//
// write example file
//
std::string ReactionParser::writeExample()
{
    std::stringstream stream;
    
    // directive name
    stream << "[name]" << '\n'
           << "example reaction\n"
           << '\n';
    
    // directive reactants
    stream << "[reactants]" << '\n'
           << "# mandatory: you have to specify at least one reacting molecule!\n"
           << "# molID      molName     atomName    atomID\n"
           << "  1          MOL         CM          1\n"
           << "  1          MOL         HM          2\n"
           << "  1          MOL         HM          3\n"
           << "  1          MOL         HM          4\n"
           << "\n"
           << "  2          MOL         CM          1\n"
           << "  2          MOL         HM          2\n"
           << "  2          MOL         HM          3\n"
           << "  2          MOL         HM          4\n"
           << '\n';
    
    // directive products
    stream << "[products]" << '\n'
           << "# mandatory: you have to specify at least one product molecule!\n"
           << "# molID      molName     atomName   atomID     origin->molID      origin->atomID\n"
           << "  1          NEW         CE          1         1                  1\n"
           << "  1          NEW         HE          2         1                  2\n"
           << "  1          NEW         HE          3         1                  3\n"
           << "  1          NEW         HE          4         1                  4\n"
           << "  1          NEW         CE          5         2                  1\n"
           << "  1          NEW         HE          6         2                  2\n"
           << "  1          NEW         HE          7         2                  3\n"
           << "  1          NEW         HE          8         2                  4\n"
           << '\n';

    // directive criteria
    stream << "[criteria]" << '\n'
           << "# mandatory: you have to specify at least one distance criterion!\n"
           << "# one criterion per line, beginning with a type specification (dist / ang / dih)\n"
           << "# type   molID   atomID   molID   atomID   minValue   maxValue\n" 
           << "  dist   1       1        2       1        0.0        4.0\n"
           << "# type   molID   atomID   molID   atomID   molID   atomID   minValue   maxValue\n" 
           << "  ang    1       2        1       1        2       1        110        150\n"
           << "# type   molID   atomID   molID   atomID   molID   atomID   molID   atomID   minValue   maxValue\n" 
           << "  dih    1       2        1       1        2       1        2       2        -20        20\n"
           << '\n';
    
    // directive translations
    stream << "[translations]" << '\n'
           << "# optional: you can define translations for specific atoms along an atom-atom connection axis\n"
           << "#           which will be performed during the transformation reactants -> products.\n"
           << "# one translation per line. a positive value moves the first atom towards the second atom, negative values move away.\n"
           << "# molID  atomID  molID   atomID  value\n"
           << "  1      1       2       2       1.0\n"
           << '\n';

    // directive energy
    stream << "[energy]" << '\n'
           << "# optional: only required if a Metropolis MC acceptance criterion should be used!\n"
           << "# correction term for reaction energy of this reaction.\n"
           << "# value\n"
           << "  -66.0 \n"
           << '\n';
   
    // directive activation energy
    stream << "[activation]" << '\n'
           << "# optional: only used if a Metropolis MC acceptance criterion should be used!\n"
           << "# energy of activation for this reaction."
           << "# value\n"
           << "  10.0 \n"
           << '\n';

    // directive rate
    stream << "[rate]" << '\n'
           << "# optional: only required if a rate-based acceptance criterion should be used!\n"
           << "# attention: given distances are assumed to correspond to the first distance criterion that is given!\n"
           << "# distance   rate value\n"
           << "  0.35       0.2\n"
           << "  0.40       0.06\n"
           << "  0.50       0.04\n"
           << '\n';

    stream << "# notes: \n";
    stream << "# - the # symbol marks the beginning of a comment \n";
    stream << "# - comments or empty lines are ignored\n";
    stream << "# - white spaces suffice to separate columns\n";
    stream << "# - all units are set to same units as in the corresponding md engine\n";
    stream << "# - atoms in directives [reactants] and [products] must appear consecutively numbered, \n"
           << "#   in the same order as in the corresponding topology files\n"
           << "#   and named accordingly\n";

    return stream.str();
}
