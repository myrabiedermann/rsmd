
#include "parser/topologyParserGMX.hpp"



void TopologyParserGMX::read( Topology& topology, const std::size_t& cycle )
{
    // convert filenames
    std::stringstream topFile {};
    std::stringstream coordFile {};
    topFile << cycle << ".top";
    coordFile << cycle << "-md.gro";

    // read topology
    auto topologyMap = read_top( topFile.str() );
    read_gro( coordFile.str(), topology );

    // some consistency checks:
    unsigned int atomCounter = 0;
    for( const auto& moleculetype: topologyMap )
    {
        auto foundMolecules = topology.getMolecules( moleculetype.first );
        if( foundMolecules.size() != moleculetype.second )   
            rsmdWARNING(".top and .gro don't match (# molecules of type " << moleculetype.first << " " << moleculetype.second << " vs. " << foundMolecules.size() << ")") 
        atomCounter += foundMolecules.size();
    }
    if( atomCounter != topology.size() )
        rsmdWARNING( " total number of molecules in .gro and .top doesn't match" << "(" << atomCounter << " vs. " << topology.size() << ")" )
}

void TopologyParserGMX::readRelaxed( Topology& topology, const std::size_t& cycle )
{
    // convert filenames
    std::stringstream topFile {};
    std::stringstream coordFile {};
    topFile << cycle << ".top";
    coordFile << cycle << "-rs.gro";

    // read topology
    auto topologyMap = read_top( topFile.str() );
    read_gro( coordFile.str(), topology );

    // some consistency checks:
    unsigned int atomCounter = 0;
    for( const auto& moleculetype: topologyMap )
    {
        auto foundMolecules = topology.getMolecules( moleculetype.first );
        if( foundMolecules.size() != moleculetype.second )   
            rsmdWARNING(".top and .gro don't match (# molecules of type " << moleculetype.first << " " << moleculetype.second << " vs. " << foundMolecules.size() << ")") 
        atomCounter += foundMolecules.size();
    }
    if( atomCounter != topology.size() )
        rsmdWARNING( " total number of molecules in .gro and .top doesn't match" << "(" << atomCounter << " vs. " << topology.size() << ")" )
}


void TopologyParserGMX::write(Topology& top, const std::size_t& currentCycle)
{
    rsmdDEBUG(__PRETTY_FUNCTION__);

    // convert filenames
    std::stringstream cycle {};
    cycle << currentCycle;

    // write topology
    write_top( cycle.str() + ".top", top );
    write_gro( cycle.str() + "-rs.gro", top );
    write_index( cycle.str() + ".reactants.ndx", cycle.str() + ".products.ndx", top );
}



std::map<std::string, unsigned int> TopologyParserGMX::read_top( const std::string& topFile )
{
    std::map<std::string, unsigned int> topologyMap {};

    bool readFileContent = (topologyFileContent.empty() ? true : false);

	// open topology file
	std::ifstream FILE( topFile );
	if( ! FILE ){   // safety check
		rsmdCRITICAL( topFile << " doesn't exist, cannot read topology" )
	} else {
        bool directiveMolecules {false};
        bool directiveSystem    {false};
        while( FILE.good() )
        {
            // get next line
            std::string line;
            std::getline(FILE, line, '\n');
            std::stringstream linestream(line);

            // case: [system] --> read and save systemName
            if( directiveSystem )
            {
                if( line.empty() ) continue;
                systemName = enhance::trimString(line);
                directiveSystem = false;
            }
            // case: [molecules] --> read and save 
            else if( directiveMolecules )
            {
                if( line.empty() ) continue;
                std::string moltype {};
                int nMolecules = 0;
                linestream >> moltype >> nMolecules;
                topologyMap[moltype] = nMolecules;
            }
            // case else --> check if line containes a directive
            else
            {   
                // any kind of directive []
                if( line.find('[') != std::string::npos )
                {
                    // seach for directive [system]
                    auto pos = line.find("system");
                    if( pos != std::string::npos )
                    {   
                        directiveSystem = true;
                        directiveMolecules = false;
                    }
                    // seach for directive [molecules]
                    pos = line.find("molecules");
                    if( pos != std::string::npos )
                    {   
                        directiveSystem = false;
                        directiveMolecules = true;
                    }
                }
                // read and save line to topologyFileContent
                if( readFileContent )
                {
                    topologyFileContent.emplace_back(line);
                }
            }
        }
    }   
	FILE.close();

    return topologyMap;
}



void TopologyParserGMX::read_gro( const std::string& groFile, Topology& top )
{	
	int totNrOfAtoms = 0;
	
	// open .gro file
	std::ifstream FILE( groFile );
	
	if( ! FILE ){   // check if file exists
       rsmdCRITICAL(groFile << " doesn't exist, cannot read structure")
	} else {			
		
		// first line: system name
		std::string line;
		std::getline(FILE, line, '\n');
        line = enhance::trimString(line);
        if( line != systemName )    rsmdWARNING("system names don't agree (" << systemName << " vs. " << line << ")")

        // second line: number of atoms
        std::getline(FILE, line, '\n');
        std::stringstream linestream(line);
        linestream >> totNrOfAtoms;

        // read atom descriptions
        int counter = 0;
        while( counter < totNrOfAtoms )
        {
            std::getline(FILE, line, '\n');
            if(line.size() == 44)
            {
                line.append(std::string("  0.0000"));   // append empty velocities
                line.append(std::string("  0.0000"));
                line.append(std::string("  0.0000"));
            }
           
            // molecule related information
            int resid           = std::stoi( line.substr(0,5) );
            std::string resname = line.substr(5,5);
            resname.erase(std::remove_if( resname.begin(), resname.end(), ::isspace), resname.end());
           
            // atom related information
            Atom atom;
            atom.name = line.substr(10,5);
            atom.name.erase(std::remove_if( atom.name.begin(), atom.name.end(), ::isspace), atom.name.end());
            atom.id = std::stoi( line.substr(15,5) );
            atom.position(0) = std::stof( line.substr(20,8) );
            atom.position(1) = std::stof( line.substr(28,8) );
            atom.position(2) = std::stof( line.substr(36,8) );
            atom.velocity(0)= std::stof( line.substr(44,8) );
            atom.velocity(1) = std::stof( line.substr(52,8) );
            atom.velocity(2) = std::stof( line.substr(60,8) );

            // add atom and all infos to topology:
            auto& mol = top.getAddMolecule( resid, resname );
            mol.addAtom(atom);

            counter ++;
        }

        // last line: box vector
        std::getline(FILE, line, '\n');
        std::stringstream tmpstream(line);
        REALVEC box;
        tmpstream >> box(0) >> box(1) >> box(2);
        top.setDimensions(box);

	}
	
	// close file
	FILE.close();
	
}



void TopologyParserGMX::write_top( const std::string& topFile, Topology& top )
{
    std::ofstream FILE( topFile );
    if( FILE.bad() ) rsmdCRITICAL("something went wrong with outstream to " << topFile);
    
    for(auto& line: topologyFileContent)
    {
        // if line contains [system]
        if( line.find('[') != std::string::npos && line.find("system") != std::string::npos )
        {
            FILE << line << '\n';
            FILE << systemName << '\n';
        }
        // if line contains [molecules]
        else if( line.find('[') != std::string::npos && line.find("molecules") != std::string::npos )
        {
            FILE << line << '\n';
            for(auto& mt: top.getMoleculetypes() )
            {
                auto number = top.getMolecules( mt ).size();       
                FILE << std::setw(5) << std::left << mt << number << '\n';
            }
        }  
        else
        {
            FILE << line << '\n';
        }
    }
    
    FILE.close();
}


void TopologyParserGMX::write_gro( const std::string& groFile, Topology& top )
{

    std::ofstream FILE( groFile );
    if( FILE.bad() ) rsmdCRITICAL("something went wrong with outstream to " << groFile);

    // first two lines: system name / other info and # of atoms
    FILE << systemName << " (created by reactiveMD)" << '\n';
    FILE << std::setw(6) << top.getNAtoms() << '\n';
    
    // assumes that topology has been sorted beforehand
    // (gromacs needs molecules sorted according to types and this has to match the sequence in .top !)
    for( const auto& mol: top )
    {
        for(const auto& atom: mol)
        {
            FILE << std::setw(5) << std::right << mol.getID() 
                 << std::setw(5) << std::left  << mol.getName()
                 << std::setw(5) << std::right << atom.name
                 << std::setw(5) << std::right << atom.id;
            for( const auto& p: atom.position )
                FILE << std::fixed << std::right << std::setprecision(3) << std::setw(8) << p;
            for( const auto& v: atom.velocity )
                FILE << std::fixed << std::right << std::setprecision(4) << std::setw(8) << v;
            FILE << '\n';
        }
    }

    // box dimensions
    for( const auto& d: top.getDimensions() )
        FILE << std::setw(10) << std::setprecision(6) << d;
    FILE << '\n';
    
    FILE.close();
    
}



void TopologyParserGMX::write_index(const std::string& reactants, const std::string& products, Topology& top) 
{
    std::ofstream REACTANTS( reactants );
    if( REACTANTS.bad() ) rsmdCRITICAL("something went wrong with outstream to 'reactants.ndx'");
    std::ofstream PRODUCTS( products );
    if( PRODUCTS.bad() ) rsmdCRITICAL("something went wrong with outstream to 'products.ndx'");

    REACTANTS << "[xxx]\n";
    PRODUCTS << "[xxx]\n";

    for( const auto& idpair: top.getReactionRecordsAtoms() )
    {
        REACTANTS << idpair.first << " ";
        PRODUCTS << idpair.second << " ";
    }
    REACTANTS << '\n';
    PRODUCTS << '\n';

    REACTANTS.close();
    PRODUCTS.close();
}