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

#include "engine/engineGMX.hpp"


void EngineGMX::setup(const Parameters& parameters)
{
    // set parameters:
    executablePath = parameters.getOption("simulation.engine").as<std::string>();

    mdp_file =            parameters.getOption("gromacs.mdp").as<std::string>();
    mdp_file_energy =     parameters.getOption("gromacs.mdp.energy").as<std::string>();
    mdp_file_relaxation = parameters.getOption("gromacs.mdp.relaxation").as<std::string>();
    if( parameters.getOption("reaction.mc").as<bool>() ) 
    {
        // set description for energy computation
        if( parameters.getOption("reaction.computeLocalPotentialEnergy").as<bool>() )
        {
            computeLocalPotentialEnergies = true;
            if( parameters.getOption("reaction.computeSolvationPotentialEnergy").as<bool>() )
            {
                computeSolvationPotentialEnergies = true;
            }
        }
        if( parameters.getOption("reaction.averagePotentialEnergy").as<REAL>() != 0 )
        {
            averagePotentialEnergies = true;
        }
    }

    // set extension time for appending simulations
    read_mdp( mdp_file );
    extensionTime_str = std::to_string(extensionTime);

    // get usable number of threads:
    int nt = parameters.getOption("gromacs.nt").as<int>();
    int ntmpi = parameters.getOption("gromacs.ntmpi").as<int>();
    int ntomp = parameters.getOption("gromacs.ntomp").as<int>();
    if( nt == 0 && ntmpi == 0 && ntomp == 0 )
    {
        rsmdLOG( "gromacs.nt, gromacs.tmpi and gromacs.ntomp are all set to zero." );
        nt = std::thread::hardware_concurrency();
        rsmdLOG( "... detected " << nt << " threads on this machine, setting gromacs.nt to " << nt );
    }
    nt_as_str = std::to_string(nt);
    ntmpi_as_str = std::to_string(ntmpi);
    ntomp_as_str = std::to_string(ntomp);

    // check what to do in cleanup() after rs was rejected:
    saveRejectedFiles = parameters.getOption("reaction.saveRejected").as<bool>();
    rejectedFilekeys = {".top", "-rs.tpr", "-rs.gro", "-rs.log", "-rs.edr", "-rs.cpt", "-rs.xtc", "-rs-mdpout.mdp", ".reactants.ndx", ".products.ndx"};
    if( parameters.getOption("reaction.mc").as<bool>() )    rejectedFilekeys.emplace_back("-rs.xvg");

    // set backup policy
    if( parameters.getOption("gromacs.backup").as<bool>() )
    {
        backupPolicy = "-backup";
    }

    // check that topology/coordinate files are present
    std::string topologyFile = parameters.getOption("gromacs.topology").as<std::string>();
    std::string coordinatesFile = parameters.getOption("gromacs.coordinates").as<std::string>();
    switch( parameters.getSimulationMode() )
    {
        case SIMMODE::NEW:
            if( ! std::filesystem::exists(std::filesystem::current_path()/"0.top") )
            {
                rsmdLOG( "... copying '" << topologyFile << "' -> '0.top'" );
                std::filesystem::copy_file(std::filesystem::current_path()/topologyFile, std::filesystem::current_path()/"0.top");
            }
            if( ! std::filesystem::exists(std::filesystem::current_path()/"0-md.gro") )
            {
                rsmdLOG( "... copying '" << coordinatesFile << "' -> '0-md.gro'" );
                std::filesystem::copy_file(std::filesystem::current_path()/coordinatesFile, std::filesystem::current_path()/"0-md.gro");
            }
            break;
        
        case SIMMODE::RESTART:
            if( ! std::filesystem::exists( std::filesystem::current_path()/topologyFile ) )
            {
                rsmdCRITICAL("existence of topology file '" << topologyFile << "' is mandatory in order to restart the simulation" );
            }
            if( ! std::filesystem::exists( std::filesystem::current_path()/coordinatesFile ) )
            {
                rsmdCRITICAL("existence of coordinates file '" << coordinatesFile << "' is mandatory in order to restart the simulation" );
            }
            break;
    }


    // check given executable:
    verifyExecutable();
}



void EngineGMX::verifyExecutable() 
{
    rsmdLOG( "... checking simulation.engine ..." );
    bool okay {true};
    
    try
    {
        execute(executablePath.c_str(), executablePath.c_str(), "-version", "-quiet", "-nocopyright");
    }
    catch(const std::exception& e)
    {
       rsmdCRITICAL( e.what() );
       okay = false;
    }
    
    if( okay )
    {
        rsmdLOG( "... simulation.engine seems to be okay" );
    }
    else
    {
        rsmdCRITICAL(  "... simulation.engine failed!" );
    }
    
}



// md           in: cycle = X
//              grompp -f normal.mdp -t X-rs.cpt -c X-rs.gro -p X-rs.top -o X-md.tpr
//              mdrun  -s X-md.tpr -deffnm X-md
void EngineGMX::runMD( const std::size_t& cycle)
{
    std::stringstream keyIn, keyOut, key {};
    keyIn << cycle << "-rs";
    keyOut << cycle << "-md";
    key << cycle;
    try
    {
        // run grompp -f mdp.mdp -p top -c gro.gro -o tpr.tpr -t cpt.cpt
        // void EngineGMX::grompp( const std::string& mdp, const std::string& top, const std::string& gro, const std::string& tpr )
        grompp( mdp_file, key.str(), keyIn.str(), keyOut.str());

        // run mdrun -s tpr.tpr -deffnm tpr
        // void EngineGMX::mdrun( const std::string& tpr )
        mdrun( keyOut.str() );
    }
    catch(const std::exception& e)
    {
        rsmdCRITICAL( "catched expection in EngineGMX::runMD(): " << e.what() );
    }
    
}



// md           in: cycle = 0
//              grompp -f normal.mdp -c 0-md.gro -p 0.top -o 0-md.tpr
//              mdrun  -s 0-md.tpr -deffnm 0-md
void EngineGMX::runMDInitial()
{
    try
    {
        grompp( mdp_file, "0", "0-md", "0-md");
        mdrun( "0-md" );
    }
    catch(const std::exception& e)
    {
        rsmdCRITICAL( "caught expection in EngineGMX::runMDInitial(): " << e.what() );
    }
}


// mdAppending  in: cycle = X, lastReactiveCycle = Y, time = T
//              convert-tpr -s (X-1)-md.tpr -f Y-md.cpt -o X-md.tpr -extend T
//              mdrun  -s X-md.tpr -cpi Y-md.cpt -append -deffnm Y-md
void EngineGMX::runMDAppending( const std::size_t& cycle, const std::size_t& lastReactiveCycle)
{
    std::stringstream tprOld, tpr, key {};
    tprOld << (cycle - 1)  << "-md";
    tpr << cycle << "-md";
    key << lastReactiveCycle << "-md";

    try
    {
        // run convert-tpr -s tpr.tpr -f cpt.cpt -o tpr_new.tpr -extend time
        // void EngineGMX::convert_tpr( const std::string& tpr, const std::string& tpr_new)
        convert_tpr( tprOld.str(), tpr.str()); 
        
        // run mdrun -s tpr.tpr -deffnm fnm -cpi cpt.cpt -append
        // void EngineGMX::mdrun( const std::string& tpr, const std::string& fnm, const std::string& cpt )
        mdrun( tpr.str(), key.str(), key.str() );
    }
    catch(const std::exception& e)
    {
        rsmdCRITICAL( "caught expection in EngineGMX::runMDAppending(): " << e.what() );
    }
    
}



// rs / relax   in: cycle = X 
//              grompp -f relax.mdp -c X-rs.gro -p X-rs.top -o X-rs.tpr
//              mdrun  -s X-rs.tpr -deffnm X-rs
bool EngineGMX::runRelaxation( const std::size_t& cycle )
{
    std::stringstream keyOut, key {};
    keyOut << cycle << "-rs";
    key << cycle;
    bool statusRelaxation = true;

    try
    {
        // run grompp -f mdp.mdp -p top -c gro.gro -o tpr.tpr
        // void EngineGMX::grompp( const std::string& mdp, const std::string& top, const std::string& gro, const std::string& tpr )
        grompp( mdp_file_relaxation, key.str(), keyOut.str(), keyOut.str() );

        // run mdrun -s tpr.tpr -deffnm tpr
        // void EngineGMX::mdrun( const std::string& tpr )
        mdrun( keyOut.str() );
    }
    catch(const std::exception& e)
    {
        rsmdWARNING( "caught expection in EngineGMX::runRelaxation(): " << e.what() );
        statusRelaxation = false;
    }
    return statusRelaxation;
}



// energy   in: cycle = X, lastReactiveCycle = Y 
//          energy -f X-rs.edr -o X-rs.xvg
//          energy -f Y-md.edr -o Y-md.xvg
void EngineGMX::runEnergyComputation( const std::size_t& currentCycle, const std::size_t& lastReactiveCycle )
{
    //      energy -f edr.edr -o xvg.xvg
    std::stringstream before, after, cycle, cycleBefore {};
    before << lastReactiveCycle << "-md"; 
    after << currentCycle << "-rs";
    cycle << currentCycle;
    cycleBefore << lastReactiveCycle;

    try
    {
        if( computeLocalPotentialEnergies )
        {
            std::string backup = backupPolicy;
            backupPolicy = "-nobackup";

            if( computeSolvationPotentialEnergies )
            {
                // create tpr files
                // ... for only reactant/product atoms
                convert_tpr( before.str(), "reactants", cycle.str()+".reactants" );
                convert_tpr( after.str(), "products", cycle.str()+".products" );
                // ... for solvation group
                grompp(mdp_file_energy, cycleBefore.str(), before.str(), "reactants_solvation", cycle.str()+".reactants");
                grompp(mdp_file_energy, cycle.str(), after.str(), "products_solvation", cycle.str()+".products");

                if( averagePotentialEnergies )
                {
                    // create .gro file for only reactant/product atoms
                    trjconv( before.str(), cycle.str()+".reactants", before.str()+".xtc", "reactants.xtc" );
                    trjconv( after.str(), cycle.str()+".products", after.str()+".xtc", "products.xtc" );
                    // ... and mdrun rerun to create .edr file
                    mdrunRerun("reactants", "reactants.xtc", "reactants");
                    mdrunRerun("products", "products.xtc", "products");
                    // mdrun rerun for solvation group
                    mdrunRerun("reactants_solvation", before.str()+".xtc", "reactants_solvation");
                    mdrunRerun("products_solvation", after.str()+".xtc", "products_solvation");
                }
                else
                {
                    // create .gro file for only reactant/product atoms
                    trjconv( before.str(), cycle.str()+".reactants", before.str()+".gro", "reactants.gro" );
                    trjconv( after.str(), cycle.str()+".products", after.str()+".gro", "products.gro" );
                    // ... and mdrun rerun to create .edr file
                    mdrunRerun("reactants", "reactants.gro", "reactants");
                    mdrunRerun("products", "products.gro", "products");
                    // mdrun rerun for solvation group
                    mdrunRerun("reactants_solvation", before.str()+".gro", "reactants_solvation");
                    mdrunRerun("products_solvation", after.str()+".gro", "products_solvation");
                }

                // convert to .xvg files
                energy( "reactants", before.str() );
                energy( "products", after.str() );
                energySolvation( "reactants_solvation", "reactants_solvation" );
                energySolvation( "products_solvation", "products_solvation" );
            }
            else
            {
                // first: create .tpr file for only reactant/product atoms
                convert_tpr( before.str(), "reactants", cycle.str()+".reactants" );
                convert_tpr( after.str(), "products", cycle.str()+".products" );
                if( averagePotentialEnergies )
                {
                    // second: create .gro file or .xtc file for only reactant/product atoms
                    trjconv( before.str(), cycle.str()+".reactants", before.str()+".xtc", "reactants.xtc" );
                    trjconv( after.str(), cycle.str()+".products", after.str()+".xtc", "products.xtc" );
                    // third: mdrun rerun to create .edr file
                    mdrunRerun("reactants", "reactants.xtc", "reactants");
                    mdrunRerun("products", "products.xtc", "products");
                }
                else
                {
                    // second: create .xtc file for only reactant/product atoms
                    trjconv( before.str(), cycle.str()+".reactants", before.str()+".gro", "reactants.gro" );
                    trjconv( after.str(), cycle.str()+".products", after.str()+".gro", "products.gro" );
                    // third: mdrun rerun .xtc to create .edr file
                    mdrunRerun("reactants", "reactants.gro", "reactants");
                    mdrunRerun("products", "products.gro", "products");
                }
                // forth: convert to .xvg file
                energy( "reactants", before.str() );
                energy( "products", after.str() );
            }
            backupPolicy = backup;
        }
        else
        {
            energy( before.str(), before.str() );
            energy( after.str(), after.str() );
        }
    }
    catch(const std::exception& e)
    {
        rsmdCRITICAL( "caught expection in EngineGMX::runEnergyComputation(): " << e.what() );
    }
}


//
// cleanup: rename or delete all files produced during the rejected reactive step
//
void EngineGMX::cleanup( const std::size_t& cycle )
{
    std::string key = std::to_string(cycle);
    std::filesystem::path thisPath = std::filesystem::current_path();

    if( saveRejectedFiles )
    {
        rsmdDEBUG("... moving files from rejected reactive step");
        for( auto filename: rejectedFilekeys )
        {
            try
            {
                std::filesystem::rename( thisPath/(key+filename), thisPath/("rejected-"+key+filename) ); 
            }
            catch(const std::exception& e)
            {
                rsmdWARNING( "   caught exception while trying to rename " << thisPath/(key+filename) << ": " << e.what() );
            }
        }
    }
    else
    {
        rsmdDEBUG("... deleting files from rejected reactive step");
        for( auto filename: rejectedFilekeys )
        {
            try
            {
                std::filesystem::remove(thisPath/(key+filename));
            }
            catch(const std::exception& e)
            {
                rsmdWARNING( "   caught exception while trying to delete " << thisPath/(key+filename) << ": " << e.what() );
            }
            
        }
    }
}



//
// helper functions
//
//     grompp -f mdp.mdp -c gro.gro -p top.top -o tpr.tpr
void EngineGMX::grompp( const std::string& mdp, const std::string& top, const std::string& gro, const std::string& tpr )
{
    execute( executablePath.c_str(), executablePath.c_str(), "grompp", 
            "-f", mdp.c_str(), 
            "-p", (top + ".top").c_str(), 
            "-c", (gro + ".gro").c_str(), 
            "-o", (tpr + ".tpr").c_str(), 
            "-po", (tpr + "-mdpout.mdp").c_str(), 
            "-quiet", "-nocopyright", backupPolicy.c_str() );
}

//     grompp -f mdp.mdp -c gro.gro -p top.top -o tpr.tpr -n ndx.ndx
void EngineGMX::grompp( const std::string& mdp, const std::string& top, const std::string& gro, const std::string& tpr, const std::string& ndx )
{
    execute( executablePath.c_str(), executablePath.c_str(), "grompp", 
            "-f", mdp.c_str(), 
            "-p", (top + ".top").c_str(), 
            "-c", (gro + ".gro").c_str(), 
            "-o", (tpr + ".tpr").c_str(), 
            "-n", (ndx + ".ndx").c_str(),
            "-quiet", "-nocopyright", backupPolicy.c_str() );       
}

//     convert-tpr -s tpr.tpr -o tpr_new.tpr -extend time
void EngineGMX::convert_tpr( const std::string& tpr, const std::string& tpr_new )
{
    execute( executablePath.c_str(), executablePath.c_str(), "convert-tpr", 
            "-s", (tpr + ".tpr").c_str(), 
            "-o", (tpr_new + ".tpr").c_str(), 
            "-extend", extensionTime_str.c_str(), 
            "-quiet", "-nocopyright", backupPolicy.c_str() );       
}

//     convert-tpr -s tpr.tpr -o tpr_new.tpr -n ndx.ndx
void EngineGMX::convert_tpr( const std::string& tpr, const std::string& tpr_new, const std::string& ndx )
{
    execute( executablePath.c_str(), executablePath.c_str(), "convert-tpr", 
            "-s", (tpr + ".tpr").c_str(), 
            "-o", (tpr_new + ".tpr").c_str(), 
            "-n", (ndx + ".ndx").c_str(), 
            "-quiet", "-nocopyright", backupPolicy.c_str() );       
}

//     trjconv -s tpr.tpr -n ndx.ndx -f gro.gro -o gro_new.gro 
void EngineGMX::trjconv( const std::string& tpr, const std::string& ndx, const std::string& trj_old, const std::string& trj_new )
{
    execute( executablePath.c_str(), executablePath.c_str(), "trjconv", 
            "-s", (tpr + ".tpr").c_str(), 
            "-n", (ndx + ".ndx").c_str(), 
            "-f", trj_old.c_str(), 
            "-o", trj_new.c_str(), 
            "-quiet", "-nocopyright", backupPolicy.c_str() );       
}

//      mdrun -s tpr.tpr -deffnm tpr
void EngineGMX::mdrun( const std::string& tpr )
{
    execute( executablePath.c_str(), executablePath.c_str(), "mdrun", 
            "-nt", nt_as_str.c_str(), 
            "-ntmpi", ntmpi_as_str.c_str(), 
            "-ntomp", ntomp_as_str.c_str(), 
            "-s", (tpr + ".tpr").c_str(), 
            "-deffnm", tpr.c_str(),
            "-quiet", "-nocopyright", backupPolicy.c_str() );
}

//      mdrun -s tpr.tpr -deffnm fnm -cpi cpt.cpt -append
void EngineGMX::mdrun( const std::string& tpr, const std::string& fnm, const std::string& cpt )
{
    execute( executablePath.c_str(), executablePath.c_str(), "mdrun", 
            "-nt", nt_as_str.c_str(), 
            "-ntmpi", ntmpi_as_str.c_str(), 
            "-ntomp", ntomp_as_str.c_str(), 
            "-s", (tpr + ".tpr").c_str(), 
            "-deffnm", fnm.c_str(),
            "-cpi", (cpt + ".cpt").c_str(), "-append", 
            "-quiet", "-nocopyright", backupPolicy.c_str() );
}

//      mdrun -s tpr.tpr -deffnm deffnm -rerun trj
void EngineGMX::mdrunRerun( const std::string& tpr, const std::string& trj, const std::string& fnm)
{
    execute( executablePath.c_str(), executablePath.c_str(), "mdrun", 
            "-nt", nt_as_str.c_str(), 
            "-ntmpi", ntmpi_as_str.c_str(), 
            "-ntomp", ntomp_as_str.c_str(),  
            "-s", (tpr + ".tpr").c_str(), 
            "-rerun", trj.c_str(), 
            "-e", (fnm + ".edr").c_str(), 
            "-g", (fnm + ".log").c_str(), 
            "-quiet", "-nocopyright", backupPolicy.c_str() );
}


//      energy -f edr.edr -o xvg.xvg
void EngineGMX::energy( const std::string& edr, const std::string& xvg)
{
    std::string pipeIn = "Potential\n";
    execute( pipeIn, executablePath.c_str(), executablePath.c_str(), "energy", 
            "-f", (edr + ".edr").c_str(),
            "-o", (xvg + ".xvg").c_str(), 
            "-quiet", "-nocopyright", backupPolicy.c_str() );
    
}
void EngineGMX::energySolvation( const std::string& edr, const std::string& xvg)
{
    std::string pipeIn = "Coul-SR:xxx-rest\n LJ-SR:xxx-rest\n";
    execute( pipeIn, executablePath.c_str(), executablePath.c_str(), "energy", 
            "-f", (edr + ".edr").c_str(),
            "-o", (xvg + ".xvg").c_str(), 
            "-quiet", "-nocopyright", backupPolicy.c_str() );
    
}


//
// read mdp file and compare with given input, 
// give warnings if something doesn't match
//
void EngineGMX::read_mdp( const std::string& filename )
{
    std::size_t nSteps = 0;
    REAL        dt = 0;

    std::ifstream FILE( filename );
    if( ! FILE )
    {
        rsmdCRITICAL( "could not read file '" << filename << "'");
    }

    std::string line {};
    while( std::getline(FILE, line, '\n') )
    {
        if( line.empty() ) continue;

        auto trimmedLine = enhance::trimString(line);
        if( trimmedLine[0] == ';' ) continue;

        auto splitted = enhance::splitString(line, '=');
        
        // search for dt and nsteps to get trajectory length 
        if( enhance::trimString(splitted[0]) == "nsteps" )   std::stringstream(splitted[1]) >> nSteps;
        if( enhance::trimString(splitted[0]) == "dt" )       std::stringstream(splitted[1]) >> dt;

    }

    extensionTime = nSteps * dt;
    rsmdLOG( "... reading md sequence length = " << extensionTime << " ps from '" << filename << "'");

}