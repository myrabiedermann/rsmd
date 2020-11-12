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

#include "parameters/parameters.hpp"

//
// constructor
//
Parameters::Parameters(int argc, char* argv[])
{
    read(argc, argv);
    check(); 
}


//
// read program options from commandline / from file
//
void Parameters::read(int argc, char* argv[])
{
    namespace po = boost::program_options;

    parameterMap.clear();   
    
    // add all options
    // ... general options:
    po::options_description generalOptions("General options");
    generalOptions.add_options()
        ("help,h",   "produce this help")
        ("input,i",   po::value<std::string>(&configFileName), "input file from which to read program options")
        ("output,o",  po::value<std::string>()->default_value("RESTART"), "output file where program options for a restart are written to")
        ("rseed",     po::value<std::size_t>()->default_value(0), "random seed (0: true random, else: given seed)")
        ("statistics", po::value<std::string>()->default_value("statistics.data"), "output file for statistics on reactive steps")
    ;

    // ... helper options
    po::options_description helpOptions("Tip: add one of the following options if you require additional help");
    helpOptions.add_options()
        ("credits",   po::bool_switch(), "authorship etc.")
        ("reaction",  po::bool_switch(), "get help on how to write a reaction input file")
        ("gromacs",   po::bool_switch(), "produce help for GROMACS related options")
    ;

    // ... simulation setup related options:
    po::options_description simulationOptions("Simulation setup related options");
    simulationOptions.add_options()
        ("simulation.engine",  po::value<std::string>(), "path to the MD engine executable")
        ("simulation.cycles",  po::value<std::size_t>()->default_value(1), "# of cycles")
        ("simulation.restart", po::bool_switch(), "restart simulation and append to existing simulation files")
        ("simulation.restartCycle", po::value<std::size_t>(), "restart with this cycle")
        ("simulation.restartCycleFiles", po::value<std::size_t>(), "append to simulation files named according to this cycle")
    ;
    
    // ... reaction related options:
    po::options_description reactionOptions("Reaction related options");
    reactionOptions.add_options()
        ("reaction.file", po::value<std::vector<std::string>>()->multitoken(), "reaction input files (multiple args or occurrences possible)")
        ("reaction.mc",    po::bool_switch(), "use Metropolis MC acceptance criterion")
        ("reaction.rate",  po::bool_switch(), "use rate-based acceptance criterion")
        ("reaction.frequency",   po::value<REAL>(), "attempt frequency for reactive steps \n(required if reaction.rates)")
        ("reaction.temperature", po::value<REAL>(), "simulation temperature (required if reaction.mc)" )
        ("reaction.averagePotentialEnergy", po::value<REAL>()->default_value(0.0), "time interval over which to average potential energies (only if reaction.mc)" )
        ("reaction.computeLocalPotentialEnergy", po::bool_switch(), "compute local potential energies (only if reaction.mc)")
        ("reaction.computeSolvationPotentialEnergy", po::bool_switch(), "compute solvation interaction (only if reaction.mc)")
        ("reaction.saveRejected", po::bool_switch(), "save md files from failed reactive steps instead of deleting them")
    ;

    // ... md engine related options
    po::options_description gromacsOptions("GROMACS related options");
    gromacsOptions.add_options()
        ("gromacs.topology",       po::value<std::string>(), "topology file (.top)")
        ("gromacs.coordinates",    po::value<std::string>(), "coordinates file (.gro)")
        ("gromacs.mdp",            po::value<std::string>(), "md parameter file (.mdp)")
        ("gromacs.mdp.energy",     po::value<std::string>()->default_value(""), "md parameter file for energy computation with solvation interaction (.mdp)")
        ("gromacs.mdp.relaxation", po::value<std::string>(), "md parameter file for relaxation (.mdp)")
        ("gromacs.backup",         po::bool_switch(), "whether or not gromacs should backup files or overwrite them")
        ("gromacs.nt",             po::value<int>()->default_value(0), "total number of threads to start (0 is guess)")
        ("gromacs.ntmpi",          po::value<int>()->default_value(0), "number of thread-MPI ranks to start (0 is guess)")
        ("gromacs.ntomp",          po::value<int>()->default_value(0), "number of OpenMP threads per MPI rank to start (0 is guess)")
    ;


    // ... merge all options 
    po::options_description allOptions("");
    allOptions.add(generalOptions).add(helpOptions).add(simulationOptions).add(reactionOptions).add(gromacsOptions);


    // try storing + notifying the parameterMap
    try
    {
        // ... first: parse all options from command line
        po::store( po::command_line_parser(argc, argv).options(allOptions).run(), parameterMap );
        po::notify( parameterMap );

        // ... second: parse options from input file
        if( ! configFileName.empty() )
        {
            std::ifstream FILE(configFileName);
            if( !FILE )
            {
                std::cout << "error: could not open the specified input file\n";
                std::exit(EXIT_FAILURE);
            }
            else
            {
                po::store( po::parse_config_file(FILE, allOptions), parameterMap );
                po::notify( parameterMap );   
            }
            FILE.close();
        }
    }
    catch( const std::exception& e )
    {
        std::cout << "error while parsing program options: " << e.what() << '\n';
        std::exit(EXIT_FAILURE);
    }


    // if there was a call for help
    if( parameterMap.count("help") )
    {
        if( getOption("credits").as<bool>() )
        {
            std::cout << programName
                      << addDescription
                      << "\n"
                      << licensingInfo
                      << "\n"
                      << "please read and cite: <to be published>\n";
        }
        else if( getOption("reaction").as<bool>() )
        {
            ReactionParser parser;
            std::cout << programName
                      << "# a reaction input file should comply with the following syntax:\n\n";
            std::cout << parser.writeExample() << '\n';
        }
        else if( getOption("gromacs").as<bool>() )
        {
            std::cout << programName
                      << gromacsOptions;
        }
        else
        {
            std::cout << programName
                      << addDescription
                      << '\n'
                      << licensingInfo;
            std::cout << '\n' << generalOptions;
            std::cout << '\n' << simulationOptions;
            std::cout << '\n' << reactionOptions;
            std::cout << '\n' << helpOptions;
            std::cout << '\n' << "Tip: to achieve a civilised shutdown of this program, e.g. if the runtime you \n"
                              << "     allocated for the job is about to run out, send SIGUSR1.\n";
            std::cout << '\n' << "Please note: all physical quantities must be given in the same units that are\n"
                              << "             used by the md engine with which the simulations are performed.\n";
        }
        
        std::exit(EXIT_SUCCESS);
    }   

    // if not a single option was given
    if( argc == 1 )
    {
        std::cout << programName   
                  << addDescription
                  << "\n"
                  << licensingInfo
                  << "\n"
                  << "please specify some program options!\n"
                  << "(you can use the --help flag for further information)\n";
        std::exit(EXIT_FAILURE);
    }
    
}


//
// check for completeness / workability with given program options 
//
void Parameters::check()
{
    // get mdEngine
    if( parameterMap.count("simulation.engine") )
    {
        std::string tmp = getOption("simulation.engine").as<std::string>();
        if( tmp.find("gmx") != std::string::npos )
        {
            mdEngine = ENGINE::GROMACS;
        }
        else
        {
            std::cout << "error: could not recognise md engine from given program option 'simulation.engine' \n";
            std::exit(EXIT_FAILURE);
        }
    }
    else
    {
        std::cout << "error: program option 'simulation.engine' is mandatory\n";
        std::exit(EXIT_FAILURE);
    }

    // get simulation mode
    if( getOption("simulation.restart").as<bool>() )
    {
        simulationMode = SIMMODE::RESTART;
        if( ! parameterMap.count("simulation.restartCycle") || ! parameterMap.count("simulation.restartCycleFiles") )
        {
            std::cout << "error: program options 'simulation.restartCycle' and 'simulation.restartCycleFiles' are mandatory if 'simulation.restart' is set.\n";
            std::exit(EXIT_FAILURE);
        }
    }
    // get simulation algorithm
    if( getOption("reaction.mc").as<bool>() )
    {
        simulationAlgorithm = SIMALGORITHM::MC;
    }
    else if( getOption("reaction.rate").as<bool>() )
    {
        simulationAlgorithm = SIMALGORITHM::RATE;
    }
    
    // check for mandatory, depending and conflicting options
    if( ! getOption("simulation.restart").as<bool>() && parameterMap.count("simulation.restartCycle") )
    {
        std::cout << "warning: you set 'simulation.restartCycle' but simulation.restart = off. that doesn't seem right\n";
        std::exit(EXIT_FAILURE);
    }
    if( ! getOption("simulation.restart").as<bool>() && parameterMap.count("simulation.restartCycleFiles") )
    {
        std::cout << "warning: you set 'simulation.restartCycleFiles' but simulation.restart = off. that doesn't seem right\n";
        std::exit(EXIT_FAILURE);
    }
    if( ! parameterMap.count("reaction.file") )
    {
        std::cout << "error: at least one occurrence of program option 'reaction.file' is mandatory\n";
        std::exit(EXIT_FAILURE);
    }
    if( getOption("reaction.rate").as<bool>() == getOption("reaction.mc").as<bool>() )
    {
        std::cout << "error: program options 'reaction.rate' and 'reaction.mc' are mutually exclusive, you need to set either of them\n";
        std::exit(EXIT_FAILURE);
    }
    if( getOption("reaction.rate").as<bool>() && ! parameterMap.count("reaction.frequency") )
    {
        std::cout << "error: program option 'reaction.frequency' is mandatory if 'reaction.rate' is set\n";
        std::exit(EXIT_FAILURE);
    }
    if( getOption("reaction.mc").as<bool>() && ! parameterMap.count("reaction.temperature") )
    {
        std::cout << "error: program option 'reaction.temperature' is mandatory if 'reaction.mc' is set\n";
        std::exit(EXIT_FAILURE);
    }
    if( getOption("reaction.computeSolvationPotentialEnergy").as<bool>() && ! getOption("reaction.computeLocalPotentialEnergy").as<bool>() )
    {
        std::cout << "error: computing interaction energies with solvent without setting 'reaction.computeLocalPotentialEnergy' makes no sense.\n";
        std::exit(EXIT_FAILURE);
    }

    if( mdEngine == ENGINE::GROMACS )
    {
        if( ! parameterMap.count("gromacs.topology") )
        {
            std::cout << "error: program option gromacs.topology is mandatory\n";
            std::exit(EXIT_FAILURE);
        }
        if( ! parameterMap.count("gromacs.coordinates") )
        {
            std::cout << "error: program option 'gromacs.coordinates' is mandatory\n";
            std::exit(EXIT_FAILURE);
        }
        if( ! parameterMap.count("gromacs.mdp") )
        {
            std::cout << "error: program option 'gromacs.mdp' is mandatory\n";
            std::exit(EXIT_FAILURE);
        }
        if( ! parameterMap.count("gromacs.mdp.relaxation") )
        {
            std::cout << "error: program option 'gromacs.mdp.relaxation' is mandatory\n";
            std::exit(EXIT_FAILURE);
        }
        if( getOption("reaction.computeSolvationPotentialEnergy").as<bool>() && getOption("gromacs.mdp.energy").as<std::string>().empty() )
        {
            std::cout << "error: program option 'gromacs.mdp.energy' is mandatory if 'reaction.computeSolvationPotentialEnergy' is set\n";
            std::exit(EXIT_FAILURE);
        }
    }
}



//
// print options to string
//
std::string Parameters::str() const
{
    std::stringstream stream {};

    stream << rsmdALL_formatting << "--- General options:\n";
    if( ! configFileName.empty() )
        stream << rsmdALL_formatting << formatted( "input", configFileName ) << '\n';
    else
        stream << rsmdALL_formatting << formatted( "input", "none") << '\n';
    stream << rsmdALL_formatting << formatted( "output", getOption("output").as<std::string>() ) << '\n';
    stream << rsmdALL_formatting << formatted( "statistics", getOption("statistics").as<std::string>() ) << '\n';
    stream << rsmdALL_formatting << formatted( "rseed", getOption("rseed").as<std::size_t>() ) << '\n';

    stream << rsmdALL_formatting << "--- Simulation setup related options:\n"
           << rsmdALL_formatting << formatted( "simulation.engine", getOption("simulation.engine").as<std::string>() ) << '\n'
           << rsmdALL_formatting << formatted( "simulation.cycles", getOption("simulation.cycles").as<std::size_t>() ) << '\n';
    if( getOption("simulation.restart").as<bool>() )
    {
        stream << rsmdALL_formatting << formatted( "simulation.restartCycle", getOption("simulation.restartCycle").as<std::size_t>() ) << '\n'
               << rsmdALL_formatting << formatted( "simulation.restartCycleFiles", getOption("simulation.restartCycleFiles").as<std::size_t>() ) << '\n';
    }

    stream << rsmdALL_formatting << "--- Reaction related options:\n";
    stream << rsmdALL_formatting << formatted( "reaction.file(s)", getOption("reaction.file").as<std::vector<std::string>>() ) << '\n';
    if( getOption("reaction.mc").as<bool>() )
    {
        stream << rsmdALL_formatting << formatted( "reaction.mc", getOption("reaction.mc").as<bool>() ) << '\n'
               << rsmdALL_formatting << formatted( "reaction.temperature", getOption("reaction.temperature").as<REAL>() ) << '\n'
               << rsmdALL_formatting << formatted( "reaction.averagePotentialEnergy", getOption("reaction.averagePotentialEnergy").as<REAL>() ) << '\n'
               << rsmdALL_formatting << formatted( "reaction.computeLocalPotentialEnergy", getOption("reaction.computeLocalPotentialEnergy").as<bool>() ) << '\n'
               << rsmdALL_formatting << formatted( "reaction.computeSolvationPotentialEnergy", getOption("reaction.computeSolvationPotentialEnergy").as<bool>() ) << '\n';
    }
    else if( getOption("reaction.rate").as<bool>() )
    {
        stream << rsmdALL_formatting << formatted( "reaction.rate", getOption("reaction.rate").as<bool>() ) << '\n'
               << rsmdALL_formatting << formatted( "reaction.frequency", getOption("reaction.frequency").as<REAL>() ) << '\n';
    }
    stream << rsmdALL_formatting << formatted( "saveRejected", getOption("reaction.saveRejected").as<bool>() ) << '\n';

    if( mdEngine == ENGINE::GROMACS )
    {
        stream << rsmdALL_formatting << "--- GROMACS related options:\n"
               << rsmdALL_formatting << formatted("gromacs.topology", getOption("gromacs.topology").as<std::string>() ) << '\n'
               << rsmdALL_formatting << formatted("gromacs.coordinates", getOption("gromacs.coordinates").as<std::string>() ) << '\n'
               << rsmdALL_formatting << formatted("gromacs.mdp", getOption("gromacs.mdp").as<std::string>() ) << '\n'
               << rsmdALL_formatting << formatted("gromacs.mdp.relaxation", getOption("gromacs.mdp.relaxation").as<std::string>() ) << '\n';
        if( getOption("reaction.computeSolvationPotentialEnergy").as<bool>() )
        {
            stream << rsmdALL_formatting << formatted("gromacs.mdp.energy", getOption("gromacs.mdp.energy").as<std::string>() ) << '\n';
        }
        
        stream << rsmdALL_formatting << formatted("gromacs.backup", getOption("gromacs.backup").as<bool>() ) << '\n'
               << rsmdALL_formatting << formatted("gromacs.nt", getOption("gromacs.nt").as<int>() ) << '\n'
               << rsmdALL_formatting << formatted("gromacs.ntmpi", getOption("gromacs.ntmpi").as<int>() ) << '\n'
               << rsmdALL_formatting << formatted("gromacs.ntomp", getOption("gromacs.ntomp").as<int>() ) << '\n';
    }
    
    return stream.str();
}


