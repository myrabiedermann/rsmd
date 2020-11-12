
#pragma once

#include "enhance/vector3d.hpp"


//
// some useful typedefs
//
typedef float                     REAL;
typedef enhance::Vector3d<float>  REALVEC;


//
// log errors, warnings etc.
//
#include <csignal>
#include <iostream>

static std::string rsmdALL_formatting       {"          "};
static std::string rsmdLOG_formatting       {"  [LOG]   "};
static std::string rsmdDEBUG_formatting     {" [DEBUG]  "};
static std::string rsmdWARNING_formatting   {"[WARNING] "};
static std::string rsmdCRITICAL_formatting  {" [ERROR]  "};

#ifndef NDEBUG
    #define rsmdDEBUG(x) {std::cerr << rsmdDEBUG_formatting; do { std::cerr << x; } while (0); std::cerr << '\n';}
#else
    #define rsmdDEBUG(x)
#endif
#define rsmdLOG(x)       {std::cout << rsmdLOG_formatting; do { std::cout << x; } while (0); std::cout << '\n';}
#define rsmdWARNING(x)   {std::cout << std::flush; std::clog << rsmdWARNING_formatting; do { std::clog << x; } while (0); std::clog << '\n';}
#define rsmdCRITICAL(x)  {std::cout << std::flush; std::cerr << rsmdCRITICAL_formatting << __FILE__ <<":" << __LINE__ << "  "; do { std::cerr << x; } while (0); std::cerr <<", raising SIGABRT\n"; std::raise(SIGABRT); }
#define rsmdEXIT(x)      {std::cout << std::flush; std::cerr << rsmdCRITICAL_formatting; do { std::cerr << x; } while (0);  std::cout << '\n'; std::exit(EXIT_FAILURE); }


