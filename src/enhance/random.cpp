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

#include "random.hpp"



enhance::RandomEngineInit::RandomEngineInit()
    : seed( true_engine() )
{
    pseudo_engine.seed(seed);
}


//
// uniform_real_distribution returns random real from [a,b)
// uniform_int_distribution  returns random int from [a,b]
//
template<>
float enhance::random(const float& a, const float& b)
{
    std::uniform_real_distribution<float> dist(a,b); 
    return dist(enhance::RandomEngine.pseudo_engine);
}

template<>
double enhance::random(const double& a, const double& b)
{
    std::uniform_real_distribution<double> dist(a,b);  
    return dist(enhance::RandomEngine.pseudo_engine);
}

template<>
unsigned int enhance::random(const unsigned int& a, const unsigned int& b)
{
    std::uniform_int_distribution<unsigned int> dist(a,b);  
    return dist(enhance::RandomEngine.pseudo_engine);
}

template<>
int enhance::random(const int& a, const int& b)
{
    std::uniform_int_distribution<int> dist(a,b);  
    return dist(enhance::RandomEngine.pseudo_engine);
}

template<>
std::size_t enhance::random(const std::size_t& a, const std::size_t& b)
{
    std::uniform_int_distribution<std::size_t> dist(a,b);  
    return dist(enhance::RandomEngine.pseudo_engine);
}