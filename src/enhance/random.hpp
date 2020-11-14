/************************************************
 *                                              *
 *                rs@md                         *
 *    (reactive steps @ molecular dynamics )    *
 *                                              *
 ************************************************/
/* 
 Copyright 2020 Myra Biedermann
 Licensed under the Apache License, Version 2.0 

 Parts of the code within this file was modified from "random.hpp"
 within github repository https://github.com/simonraschke/vesicle2.git,
 licensed under Apache License Version 2.0

 Myra Biedermann thankfully acknowledges support 
 from Simon Raschke.
*/
#pragma once

#include <random>
#include <algorithm>
#include <iostream>

// 
// random number generator and random iterator utility
//

namespace enhance
{
    static struct RandomEngineInit
    {
        RandomEngineInit();
        auto getSeed()         const { return seed; };
        void setSeed(unsigned int s) { seed = s; pseudo_engine.seed(seed); };
        std::mt19937_64 pseudo_engine {};

      private:
        std::random_device true_engine {};
        unsigned int seed {};

    }RandomEngine;


    // call these functions to get random number
    template<typename T>
    T random(const T& a, const T& b);



    // shuffle randomly
    template<class D>
    void shuffle(D first, D last)
    {
        std::shuffle(first, last, enhance::RandomEngine.pseudo_engine);
    }

    // shuffle randomly with associated weights
    template<class D, class W>
    void weighted_shuffle(D first, D last, W first_weight, W last_weight)
    {
        while (first != last and first_weight != last_weight)
        {
            std::discrete_distribution<> dd(first_weight, last_weight);
            auto i = dd(enhance::RandomEngine.pseudo_engine);
            if( i )
            {
                std::iter_swap(first, std::next(first, i));
                std::iter_swap(first_weight, std::next(first_weight, i));
            }
            ++first;
            ++first_weight;
        }
    }


    // random pick from a sequence
    template<class D>
    D random_choice(D first, D last)
    {
        int nChoices = std::distance(first, last);
        auto choice = enhance::random(0, nChoices-1);
        return std::next(first, choice);
    } 

    // random weighted pick from a sequence
    template<class D, class W>
    D random_weighted_choice(D first, W first_weight, W last_weight)
    {
        std::discrete_distribution<int> dd(first_weight, last_weight);
        auto choice = dd(enhance::RandomEngine.pseudo_engine);

        return std::next(first, choice);
    }


    // Type aliasing
    template<typename Container>
    using IteratorCategoryOf = typename std::iterator_traits<typename Container::iterator>::iterator_category;


    // Compile time iterator check
    template<typename Container>
    using HaveRandomAccessIterator = std::is_base_of< std::random_access_iterator_tag, IteratorCategoryOf<Container>>;


    // random iterator
    static struct _randomIterator
    {
        template<typename T>
        constexpr inline auto operator() ( const T& _container ) const -> typename T::const_iterator
        {
            static_assert( HaveRandomAccessIterator<T>::value, "T has no std::random_access_iterator_tag in __random_iterator::operator()" );
            std::uniform_int_distribution<std::size_t> dist(0,_container.size()-1);
            return std::cbegin(_container) + dist(RandomEngine.pseudo_engine);
        }

        template<typename T>
        constexpr inline auto operator() ( T& _container ) -> typename T::iterator
        {
            static_assert( HaveRandomAccessIterator<T>::value, "T has no std::random_access_iterator_tag in __random_iterator::operator()" );
            std::uniform_int_distribution<std::size_t> dist(0,_container.size()-1);
            return std::begin(_container) + dist(RandomEngine.pseudo_engine);
        }
    } randomIterator __attribute__((unused));

}