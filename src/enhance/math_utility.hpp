/************************************************
 *                                              *
 *                rs@md                         *
 *    (reactive steps @ molecular dynamics )    *
 *                                              *
 ************************************************/
/* 
 Copyright 2020 Myra Biedermann
 Licensed under the Apache License, Version 2.0 

 Parts of the code within this file was modified from "math_utility.hpp"
 within github repository https://github.com/simonraschke/vesicle2.git,
 licensed under Apache License Version 2.0

 Myra Biedermann thankfully acknowledges support 
 from Simon Raschke.
*/

#pragma once

#include "definitions.hpp"
#include "container/atom.hpp"

#include <cmath>

// 
// some useful math functions
//

namespace enhance
{

    //
    // convert rad to deg and vice versa
    // expects floating point type
    //
    template<typename T, typename ENABLER = typename std::enable_if<std::is_floating_point<T>::value>::type>
    constexpr T deg2rad(const T& __deg) noexcept
    {
        return __deg*M_PI/180;
    }

    template<typename T, typename ENABLER = typename std::enable_if<std::is_floating_point<T>::value>::type>
    constexpr T rad2deg(const T& __rad) noexcept
    {
        return __rad/M_PI*180;
    }


    // 
    // compute the (unit) normal vector to two vectors
    //
    REALVEC normalVector(const REALVEC& v1, const REALVEC& v2);

    //
    // calculate (pbc-corrected) distance between two points / two atoms
    //
    REALVEC distanceVector(const REALVEC& v1, const REALVEC& v2, const REALVEC& box);
    REALVEC distanceVector(const Atom& a1, const Atom& a2, const REALVEC& box);

    REAL distance(const REALVEC& v1, const REALVEC& v2, const REALVEC& box);
    REAL distance(const Atom& a1, const Atom& a2, const REALVEC& box);


    //
    // calculate angle between two vectors 
    // (assumes that they are pbc-corrected)
    //
    REAL angle(const REALVEC& v1, const REALVEC& v2);


    //
    // calculate (pbc-corrected) angles between three points / three atoms
    //      1 -- 2
    //            \. 
    //             3
    //
    REAL angle(const REALVEC& p1, const REALVEC& p2, const REALVEC& p3, const REALVEC& box);
    REAL angle(const Atom& a1, const Atom& a2, const Atom& a3, const REALVEC& box);


    //
    // calculate (pbc-corrected) dihedral angle between four points / four atoms
    //      1 -- 2
    //            \.
    //             3 -- 4
    // (source: https://en.wikipedia.org/wiki/Dihedral_angle)
    //
    REAL dihedral(const REALVEC& p1, const REALVEC& p2, const REALVEC& p3, const REALVEC& p4, const REALVEC& box);
    REAL dihedral(const Atom& a1, const Atom& a2, const Atom& a3, const Atom& a4, const REALVEC& box);

}