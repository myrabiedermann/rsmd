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

#include "math_utility.hpp"

// 
// compute the (unit) normal vector to two vectors
//
REALVEC enhance::normalVector(const REALVEC& v1, const REALVEC& v2)
{
    REALVEC normal = v1.cross(v2);
    REAL n = normal.norm();

    if( n == 0 )    // return (0,0,0) instead of (nan, nan, nan)
        return REALVEC(0);
    return (normal / n);
}


//
// calculate (pbc-corrected) distance between two points / two atoms
//
REALVEC enhance::distanceVector(const REALVEC& v1, const REALVEC& v2, const REALVEC& box)
{
    // outDEBUG( (box.isZero() ? "warning: given pbc dimensions are zero: " : "computing pbc-corrected distance with pbc dimensions ") << box );
    #ifndef NDEBUG
        if( box.isZero() )  rsmdDEBUG( "warning: given pbx dimensions are zero" );
    #endif

    REALVEC distance = v2 - v1;
    distance(0) = distance(0) - box(0) * std::round( distance(0)/box(0) );
    distance(1) = distance(1) - box(1) * std::round( distance(1)/box(1) );
    distance(2) = distance(2) - box(2) * std::round( distance(2)/box(2) );

    return distance;
}

REALVEC enhance::distanceVector(const Atom& a1, const Atom& a2, const REALVEC& box)
{
    return distanceVector(a1.position, a2.position, box);
}

REAL enhance::distance(const REALVEC& v1, const REALVEC& v2, const REALVEC& box)
{
    REALVEC distance = enhance::distanceVector(v1, v2, box);
    return distance.norm();
}

REAL enhance::distance(const Atom& a1, const Atom& a2, const REALVEC& box)
{
    return distance(a1.position, a2.position, box);
}


//
// calculate angle between two vectors 
// (assumes that they are pbc-corrected)
//
REAL enhance::angle(const REALVEC& v1, const REALVEC& v2)
{
    REAL angle = std::acos( v1.dot(v2) / (v1.norm() * v2.norm() ) );
    // return angle;
    return rad2deg(angle);
}


//
// calculate (pbc-corrected) angles between three points / three atoms
//      1 -- 2
//            \. 
//             3
//
REAL enhance::angle(const REALVEC& p1, const REALVEC& p2, const REALVEC& p3, const REALVEC& box)
{
    REALVEC vector1 = distanceVector(p1, p2, box);
    REALVEC vector2 = distanceVector(p2, p3, box);
    return angle(vector1, vector2);
}

REAL enhance::angle(const Atom& a1, const Atom& a2, const Atom& a3, const REALVEC& box)
{
    return angle(a1.position, a2.position, a3.position, box);
}


//
// calculate (pbc-corrected) dihedral angle between four points / four atoms
//      1 -- 2
//            \.
//             3 -- 4
// (source: https://en.wikipedia.org/wiki/Dihedral_angle)
//
REAL enhance::dihedral(const REALVEC& p1, const REALVEC& p2, const REALVEC& p3, const REALVEC& p4, const REALVEC& box)
{
    REALVEC vector1 = distanceVector(p1, p2, box);
    REALVEC vector2 = distanceVector(p2, p3, box);
    REALVEC vector3 = distanceVector(p3, p4, box);

    REALVEC n1 = normalVector(vector1, vector2);
    REALVEC n2 = normalVector(vector2, vector3);

    REAL x = ( n1.cross(n2) ).dot( vector2/vector2.norm() );
    REAL y = n1.dot(n2);

    REAL angle = std::atan2( x, y );

    return rad2deg(angle);
}

REAL enhance::dihedral(const Atom& a1, const Atom& a2, const Atom& a3, const Atom& a4, const REALVEC& box)
{
    return dihedral(a1.position, a2.position, a3.position, a4.position, box);
}
