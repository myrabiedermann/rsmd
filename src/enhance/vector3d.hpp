
#pragma once

#include <array>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>


//
// a custom 3d vector class
//
// inspired by Eigen::Matrix<REAL, 3, 1> 
//

namespace enhance
{
    template<typename T = float> 
    class Vector3d
    {
      protected:
        std::array<T,3> data;

      public:
        // 
        // constructors
        //
        Vector3d();
        Vector3d(const T);
        Vector3d(const T, const T, const T);
        
        Vector3d(const Vector3d<T>& other) = default;
        Vector3d(Vector3d<T>&& other) = default;

        template<typename O>
        Vector3d(const Vector3d<O>&);
        
        //
        // destructor
        //
        ~Vector3d<T>() = default;
        
        //
        // operators
        //
        T& operator()(std::size_t i);
        constexpr T operator()(std::size_t i) const;
        
        T& operator[](std::size_t i);
        constexpr T operator[](std::size_t i) const;

        Vector3d<T>& operator=(const Vector3d<T>&);
        Vector3d<T>& operator=(Vector3d<T>&&);

        template<typename O>
        const Vector3d<T> operator*(const O&) const;

        template<typename O>
        Vector3d<T> operator*(O&&) const;

        template<typename O>
        const Vector3d<T> operator/(const O&) const;

        template<typename O>
        Vector3d<T> operator/(O&&) const;

        template<typename O>
        const Vector3d<T> operator+(const Vector3d<O>&) const ;
       
        template<typename O>
        Vector3d<T> operator+(Vector3d<O>&&) const;

        template<typename O>
        const Vector3d<T> operator-(const Vector3d<O>&) const;
     
        template<typename O>
        Vector3d<T> operator-(Vector3d<O>&&) const;

        const Vector3d<T> operator-() const;

        template<typename O>
        Vector3d<T> operator*=(const O&);
        
        template<typename O>
        Vector3d<T> operator/=(const O&);

        template<typename O>
        Vector3d<T> operator+=(const Vector3d<O>&);

        template<typename O>
        Vector3d<T> operator-=(const Vector3d<O>&);

        bool operator==(const Vector3d<T>&);
        bool operator!=(const Vector3d<T>&);

        //
        // iterators
        //
        inline auto begin()         { return std::begin(data); };
        inline auto end()           { return std::end(data); };
    
        inline auto begin()   const { return std::begin(data); };
        inline auto end()     const { return std::end(data); };
    
        inline auto cbegin()  const { return std::cbegin(data); };
        inline auto cend()    const { return std::cend(data); };
    
        inline auto rbegin()        { return std::rbegin(data); };
        inline auto rend()          { return std::rend(data); };
    
        inline auto rbegin()  const { return std::rbegin(data); };
        inline auto rend()    const { return std::rend(data); };
    
        inline auto crbegin() const { return std::crbegin(data); };
        inline auto crend()   const { return std::crend(data); };

        //
        // some other useful member functions
        //
        float norm() const;

        template<typename O>
        T dot(const Vector3d<O>&) const;

        template<typename O>
        Vector3d<T> cross(const Vector3d<O>&) const;

        void setZero();
        bool isZero() const;

        // 
        // befriend operator to be able to do scalar * vector
        //
        template<typename O>
        friend inline Vector3d<T> operator*(O scalar, Vector3d<T> vector) 
        {
            return (vector * scalar);
        }

        //
        // write to stream
        //
        friend inline std::ostream& operator<<(std::ostream& os, const Vector3d<T>& vec)
        {
            os << "[" 
            << vec(0) << ", "
            << vec(1) << ", "
            << vec(2) << "]";
            
            return os;
        }
    };



    // 
    // constructors
    //
    template<typename T>
    Vector3d<T>::Vector3d()
        : data( {static_cast<T>(0), static_cast<T>(0), static_cast<T>(0)} )
    {}

    template<typename T>
    Vector3d<T>::Vector3d(const T v)
        : data( {v, v, v} )
    {}
    
    template<typename T>
    Vector3d<T>::Vector3d(const T v1, const T v2, const T v3)
        : data( {v1, v2, v3} )
    {}
    
    template<typename T>
    template<typename O>
    Vector3d<T>::Vector3d(const Vector3d<O>& other)
        : data( {static_cast<T>(other(0)), static_cast<T>(other(1)), static_cast<T>(other(2))} )
    {}



    // 
    // operators
    //
    template<typename T>
    T& Vector3d<T>::operator()(std::size_t i)
    {
        return data[i];
    }

    template<typename T>
    constexpr T Vector3d<T>::operator()(std::size_t i) const
    {
        return data[i];
    }

    template<typename T>
    T& Vector3d<T>::operator[](std::size_t i)
    {
        return data[i];
    }

    template<typename T>
    constexpr T Vector3d<T>::operator[](std::size_t i) const
    {
        return data[i];
    }

    template<typename T>
    Vector3d<T>& Vector3d<T>::operator=(const Vector3d<T>& other)
    {   
        data = other.data;
        return *this;
    }

    template<typename T>
    Vector3d<T>& Vector3d<T>::operator=(Vector3d<T>&& other)
    {
        data = std::move(other.data);
        return *this;
    }

    template<typename T>
    template<typename O>
    const Vector3d<T> Vector3d<T>::operator*(const O& scalar) const
    {
        Vector3d<T> copy(0);
        std::transform(std::begin(data), std::end(data), std::begin(copy), [&scalar](const T x){ return x * scalar; });
        return copy;
    }

    template<typename T>
    template<typename O>
    Vector3d<T> Vector3d<T>::operator*(O&& scalar) const
    {
        Vector3d<T> copy(0);
        std::transform(std::begin(data), std::end(data), std::begin(copy), [&scalar](const T x){ return x * scalar; });
        return copy;
    }

    template<typename T>
    template<typename O>
    const Vector3d<T> Vector3d<T>::operator/(const O& scalar) const
    {
        Vector3d<T> copy(0);
        std::transform(std::begin(data), std::end(data), std::begin(copy), [&scalar](const T x){ return x / scalar; });
        return copy;
    }

    template<typename T>
    template<typename O>
    Vector3d<T> Vector3d<T>::operator/(O&& scalar) const
    {
        Vector3d<T> copy(0);
        std::transform(std::begin(data), std::end(data), std::begin(copy), [&scalar](const T x){ return x / scalar; });
        return copy;
    }

    template<typename T>
    template<typename O>
    const Vector3d<T> Vector3d<T>::operator+(const Vector3d<O>& other) const
    {
        Vector3d<T> copy(0);
        std::transform(std::begin(data), std::end(data), std::begin(other), std::begin(copy), std::plus<T>());
        return copy;
    }

    template<typename T>
    template<typename O>
    Vector3d<T> Vector3d<T>::operator+(Vector3d<O>&& other) const
    {
        Vector3d<T> copy(0);
        std::transform(std::begin(data), std::end(data), std::begin(other), std::begin(copy), std::plus<T>());
        return copy;
    }

    template<typename T>
    template<typename O>
    const Vector3d<T> Vector3d<T>::operator-(const Vector3d<O>& other) const
    {   
        Vector3d<T> copy(0);
        std::transform(std::begin(data), std::end(data), std::begin(other), std::begin(copy), std::minus<T>());
        return copy;
    }
  
    template<typename T>
    template<typename O>
    Vector3d<T> Vector3d<T>::operator-(Vector3d<O>&& other) const
    {   
        Vector3d<T> copy(0);
        std::transform(std::begin(data), std::end(data), std::begin(other), std::begin(copy), std::minus<T>());
        return copy;
    }

    template<typename T>
    const Vector3d<T> Vector3d<T>::operator-() const
    {
        Vector3d<T> copy(*this);
        std::transform(std::begin(data), std::end(data), std::begin(copy.data), std::negate<T>());
        return copy;
    }

    template<typename T>
    template<typename O>
    Vector3d<T> Vector3d<T>::operator*=(const O& scalar)
    {
        *this = *this * scalar;
        return *this;
    }

    template<typename T>
    template<typename O>
    Vector3d<T> Vector3d<T>::operator/=(const O& scalar)
    {
        *this = *this / scalar;
        return *this;
    }

    template<typename T>
    template<typename O>
    Vector3d<T> Vector3d<T>::operator+=(const Vector3d<O>& other)
    {
        *this = *this + other;
        return *this;
    }

    template<typename T>
    template<typename O>
    Vector3d<T> Vector3d<T>::operator-=(const Vector3d<O>& other)
    {
        *this = *this - other;
        return *this;
    }

    template<typename T>
    bool Vector3d<T>::operator==(const Vector3d<T>& other)
    {   
        return this == std::addressof(other);
    }

    template<typename T>
    bool Vector3d<T>::operator!=(const Vector3d<T>& other)
    {
        return this != std::addressof(other);
    }


    //
    // some other useful member functions
    //
    template<typename T>
    float Vector3d<T>::norm() const
    {   
        float squared = std::inner_product(std::begin(data), std::end(data), std::begin(data), static_cast<float>(0));
        return std::sqrt(squared); 
    }

    template<typename T>
    template<typename O>
    T Vector3d<T>::dot(const Vector3d<O>& other) const
    {
        return std::inner_product(std::begin(data), std::end(data), std::begin(other), static_cast<T>(0)); 
    }

    template<typename T>
    template<typename O>
    Vector3d<T> Vector3d<T>::cross(const Vector3d<O>& other) const
    {
        Vector3d<T> copy(0);
        copy(0) = data[1] * other[2] - data[2] * other[1];
        copy(1) = data[2] * other[0] - data[0] * other[2];
        copy(2) = data[0] * other[1] - data[1] * other[0];
        return copy;
    }

    template<typename T>
    void Vector3d<T>::setZero()
    {
        std::fill(std::begin(data), std::end(data), static_cast<T>(0));
    }

    template<typename T>
    bool Vector3d<T>::isZero() const
    {
        return std::all_of(std::begin(data), std::end(data), [](const T& x){ return x == static_cast<T>(0); });
    }

}

