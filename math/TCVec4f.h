//
//  TCVec4f.h
//  threecpp
//
//  Created by CastingJ on 16/9/21.
//
//

#ifndef __threecpp__TCVec4f__
#define __threecpp__TCVec4f__

#include "TCVec3f.h"

class TCVec4f
{
public:
    
    /** Data type of vector components.*/
    typedef float value_type;
    
    /** Number of vector components. */
    enum { num_components = 4 };
    
    /** Vec member variable. */
    value_type _v[4];
    
    // Methods are defined here so that they are implicitly inlined
    
    /** Constructor that sets all components of the vector to zero */
    TCVec4f() { _v[0]=0.0f; _v[1]=0.0f; _v[2]=0.0f; _v[3]=0.0f;}
    
    TCVec4f(value_type x, value_type y, value_type z, value_type w)
    {
        _v[0]=x;
        _v[1]=y;
        _v[2]=z;
        _v[3]=w;
    }
    
    TCVec4f(const TCVec3f& v3,value_type w)
    {
        _v[0]=v3[0];
        _v[1]=v3[1];
        _v[2]=v3[2];
        _v[3]=w;
    }
    
    inline bool operator == (const TCVec4f& v) const { return _v[0]==v._v[0] && _v[1]==v._v[1] && _v[2]==v._v[2] && _v[3]==v._v[3]; }
    
    inline bool operator != (const TCVec4f& v) const { return _v[0]!=v._v[0] || _v[1]!=v._v[1] || _v[2]!=v._v[2] || _v[3]!=v._v[3]; }
    
    inline bool operator <  (const TCVec4f& v) const
    {
        if (_v[0]<v._v[0]) return true;
        else if (_v[0]>v._v[0]) return false;
        else if (_v[1]<v._v[1]) return true;
        else if (_v[1]>v._v[1]) return false;
        else if (_v[2]<v._v[2]) return true;
        else if (_v[2]>v._v[2]) return false;
        else return (_v[3]<v._v[3]);
    }
    
    inline value_type* ptr() { return _v; }
    inline const value_type* ptr() const { return _v; }
    
    inline void set( value_type x, value_type y, value_type z, value_type w)
    {
        _v[0]=x; _v[1]=y; _v[2]=z; _v[3]=w;
    }
    
    inline value_type& operator [] (unsigned int i) { return _v[i]; }
    inline value_type  operator [] (unsigned int i) const { return _v[i]; }
    
    inline value_type& x() { return _v[0]; }
    inline value_type& y() { return _v[1]; }
    inline value_type& z() { return _v[2]; }
    inline value_type& w() { return _v[3]; }
    
    inline value_type x() const { return _v[0]; }
    inline value_type y() const { return _v[1]; }
    inline value_type z() const { return _v[2]; }
    inline value_type w() const { return _v[3]; }
    
    inline value_type& r() { return _v[0]; }
    inline value_type& g() { return _v[1]; }
    inline value_type& b() { return _v[2]; }
    inline value_type& a() { return _v[3]; }
    
    inline value_type r() const { return _v[0]; }
    inline value_type g() const { return _v[1]; }
    inline value_type b() const { return _v[2]; }
    inline value_type a() const { return _v[3]; }
    
    inline unsigned int asABGR() const
    {
        return (unsigned int)clampTo((_v[0]*255.0f),0.0f,255.0f)<<24 |
        (unsigned int)clampTo((_v[1]*255.0f),0.0f,255.0f)<<16 |
        (unsigned int)clampTo((_v[2]*255.0f),0.0f,255.0f)<<8  |
        (unsigned int)clampTo((_v[3]*255.0f),0.0f,255.0f);
    }
    
    inline unsigned int asRGBA() const
    {
        return (unsigned int)clampTo((_v[3]*255.0f),0.0f,255.0f)<<24 |
        (unsigned int)clampTo((_v[2]*255.0f),0.0f,255.0f)<<16 |
        (unsigned int)clampTo((_v[1]*255.0f),0.0f,255.0f)<<8  |
        (unsigned int)clampTo((_v[0]*255.0f),0.0f,255.0f);
    }
    
    /** Returns true if all components have values that are not NaN. */
    inline bool valid() const { return !(isNaN(_v[0]) || isNaN(_v[1]) || isNaN(_v[2]) || isNaN(_v[3])); }
    
    /** Dot product. */
    inline value_type operator * (const TCVec4f& rhs) const
    {
        return _v[0]*rhs._v[0]+
        _v[1]*rhs._v[1]+
        _v[2]*rhs._v[2]+
        _v[3]*rhs._v[3] ;
    }
    
    /** Multiply by scalar. */
    inline TCVec4f operator * (value_type rhs) const
    {
        return TCVec4f(_v[0]*rhs, _v[1]*rhs, _v[2]*rhs, _v[3]*rhs);
    }
    
    /** Unary multiply by scalar. */
    inline TCVec4f& operator *= (value_type rhs)
    {
        _v[0]*=rhs;
        _v[1]*=rhs;
        _v[2]*=rhs;
        _v[3]*=rhs;
        return *this;
    }
    
    /** Divide by scalar. */
    inline TCVec4f operator / (value_type rhs) const
    {
        return TCVec4f(_v[0]/rhs, _v[1]/rhs, _v[2]/rhs, _v[3]/rhs);
    }
    
    /** Unary divide by scalar. */
    inline TCVec4f& operator /= (value_type rhs)
    {
        _v[0]/=rhs;
        _v[1]/=rhs;
        _v[2]/=rhs;
        _v[3]/=rhs;
        return *this;
    }
    
    /** Binary vector add. */
    inline TCVec4f operator + (const TCVec4f& rhs) const
    {
        return TCVec4f(_v[0]+rhs._v[0], _v[1]+rhs._v[1],
                     _v[2]+rhs._v[2], _v[3]+rhs._v[3]);
    }
    
    /** Unary vector add. Slightly more efficient because no temporary
     * intermediate object.
     */
    inline TCVec4f& operator += (const TCVec4f& rhs)
    {
        _v[0] += rhs._v[0];
        _v[1] += rhs._v[1];
        _v[2] += rhs._v[2];
        _v[3] += rhs._v[3];
        return *this;
    }
    
    /** Binary vector subtract. */
    inline TCVec4f operator - (const TCVec4f& rhs) const
    {
        return TCVec4f(_v[0]-rhs._v[0], _v[1]-rhs._v[1],
                     _v[2]-rhs._v[2], _v[3]-rhs._v[3] );
    }
    
    /** Unary vector subtract. */
    inline TCVec4f& operator -= (const TCVec4f& rhs)
    {
        _v[0]-=rhs._v[0];
        _v[1]-=rhs._v[1];
        _v[2]-=rhs._v[2];
        _v[3]-=rhs._v[3];
        return *this;
    }
    
    /** Negation operator. Returns the negative of the TCVec4f. */
    inline const TCVec4f operator - () const
    {
        return TCVec4f (-_v[0], -_v[1], -_v[2], -_v[3]);
    }
    
    /** Length of the vector = sqrt( vec . vec ) */
    inline value_type length() const
    {
        return sqrtf( _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3]);
    }
    
    /** Length squared of the vector = vec . vec */
    inline value_type length2() const
    {
        return _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3];
    }
    
    /** Normalize the vector so that it has length unity.
     * Returns the previous length of the vector.
     */
    inline value_type normalize()
    {
        value_type norm = TCVec4f::length();
        if (norm>0.0f)
        {
            value_type inv = 1.0f/norm;
            _v[0] *= inv;
            _v[1] *= inv;
            _v[2] *= inv;
            _v[3] *= inv;
        }
        return( norm );
    }
    
};    // end of class TCVec4f

/** Compute the dot product of a (Vec3,1.0) and a TCVec4f. */
inline TCVec4f::value_type operator * (const TCVec3f& lhs,const TCVec4f& rhs)
{
    return lhs[0]*rhs[0]+lhs[1]*rhs[1]+lhs[2]*rhs[2]+rhs[3];
}

/** Compute the dot product of a TCVec4f and a (Vec3,1.0). */
inline TCVec4f::value_type operator * (const TCVec4f& lhs,const TCVec3f& rhs)
{
    return lhs[0]*rhs[0]+lhs[1]*rhs[1]+lhs[2]*rhs[2]+lhs[3];
}

/** multiply by vector components. */
inline TCVec4f componentMultiply(const TCVec4f& lhs, const TCVec4f& rhs)
{
    return TCVec4f(lhs[0]*rhs[0], lhs[1]*rhs[1], lhs[2]*rhs[2], lhs[3]*rhs[3]);
}

/** divide rhs components by rhs vector components. */
inline TCVec4f componentDivide(const TCVec4f& lhs, const TCVec4f& rhs)
{
    return TCVec4f(lhs[0]/rhs[0], lhs[1]/rhs[1], lhs[2]/rhs[2], lhs[3]/rhs[3]);
}


#endif /* defined(__threecpp__TCVec4f__) */
