//
//  TCVec2f.h
//  threecpp
//
//  Created by CastingJ on 16/9/21.
//
//

#ifndef __threecpp__TCVec2f__
#define __threecpp__TCVec2f__

#include "TCMath.h"

class TCVec2f
{
public:
    
    /** Data type of vector components.*/
    typedef float value_type;
    
    /** Number of vector components. */
    enum { num_components = 2 };
    
    /** Vec member variable. */
    value_type _v[2];
    
    
    /** Constructor that sets all components of the vector to zero */
    TCVec2f() {_v[0]=0.0; _v[1]=0.0;}
    TCVec2f(value_type x,value_type y) { _v[0]=x; _v[1]=y; }
    
    
    inline bool operator == (const TCVec2f& v) const { return _v[0]==v._v[0] && _v[1]==v._v[1]; }
    
    inline bool operator != (const TCVec2f& v) const { return _v[0]!=v._v[0] || _v[1]!=v._v[1]; }
    
    inline bool operator <  (const TCVec2f& v) const
    {
        if (_v[0]<v._v[0]) return true;
        else if (_v[0]>v._v[0]) return false;
        else return (_v[1]<v._v[1]);
    }
    
    inline value_type * ptr() { return _v; }
    inline const value_type * ptr() const { return _v; }
    
    inline void set( value_type x, value_type y ) { _v[0]=x; _v[1]=y; }
    
    inline value_type & operator [] (int i) { return _v[i]; }
    inline value_type operator [] (int i) const { return _v[i]; }
    
    inline value_type & x() { return _v[0]; }
    inline value_type & y() { return _v[1]; }
    
    inline value_type x() const { return _v[0]; }
    inline value_type y() const { return _v[1]; }
    
    /** Returns true if all components have values that are not NaN. */
    inline bool valid() const { return !(isNaN(_v[0]) || isNaN(_v[1])); }
    
    /** Dot product. */
    inline value_type operator * (const TCVec2f& rhs) const
    {
        return _v[0]*rhs._v[0]+_v[1]*rhs._v[1];
    }
    
    /** Multiply by scalar. */
    inline const TCVec2f operator * (value_type rhs) const
    {
        return TCVec2f(_v[0]*rhs, _v[1]*rhs);
    }
    
    /** Unary multiply by scalar. */
    inline TCVec2f& operator *= (value_type rhs)
    {
        _v[0]*=rhs;
        _v[1]*=rhs;
        return *this;
    }
    
    /** Divide by scalar. */
    inline const TCVec2f operator / (value_type rhs) const
    {
        return TCVec2f(_v[0]/rhs, _v[1]/rhs);
    }
    
    /** Unary divide by scalar. */
    inline TCVec2f& operator /= (value_type rhs)
    {
        _v[0]/=rhs;
        _v[1]/=rhs;
        return *this;
    }
    
    /** Binary vector add. */
    inline const TCVec2f operator + (const TCVec2f& rhs) const
    {
        return TCVec2f(_v[0]+rhs._v[0], _v[1]+rhs._v[1]);
    }
    
    /** Unary vector add. Slightly more efficient because no temporary
     * intermediate object.
     */
    inline TCVec2f& operator += (const TCVec2f& rhs)
    {
        _v[0] += rhs._v[0];
        _v[1] += rhs._v[1];
        return *this;
    }
    
    /** Binary vector subtract. */
    inline const TCVec2f operator - (const TCVec2f& rhs) const
    {
        return TCVec2f(_v[0]-rhs._v[0], _v[1]-rhs._v[1]);
    }
    
    /** Unary vector subtract. */
    inline TCVec2f& operator -= (const TCVec2f& rhs)
    {
        _v[0]-=rhs._v[0];
        _v[1]-=rhs._v[1];
        return *this;
    }
    
    /** Negation operator. Returns the negative of the TCVec2f. */
    inline const TCVec2f operator - () const
    {
        return TCVec2f (-_v[0], -_v[1]);
    }
    
    /** Length of the vector = sqrt( vec . vec ) */
    inline value_type length() const
    {
        return sqrtf( _v[0]*_v[0] + _v[1]*_v[1] );
    }
    
    /** Length squared of the vector = vec . vec */
    inline value_type length2( void ) const
    {
        return _v[0]*_v[0] + _v[1]*_v[1];
    }
    
    /** Normalize the vector so that it has length unity.
     * Returns the previous length of the vector.
     */
    inline value_type normalize()
    {
        value_type norm = TCVec2f::length();
        if (norm>0.0)
        {
            value_type inv = 1.0f/norm;
            _v[0] *= inv;
            _v[1] *= inv;
        }
        return( norm );
    }
    
};    // end of class TCVec2f

/** multiply by vector components. */
inline TCVec2f componentMultiply(const TCVec2f& lhs, const TCVec2f& rhs)
{
    return TCVec2f(lhs[0]*rhs[0], lhs[1]*rhs[1]);
}

/** divide rhs components by rhs vector components. */
inline TCVec2f componentDivide(const TCVec2f& lhs, const TCVec2f& rhs)
{
    return TCVec2f(lhs[0]/rhs[0], lhs[1]/rhs[1]);
}


#endif /* defined(__threecpp__TCVec2f__) */
