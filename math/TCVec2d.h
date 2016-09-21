//
//  TCVec2d.h
//  threecpp
//
//  Created by CastingJ on 16/9/21.
//
//

#ifndef __threecpp__TCVec2d__
#define __threecpp__TCVec2d__

#include "TCVec2f.h"

class TCVec2d
{
public:
    typedef double value_type;
    
    /** Number of vector components. */
    enum { num_components = 2 };
    
    /** Vec member variable. */
    value_type _v[2];
    
    
    /** Constructor that sets all components of the vector to zero */
    TCVec2d() {_v[0]=0.0; _v[1]=0.0;}
    
    inline TCVec2d(const TCVec2f& vec) { _v[0]=vec._v[0]; _v[1]=vec._v[1];}
    
    inline operator TCVec2f() const { return TCVec2f(static_cast<float>(_v[0]),static_cast<float>(_v[1]));}
    
    TCVec2d(value_type x,value_type y) { _v[0]=x; _v[1]=y; }
    
    
    inline bool operator == (const TCVec2d& v) const { return _v[0]==v._v[0] && _v[1]==v._v[1]; }
    
    inline bool operator != (const TCVec2d& v) const { return _v[0]!=v._v[0] || _v[1]!=v._v[1]; }
    
    inline bool operator <  (const TCVec2d& v) const
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
    inline value_type operator * (const TCVec2d& rhs) const
    {
        return _v[0]*rhs._v[0]+_v[1]*rhs._v[1];
    }
    
    /** Multiply by scalar. */
    inline const TCVec2d operator * (value_type rhs) const
    {
        return TCVec2d(_v[0]*rhs, _v[1]*rhs);
    }
    
    /** Unary multiply by scalar. */
    inline TCVec2d& operator *= (value_type rhs)
    {
        _v[0]*=rhs;
        _v[1]*=rhs;
        return *this;
    }
    
    /** Divide by scalar. */
    inline const TCVec2d operator / (value_type rhs) const
    {
        return TCVec2d(_v[0]/rhs, _v[1]/rhs);
    }
    
    /** Unary divide by scalar. */
    inline TCVec2d& operator /= (value_type rhs)
    {
        _v[0]/=rhs;
        _v[1]/=rhs;
        return *this;
    }
    
    /** Binary vector add. */
    inline const TCVec2d operator + (const TCVec2d& rhs) const
    {
        return TCVec2d(_v[0]+rhs._v[0], _v[1]+rhs._v[1]);
    }
    
    /** Unary vector add. Slightly more efficient because no temporary
     * intermediate object.
     */
    inline TCVec2d& operator += (const TCVec2d& rhs)
    {
        _v[0] += rhs._v[0];
        _v[1] += rhs._v[1];
        return *this;
    }
    
    /** Binary vector subtract. */
    inline const TCVec2d operator - (const TCVec2d& rhs) const
    {
        return TCVec2d(_v[0]-rhs._v[0], _v[1]-rhs._v[1]);
    }
    
    /** Unary vector subtract. */
    inline TCVec2d& operator -= (const TCVec2d& rhs)
    {
        _v[0]-=rhs._v[0];
        _v[1]-=rhs._v[1];
        return *this;
    }
    
    /** Negation operator. Returns the negative of the TCVec2d. */
    inline const TCVec2d operator - () const
    {
        return TCVec2d (-_v[0], -_v[1]);
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
        value_type norm = TCVec2d::length();
        if (norm>0.0)
        {
            value_type inv = 1.0f/norm;
            _v[0] *= inv;
            _v[1] *= inv;
        }
        return( norm );
    }
    
};    // end of class TCVec2d

/** multiply by vector components. */
inline TCVec2d componentMultiply(const TCVec2d& lhs, const TCVec2d& rhs)
{
    return TCVec2d(lhs[0]*rhs[0], lhs[1]*rhs[1]);
}

/** divide rhs components by rhs vector components. */
inline TCVec2d componentDivide(const TCVec2d& lhs, const TCVec2d& rhs)
{
    return TCVec2d(lhs[0]/rhs[0], lhs[1]/rhs[1]);
}

#endif /* defined(__threecpp__TCVec2d__) */
