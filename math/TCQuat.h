//
//  TCQuat.h
//  threecpp
//
//  Created by CastingJ on 16/9/21.
//
//

#ifndef __threecpp__TCQuat__
#define __threecpp__TCQuat__

#include "TCVec4f.h"
#include "TCVec4d.h"

class TCMatrixf;
class TCMatrixd;

/** A TCQuaternion class. It can be used to represent an orientation in 3D space.*/
class TCQuat
{
    
public:
    
    typedef double value_type;
    
    value_type  _v[4];    // a four-vector
    
    inline TCQuat() { _v[0]=0.0; _v[1]=0.0; _v[2]=0.0; _v[3]=1.0; }
    
    inline TCQuat( value_type x, value_type y, value_type z, value_type w )
    {
        _v[0]=x;
        _v[1]=y;
        _v[2]=z;
        _v[3]=w;
    }
    
    inline TCQuat( const TCVec4f& v )
    {
        _v[0]=v.x();
        _v[1]=v.y();
        _v[2]=v.z();
        _v[3]=v.w();
    }
    
    inline TCQuat( const TCVec4d& v )
    {
        _v[0]=v.x();
        _v[1]=v.y();
        _v[2]=v.z();
        _v[3]=v.w();
    }
    
    inline TCQuat( value_type angle, const TCVec3f& axis)
    {
        makeRotate(angle,axis);
    }
    inline TCQuat( value_type angle, const TCVec3d& axis)
    {
        makeRotate(angle,axis);
    }
    
    inline TCQuat( value_type angle1, const TCVec3f& axis1,
                value_type angle2, const TCVec3f& axis2,
                value_type angle3, const TCVec3f& axis3)
    {
        makeRotate(angle1,axis1,angle2,axis2,angle3,axis3);
    }
    
    inline TCQuat( value_type angle1, const TCVec3d& axis1,
                value_type angle2, const TCVec3d& axis2,
                value_type angle3, const TCVec3d& axis3)
    {
        makeRotate(angle1,axis1,angle2,axis2,angle3,axis3);
    }
    
    inline TCQuat& operator = (const TCQuat& v) { _v[0]=v._v[0];  _v[1]=v._v[1]; _v[2]=v._v[2]; _v[3]=v._v[3]; return *this; }
    
    inline bool operator == (const TCQuat& v) const { return _v[0]==v._v[0] && _v[1]==v._v[1] && _v[2]==v._v[2] && _v[3]==v._v[3]; }
    
    inline bool operator != (const TCQuat& v) const { return _v[0]!=v._v[0] || _v[1]!=v._v[1] || _v[2]!=v._v[2] || _v[3]!=v._v[3]; }
    
    inline bool operator <  (const TCQuat& v) const
    {
        if (_v[0]<v._v[0]) return true;
        else if (_v[0]>v._v[0]) return false;
        else if (_v[1]<v._v[1]) return true;
        else if (_v[1]>v._v[1]) return false;
        else if (_v[2]<v._v[2]) return true;
        else if (_v[2]>v._v[2]) return false;
        else return (_v[3]<v._v[3]);
    }
    
    /* ----------------------------------
     Methods to access data members
     ---------------------------------- */
    
    inline TCVec4d asVec4() const
    {
        return TCVec4d(_v[0], _v[1], _v[2], _v[3]);
    }
    
    inline TCVec3d asVec3() const
    {
        return TCVec3d(_v[0], _v[1], _v[2]);
    }
    
    inline void set(value_type x, value_type y, value_type z, value_type w)
    {
        _v[0]=x;
        _v[1]=y;
        _v[2]=z;
        _v[3]=w;
    }
    
    inline void set(const TCVec4f& v)
    {
        _v[0]=v.x();
        _v[1]=v.y();
        _v[2]=v.z();
        _v[3]=v.w();
    }
    
    inline void set(const TCVec4d& v)
    {
        _v[0]=v.x();
        _v[1]=v.y();
        _v[2]=v.z();
        _v[3]=v.w();
    }
    
    void set(const TCMatrixf& matrix);
    
    void set(const TCMatrixd& matrix);
    
    void get(TCMatrixf& matrix) const;
    
    void get(TCMatrixd& matrix) const;
    
    
    inline value_type & operator [] (int i) { return _v[i]; }
    inline value_type   operator [] (int i) const { return _v[i]; }
    
    inline value_type & x() { return _v[0]; }
    inline value_type & y() { return _v[1]; }
    inline value_type & z() { return _v[2]; }
    inline value_type & w() { return _v[3]; }
    
    inline value_type x() const { return _v[0]; }
    inline value_type y() const { return _v[1]; }
    inline value_type z() const { return _v[2]; }
    inline value_type w() const { return _v[3]; }
    
    /** return true if the TCQuat represents a zero rotation, and therefore can be ignored in computations.*/
    bool zeroRotation() const { return _v[0]==0.0 && _v[1]==0.0 && _v[2]==0.0 && _v[3]==1.0; }
    
    
    /* -------------------------------------------------------------
     BASIC ARITHMETIC METHODS
     Implemented in terms of Vec4s.  Some Vec4 operators, e.g.
     operator* are not appropriate for TCQuaternions (as
     mathematical objects) so they are implemented differently.
     Also define methods for conjugate and the multiplicative inverse.
     ------------------------------------------------------------- */
    /// Multiply by scalar
    inline const TCQuat operator * (value_type rhs) const
    {
        return TCQuat(_v[0]*rhs, _v[1]*rhs, _v[2]*rhs, _v[3]*rhs);
    }
    
    /// Unary multiply by scalar
    inline TCQuat& operator *= (value_type rhs)
    {
        _v[0]*=rhs;
        _v[1]*=rhs;
        _v[2]*=rhs;
        _v[3]*=rhs;
        return *this;        // enable nesting
    }
    
    /// Binary multiply
    inline const TCQuat operator*(const TCQuat& rhs) const
    {
        return TCQuat( rhs._v[3]*_v[0] + rhs._v[0]*_v[3] + rhs._v[1]*_v[2] - rhs._v[2]*_v[1],
                    rhs._v[3]*_v[1] - rhs._v[0]*_v[2] + rhs._v[1]*_v[3] + rhs._v[2]*_v[0],
                    rhs._v[3]*_v[2] + rhs._v[0]*_v[1] - rhs._v[1]*_v[0] + rhs._v[2]*_v[3],
                    rhs._v[3]*_v[3] - rhs._v[0]*_v[0] - rhs._v[1]*_v[1] - rhs._v[2]*_v[2] );
    }
    
    /// Unary multiply
    inline TCQuat& operator*=(const TCQuat& rhs)
    {
        value_type x = rhs._v[3]*_v[0] + rhs._v[0]*_v[3] + rhs._v[1]*_v[2] - rhs._v[2]*_v[1];
        value_type y = rhs._v[3]*_v[1] - rhs._v[0]*_v[2] + rhs._v[1]*_v[3] + rhs._v[2]*_v[0];
        value_type z = rhs._v[3]*_v[2] + rhs._v[0]*_v[1] - rhs._v[1]*_v[0] + rhs._v[2]*_v[3];
        _v[3]   = rhs._v[3]*_v[3] - rhs._v[0]*_v[0] - rhs._v[1]*_v[1] - rhs._v[2]*_v[2];
        
        _v[2] = z;
        _v[1] = y;
        _v[0] = x;
        
        return (*this);            // enable nesting
    }
    
    /// Divide by scalar
    inline TCQuat operator / (value_type rhs) const
    {
        value_type div = 1.0/rhs;
        return TCQuat(_v[0]*div, _v[1]*div, _v[2]*div, _v[3]*div);
    }
    
    /// Unary divide by scalar
    inline TCQuat& operator /= (value_type rhs)
    {
        value_type div = 1.0/rhs;
        _v[0]*=div;
        _v[1]*=div;
        _v[2]*=div;
        _v[3]*=div;
        return *this;
    }
    
    /// Binary divide
    inline const TCQuat operator/(const TCQuat& denom) const
    {
        return ( (*this) * denom.inverse() );
    }
    
    /// Unary divide
    inline TCQuat& operator/=(const TCQuat& denom)
    {
        (*this) = (*this) * denom.inverse();
        return (*this);            // enable nesting
    }
    
    /// Binary addition
    inline const TCQuat operator + (const TCQuat& rhs) const
    {
        return TCQuat(_v[0]+rhs._v[0], _v[1]+rhs._v[1],
                    _v[2]+rhs._v[2], _v[3]+rhs._v[3]);
    }
    
    /// Unary addition
    inline TCQuat& operator += (const TCQuat& rhs)
    {
        _v[0] += rhs._v[0];
        _v[1] += rhs._v[1];
        _v[2] += rhs._v[2];
        _v[3] += rhs._v[3];
        return *this;            // enable nesting
    }
    
    /// Binary subtraction
    inline const TCQuat operator - (const TCQuat& rhs) const
    {
        return TCQuat(_v[0]-rhs._v[0], _v[1]-rhs._v[1],
                    _v[2]-rhs._v[2], _v[3]-rhs._v[3] );
    }
    
    /// Unary subtraction
    inline TCQuat& operator -= (const TCQuat& rhs)
    {
        _v[0]-=rhs._v[0];
        _v[1]-=rhs._v[1];
        _v[2]-=rhs._v[2];
        _v[3]-=rhs._v[3];
        return *this;            // enable nesting
    }
    
    /** Negation operator - returns the negative of the TCQuaternion.
     Basically just calls operator - () on the Vec4 */
    inline const TCQuat operator - () const
    {
        return TCQuat (-_v[0], -_v[1], -_v[2], -_v[3]);
    }
    
    /// Length of the TCQuaternion = sqrt( vec . vec )
    value_type length() const
    {
        return sqrt( _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3]);
    }
    
    /// Length of the TCQuaternion = vec . vec
    value_type length2() const
    {
        return _v[0]*_v[0] + _v[1]*_v[1] + _v[2]*_v[2] + _v[3]*_v[3];
    }
    
    /// Conjugate
    inline TCQuat conj () const
    {
        return TCQuat( -_v[0], -_v[1], -_v[2], _v[3] );
    }
    
    /// Multiplicative inverse method: q^(-1) = q^*/(q.q^*)
    inline const TCQuat inverse () const
    {
        return conj() / length2();
    }
    
    /* --------------------------------------------------------
     METHODS RELATED TO ROTATIONS
     Set a TCQuaternion which will perform a rotation of an
     angle around the axis given by the vector (x,y,z).
     Should be written to also accept an angle and a Vec3?
     
     Define Spherical Linear interpolation method also
     
     Not inlined - see the TCQuat.cpp file for implementation
     -------------------------------------------------------- */
    void makeRotate( value_type  angle,
                    value_type  x, value_type  y, value_type  z );
    void makeRotate ( value_type  angle, const TCVec3f& vec );
    void makeRotate ( value_type  angle, const TCVec3d& vec );
    
    void makeRotate ( value_type  angle1, const TCVec3f& axis1,
                     value_type  angle2, const TCVec3f& axis2,
                     value_type  angle3, const TCVec3f& axis3);
    void makeRotate ( value_type  angle1, const TCVec3d& axis1,
                     value_type  angle2, const TCVec3d& axis2,
                     value_type  angle3, const TCVec3d& axis3);
    
    /** Make a rotation TCQuat which will rotate vec1 to vec2.
     Generally take a dot product to get the angle between these
     and then use a cross product to get the rotation axis
     Watch out for the two special cases when the vectors
     are co-incident or opposite in direction.*/
    void makeRotate( const TCVec3f& vec1, const TCVec3f& vec2 );
    /** Make a rotation TCQuat which will rotate vec1 to vec2.
     Generally take a dot product to get the angle between these
     and then use a cross product to get the rotation axis
     Watch out for the two special cases of when the vectors
     are co-incident or opposite in direction.*/
    void makeRotate( const TCVec3d& vec1, const TCVec3d& vec2 );
    
    void makeRotate_original( const TCVec3d& vec1, const TCVec3d& vec2 );
    
    /** Return the angle and vector components represented by the TCQuaternion.*/
    void getRotate ( value_type & angle, value_type & x, value_type & y, value_type & z ) const;
    
    /** Return the angle and vector represented by the TCQuaternion.*/
    void getRotate ( value_type & angle, TCVec3f& vec ) const;
    
    /** Return the angle and vector represented by the TCQuaternion.*/
    void getRotate ( value_type & angle, TCVec3d& vec ) const;
    
    /** Spherical Linear Interpolation.
     As t goes from 0 to 1, the TCQuat object goes from "from" to "to". */
    void slerp   ( value_type  t, const TCQuat& from, const TCQuat& to);
    
    /** Rotate a vector by this TCQuaternion.*/
    TCVec3f operator* (const TCVec3f& v) const
    {
        // nVidia SDK implementation
        TCVec3f uv, uuv;
        TCVec3f qvec(_v[0], _v[1], _v[2]);
        uv = qvec ^ v;
        uuv = qvec ^ uv;
        uv *= ( 2.0f * _v[3] );
        uuv *= 2.0f;
        return v + uv + uuv;
    }
    
    /** Rotate a vector by this TCQuaternion.*/
    TCVec3d operator* (const TCVec3d& v) const
    {
        // nVidia SDK implementation
        TCVec3d uv, uuv;
        TCVec3d qvec(_v[0], _v[1], _v[2]);
        uv = qvec ^ v;
        uuv = qvec ^ uv;
        uv *= ( 2.0f * _v[3] );
        uuv *= 2.0f;
        return v + uv + uuv;
    }
    
protected:
    
};    // end of class prototype


#endif /* defined(__threecpp__TCQuat__) */
