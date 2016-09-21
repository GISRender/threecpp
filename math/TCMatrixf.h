//
//  TCMatrixf.h
//  threecpp
//
//  Created by CastingJ on 16/9/21.
//
//

#ifndef __threecpp__TCMatrixf__
#define __threecpp__TCMatrixf__

#include "TCVec3d.h"
#include "TCVec4d.h"
#include "TCQuat.h"

class TCMatrixd;

class TCMatrixf
{
public:
    
    typedef float value_type;
    typedef double other_value_type;
    
    inline TCMatrixf() { makeIdentity(); }
    inline TCMatrixf( const TCMatrixf& mat) { set(mat.ptr()); }
    TCMatrixf( const TCMatrixd& mat );
    inline explicit TCMatrixf( float const * const ptr ) { set(ptr); }
    inline explicit TCMatrixf( double const * const ptr ) { set(ptr); }
    inline explicit TCMatrixf( const TCQuat& TCQuat ) { makeRotate(TCQuat); }
    
    TCMatrixf( value_type a00, value_type a01, value_type a02, value_type a03,
            value_type a10, value_type a11, value_type a12, value_type a13,
            value_type a20, value_type a21, value_type a22, value_type a23,
            value_type a30, value_type a31, value_type a32, value_type a33);
    
    ~TCMatrixf() {}
    
    int compare(const TCMatrixf& m) const;
    
    bool operator < (const TCMatrixf& m) const { return compare(m)<0; }
    bool operator == (const TCMatrixf& m) const { return compare(m)==0; }
    bool operator != (const TCMatrixf& m) const { return compare(m)!=0; }
    
    inline value_type& operator()(int row, int col) { return _mat[row][col]; }
    inline value_type operator()(int row, int col) const { return _mat[row][col]; }
    
    inline bool valid() const
    {
        return !(isNaN(_mat[0][0]) || isNaN(_mat[0][1]) || isNaN(_mat[0][2]) || isNaN(_mat[0][3]) ||
                 isNaN(_mat[1][0]) || isNaN(_mat[1][1]) || isNaN(_mat[1][2]) || isNaN(_mat[1][3]) ||
                 isNaN(_mat[2][0]) || isNaN(_mat[2][1]) || isNaN(_mat[2][2]) || isNaN(_mat[2][3]) ||
                 isNaN(_mat[3][0]) || isNaN(_mat[3][1]) || isNaN(_mat[3][2]) || isNaN(_mat[3][3]));
    }
    
    
    
    inline TCMatrixf& operator = (const TCMatrixf& rhs)
    {
        if( &rhs == this ) return *this;
        set(rhs.ptr());
        return *this;
    }
    
    TCMatrixf& operator = (const TCMatrixd& other);
    
    inline void set(const TCMatrixf& rhs) { set(rhs.ptr()); }
    
    void set(const TCMatrixd& rhs);
    
    inline void set(float const * const ptr)
    {
        value_type* local_ptr = (value_type*)_mat;
        for(int i=0;i<16;++i) local_ptr[i]=(value_type)ptr[i];
    }
    
    inline void set(double const * const ptr)
    {
        value_type* local_ptr = (value_type*)_mat;
        for(int i=0;i<16;++i) local_ptr[i]=(value_type)ptr[i];
    }
    
    void set(value_type a00, value_type a01, value_type a02,value_type a03,
             value_type a10, value_type a11, value_type a12,value_type a13,
             value_type a20, value_type a21, value_type a22,value_type a23,
             value_type a30, value_type a31, value_type a32,value_type a33);
    
    value_type * ptr() { return (value_type*)_mat; }
    const value_type * ptr() const { return (const value_type *)_mat; }
    
    bool isIdentity() const
    {
        return _mat[0][0]==1.0f && _mat[0][1]==0.0f && _mat[0][2]==0.0f &&  _mat[0][3]==0.0f &&
        _mat[1][0]==0.0f && _mat[1][1]==1.0f && _mat[1][2]==0.0f &&  _mat[1][3]==0.0f &&
        _mat[2][0]==0.0f && _mat[2][1]==0.0f && _mat[2][2]==1.0f &&  _mat[2][3]==0.0f &&
        _mat[3][0]==0.0f && _mat[3][1]==0.0f && _mat[3][2]==0.0f &&  _mat[3][3]==1.0f;
    }
    
    void makeIdentity();
    
    void makeScale( const TCVec3f& );
    void makeScale( const TCVec3d& );
    void makeScale( value_type, value_type, value_type );
    
    void makeTranslate( const TCVec3f& );
    void makeTranslate( const TCVec3d& );
    void makeTranslate( value_type, value_type, value_type );
    
    void makeRotate( const TCVec3f& from, const TCVec3f& to );
    void makeRotate( const TCVec3d& from, const TCVec3d& to );
    void makeRotate( value_type angle, const TCVec3f& axis );
    void makeRotate( value_type angle, const TCVec3d& axis );
    void makeRotate( value_type angle, value_type x, value_type y, value_type z );
    void makeRotate( const TCQuat& );
    void makeRotate( value_type angle1, const TCVec3f& axis1,
                    value_type angle2, const TCVec3f& axis2,
                    value_type angle3, const TCVec3f& axis3);
    void makeRotate( value_type angle1, const TCVec3d& axis1,
                    value_type angle2, const TCVec3d& axis2,
                    value_type angle3, const TCVec3d& axis3);
    
    
    /** decompose the TCMatrix into translation, rotation, scale and scale orientation.*/
    void decompose( TCVec3f& translation,
                   TCQuat& rotation,
                   TCVec3f& scale,
                   TCQuat& so ) const;
    
    /** decompose the TCMatrix into translation, rotation, scale and scale orientation.*/
    void decompose( TCVec3d& translation,
                   TCQuat& rotation,
                   TCVec3d& scale,
                   TCQuat& so ) const;
    
    
    /** Set to an orthographic projection.
     * See glOrtho for further details.
     */
    void makeOrtho(double left,   double right,
                   double bottom, double top,
                   double zNear,  double zFar);
    
    /** Get the orthographic settings of the orthographic projection TCMatrix.
     * Note, if TCMatrix is not an orthographic TCMatrix then invalid values
     * will be returned.
     */
    bool getOrtho(double& left,   double& right,
                  double& bottom, double& top,
                  double& zNear,  double& zFar) const;
    
    /** float version of getOrtho(..) */
    bool getOrtho(float& left,   float& right,
                  float& bottom, float& top,
                  float& zNear,  float& zFar) const;
    
    
    /** Set to a 2D orthographic projection.
     * See glOrtho2D for further details.
     */
    inline void makeOrtho2D(double left,   double right,
                            double bottom, double top)
    {
        makeOrtho(left,right,bottom,top,-1.0,1.0);
    }
    
    
    /** Set to a perspective projection.
     * See glFrustum for further details.
     */
    void makeFrustum(double left,   double right,
                     double bottom, double top,
                     double zNear,  double zFar);
    
    /** Get the frustum settings of a perspective projection TCMatrix.
     * Note, if TCMatrix is not a perspective TCMatrix then invalid values
     * will be returned.
     */
    bool getFrustum(double& left,   double& right,
                    double& bottom, double& top,
                    double& zNear,  double& zFar) const;
    
    /** float version of getFrustum(..) */
    bool getFrustum(float& left,   float& right,
                    float& bottom, float& top,
                    float& zNear,  float& zFar) const;
    
    /** Set to a symmetrical perspective projection.
     * See gluPerspective for further details.
     * Aspect ratio is defined as width/height.
     */
    void makePerspective(double fovy,  double aspectRatio,
                         double zNear, double zFar);
    
    /** Get the frustum settings of a symmetric perspective projection
     * TCMatrix.
     * Return false if TCMatrix is not a perspective TCMatrix,
     * where parameter values are undefined.
     * Note, if TCMatrix is not a symmetric perspective TCMatrix then the
     * shear will be lost.
     * Asymmetric matrices occur when stereo, power walls, caves and
     * reality center display are used.
     * In these configuration one should use the AsFrustum method instead.
     */
    bool getPerspective(double& fovy,  double& aspectRatio,
                        double& zNear, double& zFar) const;
    
    /** float version of getPerspective(..) */
    bool getPerspective(float& fovy,  float& aspectRatio,
                        float& zNear, float& zFar) const;
    
    /** Set the position and orientation to be a view TCMatrix,
     * using the same convention as gluLookAt.
     */
    void makeLookAt(const TCVec3d& eye,const TCVec3d& center,const TCVec3d& up);
    
    /** Get to the position and orientation of a modelview TCMatrix,
     * using the same convention as gluLookAt.
     */
    void getLookAt(TCVec3f& eye,TCVec3f& center,TCVec3f& up,
                   value_type lookDistance=1.0f) const;
    
    /** Get to the position and orientation of a modelview TCMatrix,
     * using the same convention as gluLookAt.
     */
    void getLookAt(TCVec3d& eye,TCVec3d& center,TCVec3d& up,
                   value_type lookDistance=1.0f) const;
    
    /** invert the TCMatrix rhs, automatically select invert_4x3 or invert_4x4. */
    inline bool invert( const TCMatrixf& rhs)
    {
        bool is_4x3 = (rhs._mat[0][3]==0.0f && rhs._mat[1][3]==0.0f &&  rhs._mat[2][3]==0.0f && rhs._mat[3][3]==1.0f);
        return is_4x3 ? invert_4x3(rhs) :  invert_4x4(rhs);
    }
    
    /** 4x3 TCMatrix invert, not right hand column is assumed to be 0,0,0,1. */
    bool invert_4x3( const TCMatrixf& rhs);
    
    /** full 4x4 TCMatrix invert. */
    bool invert_4x4( const TCMatrixf& rhs);
    
    /** ortho-normalize the 3x3 rotation & scale TCMatrix */
    void orthoNormalize(const TCMatrixf& rhs);
    
    //basic utility functions to create new matrices
    inline static TCMatrixf identity( void );
    inline static TCMatrixf scale( const TCVec3f& sv);
    inline static TCMatrixf scale( const TCVec3d& sv);
    inline static TCMatrixf scale( value_type sx, value_type sy, value_type sz);
    inline static TCMatrixf translate( const TCVec3f& dv);
    inline static TCMatrixf translate( const TCVec3d& dv);
    inline static TCMatrixf translate( value_type x, value_type y, value_type z);
    inline static TCMatrixf rotate( const TCVec3f& from, const TCVec3f& to);
    inline static TCMatrixf rotate( const TCVec3d& from, const TCVec3d& to);
    inline static TCMatrixf rotate( value_type angle, value_type x, value_type y, value_type z);
    inline static TCMatrixf rotate( value_type angle, const TCVec3f& axis);
    inline static TCMatrixf rotate( value_type angle, const TCVec3d& axis);
    inline static TCMatrixf rotate( value_type angle1, const TCVec3f& axis1,
                                 value_type angle2, const TCVec3f& axis2,
                                 value_type angle3, const TCVec3f& axis3);
    inline static TCMatrixf rotate( value_type angle1, const TCVec3d& axis1,
                                 value_type angle2, const TCVec3d& axis2,
                                 value_type angle3, const TCVec3d& axis3);
    inline static TCMatrixf rotate( const TCQuat& TCQuat);
    inline static TCMatrixf inverse( const TCMatrixf& TCMatrix);
    inline static TCMatrixf orthoNormal(const TCMatrixf& TCMatrix);
    
    /** Create an orthographic projection TCMatrix.
     * See glOrtho for further details.
     */
    inline static TCMatrixf ortho(double left,   double right,
                                double bottom, double top,
                                double zNear,  double zFar);
    
    /** Create a 2D orthographic projection.
     * See glOrtho for further details.
     */
    inline static TCMatrixf ortho2D(double left,   double right,
                                  double bottom, double top);
    
    /** Create a perspective projection.
     * See glFrustum for further details.
     */
    inline static TCMatrixf frustum(double left,   double right,
                                  double bottom, double top,
                                  double zNear,  double zFar);
    
    /** Create a symmetrical perspective projection.
     * See gluPerspective for further details.
     * Aspect ratio is defined as width/height.
     */
    inline static TCMatrixf perspective(double fovy,  double aspectRatio,
                                      double zNear, double zFar);
    
    /** Create the position and orientation as per a camera,
     * using the same convention as gluLookAt.
     */
    inline static TCMatrixf lookAt(const TCVec3f& eye,
                                 const TCVec3f& center,
                                 const TCVec3f& up);
    
    /** Create the position and orientation as per a camera,
     * using the same convention as gluLookAt.
     */
    inline static TCMatrixf lookAt(const TCVec3d& eye,
                                 const TCVec3d& center,
                                 const TCVec3d& up);
    
    inline TCVec3f preMult( const TCVec3f& v ) const;
    inline TCVec3d preMult( const TCVec3d& v ) const;
    inline TCVec3f postMult( const TCVec3f& v ) const;
    inline TCVec3d postMult( const TCVec3d& v ) const;
    inline TCVec3f operator* ( const TCVec3f& v ) const;
    inline TCVec3d operator* ( const TCVec3d& v ) const;
    inline TCVec4f preMult( const TCVec4f& v ) const;
    inline TCVec4d preMult( const TCVec4d& v ) const;
    inline TCVec4f postMult( const TCVec4f& v ) const;
    inline TCVec4d postMult( const TCVec4d& v ) const;
    inline TCVec4f operator* ( const TCVec4f& v ) const;
    inline TCVec4d operator* ( const TCVec4d& v ) const;
    
#ifdef USE_DEPRECATED_API
    inline void set(const TCQuat& q) { makeRotate(q); }
    inline void get(TCQuat& q) const { q = getRotate(); }
#endif
    
    void setRotate(const TCQuat& q);
    /** Get the TCMatrix rotation as a TCQuat. Note that this function
     * assumes a non-scaled TCMatrix and will return incorrect results
     * for scaled TCMatrixces. Consider decompose() instead.
     */
    TCQuat getRotate() const;
    
    
    void setTrans( value_type tx, value_type ty, value_type tz );
    void setTrans( const TCVec3f& v );
    void setTrans( const TCVec3d& v );
    
    inline TCVec3d getTrans() const { return TCVec3d(_mat[3][0],_mat[3][1],_mat[3][2]); }
    
    inline TCVec3d getScale() const {
        TCVec3d x_vec(_mat[0][0],_mat[1][0],_mat[2][0]);
        TCVec3d y_vec(_mat[0][1],_mat[1][1],_mat[2][1]);
        TCVec3d z_vec(_mat[0][2],_mat[1][2],_mat[2][2]);
        return TCVec3d(x_vec.length(), y_vec.length(), z_vec.length());
    }
    
    /** apply a 3x3 transform of v*M[0..2,0..2]. */
    inline static TCVec3f transform3x3(const TCVec3f& v,const TCMatrixf& m);
    
    /** apply a 3x3 transform of v*M[0..2,0..2]. */
    inline static TCVec3d transform3x3(const TCVec3d& v,const TCMatrixf& m);
    
    /** apply a 3x3 transform of M[0..2,0..2]*v. */
    inline static TCVec3f transform3x3(const TCMatrixf& m,const TCVec3f& v);
    
    /** apply a 3x3 transform of M[0..2,0..2]*v. */
    inline static TCVec3d transform3x3(const TCMatrixf& m,const TCVec3d& v);
    
    // basic TCMatrixf multiplication, our workhorse methods.
    void mult( const TCMatrixf&, const TCMatrixf& );
    void preMult( const TCMatrixf& );
    void postMult( const TCMatrixf& );
    
    /** Optimized version of preMult(translate(v)); */
    inline void preMultTranslate( const TCVec3d& v );
    inline void preMultTranslate( const TCVec3f& v );
    /** Optimized version of postMult(translate(v)); */
    inline void postMultTranslate( const TCVec3d& v );
    inline void postMultTranslate( const TCVec3f& v );
    
    /** Optimized version of preMult(scale(v)); */
    inline void preMultScale( const TCVec3d& v );
    inline void preMultScale( const TCVec3f& v );
    /** Optimized version of postMult(scale(v)); */
    inline void postMultScale( const TCVec3d& v );
    inline void postMultScale( const TCVec3f& v );
    
    /** Optimized version of preMult(rotate(q)); */
    inline void preMultRotate( const TCQuat& q );
    /** Optimized version of postMult(rotate(q)); */
    inline void postMultRotate( const TCQuat& q );
    
    inline void operator *= ( const TCMatrixf& other )
    {    if( this == &other ) {
        TCMatrixf temp(other);
        postMult( temp );
    }
    else postMult( other );
    }
    
    inline TCMatrixf operator * ( const TCMatrixf &m ) const
    {
        TCMatrixf r;
        r.mult(*this,m);
        return  r;
    }
    
    /** Multiply by scalar. */
    inline TCMatrixf operator * (value_type rhs) const
    {
        return TCMatrixf(
                       _mat[0][0]*rhs, _mat[0][1]*rhs, _mat[0][2]*rhs, _mat[0][3]*rhs,
                       _mat[1][0]*rhs, _mat[1][1]*rhs, _mat[1][2]*rhs, _mat[1][3]*rhs,
                       _mat[2][0]*rhs, _mat[2][1]*rhs, _mat[2][2]*rhs, _mat[2][3]*rhs,
                       _mat[3][0]*rhs, _mat[3][1]*rhs, _mat[3][2]*rhs, _mat[3][3]*rhs);
    }
    
    /** Unary multiply by scalar. */
    inline TCMatrixf& operator *= (value_type rhs)
    {
        _mat[0][0]*=rhs;
        _mat[0][1]*=rhs;
        _mat[0][2]*=rhs;
        _mat[0][3]*=rhs;
        _mat[1][0]*=rhs;
        _mat[1][1]*=rhs;
        _mat[1][2]*=rhs;
        _mat[1][3]*=rhs;
        _mat[2][0]*=rhs;
        _mat[2][1]*=rhs;
        _mat[2][2]*=rhs;
        _mat[2][3]*=rhs;
        _mat[3][0]*=rhs;
        _mat[3][1]*=rhs;
        _mat[3][2]*=rhs;
        _mat[3][3]*=rhs;
        return *this;
    }
    
    /** Divide by scalar. */
    inline TCMatrixf operator / (value_type rhs) const
    {
        return TCMatrixf(
                       _mat[0][0]/rhs, _mat[0][1]/rhs, _mat[0][2]/rhs, _mat[0][3]/rhs,
                       _mat[1][0]/rhs, _mat[1][1]/rhs, _mat[1][2]/rhs, _mat[1][3]/rhs,
                       _mat[2][0]/rhs, _mat[2][1]/rhs, _mat[2][2]/rhs, _mat[2][3]/rhs,
                       _mat[3][0]/rhs, _mat[3][1]/rhs, _mat[3][2]/rhs, _mat[3][3]/rhs);
    }
    
    /** Unary divide by scalar. */
    inline TCMatrixf& operator /= (value_type rhs)
    {
        _mat[0][0]/=rhs;
        _mat[0][1]/=rhs;
        _mat[0][2]/=rhs;
        _mat[0][3]/=rhs;
        _mat[1][0]/=rhs;
        _mat[1][1]/=rhs;
        _mat[1][2]/=rhs;
        _mat[1][3]/=rhs;
        _mat[2][0]/=rhs;
        _mat[2][1]/=rhs;
        _mat[2][2]/=rhs;
        _mat[2][3]/=rhs;
        _mat[3][0]/=rhs;
        _mat[3][1]/=rhs;
        _mat[3][2]/=rhs;
        _mat[3][3]/=rhs;
        return *this;
    }
    
    /** Binary vector add. */
    inline TCMatrixf operator + (const TCMatrixf& rhs) const
    {
        return TCMatrixf(
                       _mat[0][0] + rhs._mat[0][0],
                       _mat[0][1] + rhs._mat[0][1],
                       _mat[0][2] + rhs._mat[0][2],
                       _mat[0][3] + rhs._mat[0][3],
                       _mat[1][0] + rhs._mat[1][0],
                       _mat[1][1] + rhs._mat[1][1],
                       _mat[1][2] + rhs._mat[1][2],
                       _mat[1][3] + rhs._mat[1][3],
                       _mat[2][0] + rhs._mat[2][0],
                       _mat[2][1] + rhs._mat[2][1],
                       _mat[2][2] + rhs._mat[2][2],
                       _mat[2][3] + rhs._mat[2][3],
                       _mat[3][0] + rhs._mat[3][0],
                       _mat[3][1] + rhs._mat[3][1],
                       _mat[3][2] + rhs._mat[3][2],
                       _mat[3][3] + rhs._mat[3][3]);
    }
    
    /** Unary vector add. Slightly more efficient because no temporary
     * intermediate object.
     */
    inline TCMatrixf& operator += (const TCMatrixf& rhs)
    {
        _mat[0][0] += rhs._mat[0][0];
        _mat[0][1] += rhs._mat[0][1];
        _mat[0][2] += rhs._mat[0][2];
        _mat[0][3] += rhs._mat[0][3];
        _mat[1][0] += rhs._mat[1][0];
        _mat[1][1] += rhs._mat[1][1];
        _mat[1][2] += rhs._mat[1][2];
        _mat[1][3] += rhs._mat[1][3];
        _mat[2][0] += rhs._mat[2][0];
        _mat[2][1] += rhs._mat[2][1];
        _mat[2][2] += rhs._mat[2][2];
        _mat[2][3] += rhs._mat[2][3];
        _mat[3][0] += rhs._mat[3][0];
        _mat[3][1] += rhs._mat[3][1];
        _mat[3][2] += rhs._mat[3][2];
        _mat[3][3] += rhs._mat[3][3];
        return *this;
    }
    
protected:
    value_type _mat[4][4];
    
};


//static utility methods
inline TCMatrixf TCMatrixf::identity(void)
{
    TCMatrixf m;
    m.makeIdentity();
    return m;
}

inline TCMatrixf TCMatrixf::scale(value_type sx, value_type sy, value_type sz)
{
    TCMatrixf m;
    m.makeScale(sx,sy,sz);
    return m;
}

inline TCMatrixf TCMatrixf::scale(const TCVec3f& v )
{
    return scale(v.x(), v.y(), v.z() );
}

inline TCMatrixf TCMatrixf::scale(const TCVec3d& v )
{
    return scale(v.x(), v.y(), v.z() );
}

inline TCMatrixf TCMatrixf::translate(value_type tx, value_type ty, value_type tz)
{
    TCMatrixf m;
    m.makeTranslate(tx,ty,tz);
    return m;
}

inline TCMatrixf TCMatrixf::translate(const TCVec3f& v )
{
    return translate(v.x(), v.y(), v.z() );
}

inline TCMatrixf TCMatrixf::translate(const TCVec3d& v )
{
    return translate(v.x(), v.y(), v.z() );
}

inline TCMatrixf TCMatrixf::rotate( const TCQuat& q )
{
    return TCMatrixf(q);
}
inline TCMatrixf TCMatrixf::rotate(value_type angle, value_type x, value_type y, value_type z )
{
    TCMatrixf m;
    m.makeRotate(angle,x,y,z);
    return m;
}
inline TCMatrixf TCMatrixf::rotate(value_type angle, const TCVec3f& axis )
{
    TCMatrixf m;
    m.makeRotate(angle,axis);
    return m;
}
inline TCMatrixf TCMatrixf::rotate(value_type angle, const TCVec3d& axis )
{
    TCMatrixf m;
    m.makeRotate(angle,axis);
    return m;
}
inline TCMatrixf TCMatrixf::rotate( value_type angle1, const TCVec3f& axis1,
                               value_type angle2, const TCVec3f& axis2,
                               value_type angle3, const TCVec3f& axis3)
{
    TCMatrixf m;
    m.makeRotate(angle1,axis1,angle2,axis2,angle3,axis3);
    return m;
}
inline TCMatrixf TCMatrixf::rotate( value_type angle1, const TCVec3d& axis1,
                               value_type angle2, const TCVec3d& axis2,
                               value_type angle3, const TCVec3d& axis3)
{
    TCMatrixf m;
    m.makeRotate(angle1,axis1,angle2,axis2,angle3,axis3);
    return m;
}
inline TCMatrixf TCMatrixf::rotate(const TCVec3f& from, const TCVec3f& to )
{
    TCMatrixf m;
    m.makeRotate(from,to);
    return m;
}
inline TCMatrixf TCMatrixf::rotate(const TCVec3d& from, const TCVec3d& to )
{
    TCMatrixf m;
    m.makeRotate(from,to);
    return m;
}

inline TCMatrixf TCMatrixf::inverse( const TCMatrixf& TCMatrix)
{
    TCMatrixf m;
    m.invert(TCMatrix);
    return m;
}

inline TCMatrixf TCMatrixf::orthoNormal(const TCMatrixf& TCMatrix)
{
    TCMatrixf m;
    m.orthoNormalize(TCMatrix);
    return m;
}

inline TCMatrixf TCMatrixf::ortho(double left, double right,
                              double bottom, double top,
                              double zNear, double zFar)
{
    TCMatrixf m;
    m.makeOrtho(left,right,bottom,top,zNear,zFar);
    return m;
}

inline TCMatrixf TCMatrixf::ortho2D(double left, double right,
                                double bottom, double top)
{
    TCMatrixf m;
    m.makeOrtho2D(left,right,bottom,top);
    return m;
}

inline TCMatrixf TCMatrixf::frustum(double left, double right,
                                double bottom, double top,
                                double zNear, double zFar)
{
    TCMatrixf m;
    m.makeFrustum(left,right,bottom,top,zNear,zFar);
    return m;
}

inline TCMatrixf TCMatrixf::perspective(double fovy,double aspectRatio,
                                    double zNear, double zFar)
{
    TCMatrixf m;
    m.makePerspective(fovy,aspectRatio,zNear,zFar);
    return m;
}

inline TCMatrixf TCMatrixf::lookAt(const TCVec3f& eye,const TCVec3f& center,const TCVec3f& up)
{
    TCMatrixf m;
    m.makeLookAt(eye,center,up);
    return m;
}

inline TCMatrixf TCMatrixf::lookAt(const TCVec3d& eye,const TCVec3d& center,const TCVec3d& up)
{
    TCMatrixf m;
    m.makeLookAt(eye,center,up);
    return m;
}

inline TCVec3f TCMatrixf::postMult( const TCVec3f& v ) const
{
    value_type d = 1.0f/(_mat[3][0]*v.x()+_mat[3][1]*v.y()+_mat[3][2]*v.z()+_mat[3][3]) ;
    return TCVec3f( (_mat[0][0]*v.x() + _mat[0][1]*v.y() + _mat[0][2]*v.z() + _mat[0][3])*d,
                 (_mat[1][0]*v.x() + _mat[1][1]*v.y() + _mat[1][2]*v.z() + _mat[1][3])*d,
                 (_mat[2][0]*v.x() + _mat[2][1]*v.y() + _mat[2][2]*v.z() + _mat[2][3])*d) ;
}
inline TCVec3d TCMatrixf::postMult( const TCVec3d& v ) const
{
    value_type d = 1.0f/(_mat[3][0]*v.x()+_mat[3][1]*v.y()+_mat[3][2]*v.z()+_mat[3][3]) ;
    return TCVec3d( (_mat[0][0]*v.x() + _mat[0][1]*v.y() + _mat[0][2]*v.z() + _mat[0][3])*d,
                 (_mat[1][0]*v.x() + _mat[1][1]*v.y() + _mat[1][2]*v.z() + _mat[1][3])*d,
                 (_mat[2][0]*v.x() + _mat[2][1]*v.y() + _mat[2][2]*v.z() + _mat[2][3])*d) ;
}

inline TCVec3f TCMatrixf::preMult( const TCVec3f& v ) const
{
    value_type d = 1.0f/(_mat[0][3]*v.x()+_mat[1][3]*v.y()+_mat[2][3]*v.z()+_mat[3][3]) ;
    return TCVec3f( (_mat[0][0]*v.x() + _mat[1][0]*v.y() + _mat[2][0]*v.z() + _mat[3][0])*d,
                 (_mat[0][1]*v.x() + _mat[1][1]*v.y() + _mat[2][1]*v.z() + _mat[3][1])*d,
                 (_mat[0][2]*v.x() + _mat[1][2]*v.y() + _mat[2][2]*v.z() + _mat[3][2])*d);
}
inline TCVec3d TCMatrixf::preMult( const TCVec3d& v ) const
{
    value_type d = 1.0f/(_mat[0][3]*v.x()+_mat[1][3]*v.y()+_mat[2][3]*v.z()+_mat[3][3]) ;
    return TCVec3d( (_mat[0][0]*v.x() + _mat[1][0]*v.y() + _mat[2][0]*v.z() + _mat[3][0])*d,
                 (_mat[0][1]*v.x() + _mat[1][1]*v.y() + _mat[2][1]*v.z() + _mat[3][1])*d,
                 (_mat[0][2]*v.x() + _mat[1][2]*v.y() + _mat[2][2]*v.z() + _mat[3][2])*d);
}

inline TCVec4f TCMatrixf::postMult( const TCVec4f& v ) const
{
    return TCVec4f( (_mat[0][0]*v.x() + _mat[0][1]*v.y() + _mat[0][2]*v.z() + _mat[0][3]*v.w()),
                 (_mat[1][0]*v.x() + _mat[1][1]*v.y() + _mat[1][2]*v.z() + _mat[1][3]*v.w()),
                 (_mat[2][0]*v.x() + _mat[2][1]*v.y() + _mat[2][2]*v.z() + _mat[2][3]*v.w()),
                 (_mat[3][0]*v.x() + _mat[3][1]*v.y() + _mat[3][2]*v.z() + _mat[3][3]*v.w())) ;
}
inline TCVec4d TCMatrixf::postMult( const TCVec4d& v ) const
{
    return TCVec4d( (_mat[0][0]*v.x() + _mat[0][1]*v.y() + _mat[0][2]*v.z() + _mat[0][3]*v.w()),
                 (_mat[1][0]*v.x() + _mat[1][1]*v.y() + _mat[1][2]*v.z() + _mat[1][3]*v.w()),
                 (_mat[2][0]*v.x() + _mat[2][1]*v.y() + _mat[2][2]*v.z() + _mat[2][3]*v.w()),
                 (_mat[3][0]*v.x() + _mat[3][1]*v.y() + _mat[3][2]*v.z() + _mat[3][3]*v.w())) ;
}

inline TCVec4f TCMatrixf::preMult( const TCVec4f& v ) const
{
    return TCVec4f( (_mat[0][0]*v.x() + _mat[1][0]*v.y() + _mat[2][0]*v.z() + _mat[3][0]*v.w()),
                 (_mat[0][1]*v.x() + _mat[1][1]*v.y() + _mat[2][1]*v.z() + _mat[3][1]*v.w()),
                 (_mat[0][2]*v.x() + _mat[1][2]*v.y() + _mat[2][2]*v.z() + _mat[3][2]*v.w()),
                 (_mat[0][3]*v.x() + _mat[1][3]*v.y() + _mat[2][3]*v.z() + _mat[3][3]*v.w()));
}
inline TCVec4d TCMatrixf::preMult( const TCVec4d& v ) const
{
    return TCVec4d( (_mat[0][0]*v.x() + _mat[1][0]*v.y() + _mat[2][0]*v.z() + _mat[3][0]*v.w()),
                 (_mat[0][1]*v.x() + _mat[1][1]*v.y() + _mat[2][1]*v.z() + _mat[3][1]*v.w()),
                 (_mat[0][2]*v.x() + _mat[1][2]*v.y() + _mat[2][2]*v.z() + _mat[3][2]*v.w()),
                 (_mat[0][3]*v.x() + _mat[1][3]*v.y() + _mat[2][3]*v.z() + _mat[3][3]*v.w()));
}
inline TCVec3f TCMatrixf::transform3x3(const TCVec3f& v,const TCMatrixf& m)
{
    return TCVec3f( (m._mat[0][0]*v.x() + m._mat[1][0]*v.y() + m._mat[2][0]*v.z()),
                 (m._mat[0][1]*v.x() + m._mat[1][1]*v.y() + m._mat[2][1]*v.z()),
                 (m._mat[0][2]*v.x() + m._mat[1][2]*v.y() + m._mat[2][2]*v.z()));
}
inline TCVec3d TCMatrixf::transform3x3(const TCVec3d& v,const TCMatrixf& m)
{
    return TCVec3d( (m._mat[0][0]*v.x() + m._mat[1][0]*v.y() + m._mat[2][0]*v.z()),
                 (m._mat[0][1]*v.x() + m._mat[1][1]*v.y() + m._mat[2][1]*v.z()),
                 (m._mat[0][2]*v.x() + m._mat[1][2]*v.y() + m._mat[2][2]*v.z()));
}

inline TCVec3f TCMatrixf::transform3x3(const TCMatrixf& m,const TCVec3f& v)
{
    return TCVec3f( (m._mat[0][0]*v.x() + m._mat[0][1]*v.y() + m._mat[0][2]*v.z()),
                 (m._mat[1][0]*v.x() + m._mat[1][1]*v.y() + m._mat[1][2]*v.z()),
                 (m._mat[2][0]*v.x() + m._mat[2][1]*v.y() + m._mat[2][2]*v.z()) ) ;
}
inline TCVec3d TCMatrixf::transform3x3(const TCMatrixf& m,const TCVec3d& v)
{
    return TCVec3d( (m._mat[0][0]*v.x() + m._mat[0][1]*v.y() + m._mat[0][2]*v.z()),
                 (m._mat[1][0]*v.x() + m._mat[1][1]*v.y() + m._mat[1][2]*v.z()),
                 (m._mat[2][0]*v.x() + m._mat[2][1]*v.y() + m._mat[2][2]*v.z()) ) ;
}

inline void TCMatrixf::preMultTranslate( const TCVec3d& v )
{
    for (unsigned i = 0; i < 3; ++i)
    {
        double tmp = v[i];
        if (tmp == 0)
            continue;
        _mat[3][0] += tmp*_mat[i][0];
        _mat[3][1] += tmp*_mat[i][1];
        _mat[3][2] += tmp*_mat[i][2];
        _mat[3][3] += tmp*_mat[i][3];
    }
}

inline void TCMatrixf::preMultTranslate( const TCVec3f& v )
{
    for (unsigned i = 0; i < 3; ++i)
    {
        float tmp = v[i];
        if (tmp == 0)
            continue;
        _mat[3][0] += tmp*_mat[i][0];
        _mat[3][1] += tmp*_mat[i][1];
        _mat[3][2] += tmp*_mat[i][2];
        _mat[3][3] += tmp*_mat[i][3];
    }
}

inline void TCMatrixf::postMultTranslate( const TCVec3d& v )
{
    for (unsigned i = 0; i < 3; ++i)
    {
        double tmp = v[i];
        if (tmp == 0)
            continue;
        _mat[0][i] += tmp*_mat[0][3];
        _mat[1][i] += tmp*_mat[1][3];
        _mat[2][i] += tmp*_mat[2][3];
        _mat[3][i] += tmp*_mat[3][3];
    }
}

inline void TCMatrixf::postMultTranslate( const TCVec3f& v )
{
    for (unsigned i = 0; i < 3; ++i)
    {
        float tmp = v[i];
        if (tmp == 0)
            continue;
        _mat[0][i] += tmp*_mat[0][3];
        _mat[1][i] += tmp*_mat[1][3];
        _mat[2][i] += tmp*_mat[2][3];
        _mat[3][i] += tmp*_mat[3][3];
    }
}

inline void TCMatrixf::preMultScale( const TCVec3d& v )
{
    _mat[0][0] *= v[0]; _mat[0][1] *= v[0]; _mat[0][2] *= v[0]; _mat[0][3] *= v[0];
    _mat[1][0] *= v[1]; _mat[1][1] *= v[1]; _mat[1][2] *= v[1]; _mat[1][3] *= v[1];
    _mat[2][0] *= v[2]; _mat[2][1] *= v[2]; _mat[2][2] *= v[2]; _mat[2][3] *= v[2];
}

inline void TCMatrixf::preMultScale( const TCVec3f& v )
{
    _mat[0][0] *= v[0]; _mat[0][1] *= v[0]; _mat[0][2] *= v[0]; _mat[0][3] *= v[0];
    _mat[1][0] *= v[1]; _mat[1][1] *= v[1]; _mat[1][2] *= v[1]; _mat[1][3] *= v[1];
    _mat[2][0] *= v[2]; _mat[2][1] *= v[2]; _mat[2][2] *= v[2]; _mat[2][3] *= v[2];
}

inline void TCMatrixf::postMultScale( const TCVec3d& v )
{
    _mat[0][0] *= v[0]; _mat[1][0] *= v[0]; _mat[2][0] *= v[0]; _mat[3][0] *= v[0];
    _mat[0][1] *= v[1]; _mat[1][1] *= v[1]; _mat[2][1] *= v[1]; _mat[3][1] *= v[1];
    _mat[0][2] *= v[2]; _mat[1][2] *= v[2]; _mat[2][2] *= v[2]; _mat[3][2] *= v[2];
}

inline void TCMatrixf::postMultScale( const TCVec3f& v )
{
    _mat[0][0] *= v[0]; _mat[1][0] *= v[0]; _mat[2][0] *= v[0]; _mat[3][0] *= v[0];
    _mat[0][1] *= v[1]; _mat[1][1] *= v[1]; _mat[2][1] *= v[1]; _mat[3][1] *= v[1];
    _mat[0][2] *= v[2]; _mat[1][2] *= v[2]; _mat[2][2] *= v[2]; _mat[3][2] *= v[2];
}


inline void TCMatrixf::preMultRotate( const TCQuat& q )
{
    if (q.zeroRotation())
        return;
    TCMatrixf r;
    r.setRotate(q);
    preMult(r);
}

inline void TCMatrixf::postMultRotate( const TCQuat& q )
{
    if (q.zeroRotation())
        return;
    TCMatrixf r;
    r.setRotate(q);
    postMult(r);
}

inline TCVec3f operator* (const TCVec3f& v, const TCMatrixf& m )
{
    return m.preMult(v);
}
inline TCVec3d operator* (const TCVec3d& v, const TCMatrixf& m )
{
    return m.preMult(v);
}
inline TCVec4f operator* (const TCVec4f& v, const TCMatrixf& m )
{
    return m.preMult(v);
}
inline TCVec4d operator* (const TCVec4d& v, const TCMatrixf& m )
{
    return m.preMult(v);
}

inline TCVec3f TCMatrixf::operator* (const TCVec3f& v) const
{
    return postMult(v);
}
inline TCVec3d TCMatrixf::operator* (const TCVec3d& v) const
{
    return postMult(v);
}
inline TCVec4f TCMatrixf::operator* (const TCVec4f& v) const
{
    return postMult(v);
}
inline TCVec4d TCMatrixf::operator* (const TCVec4d& v) const
{
    return postMult(v);
}


#endif /* defined(__threecpp__TCMatrixf__) */
