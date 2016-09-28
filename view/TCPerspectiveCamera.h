//
//  TCPerspectiveCamera.h
//  threecpp
//
//  Created by CastingJ on 16/9/28.
//
//

#ifndef __threecpp__TCPerspectiveCamera__
#define __threecpp__TCPerspectiveCamera__

#include "TCCamera.h"

class TCPerspectiveCamera : public TCCamera
{
public:
    TCPerspectiveCamera(float fovy,float aspect,float znear,float zfar);
    
    virtual ~TCPerspectiveCamera();
    
protected:
    TCPerspectiveCamera();
    
    void initProjectionMatrix();
    
private:
    float _fovy;
    float _aspect;
    float _zNear;
    float _zFar;
};

#endif /* defined(__threecpp__TCPerspectiveCamera__) */
