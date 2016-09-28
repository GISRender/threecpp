//
//  TCPerspectiveCamera.cpp
//  threecpp
//
//  Created by CastingJ on 16/9/28.
//
//

#include "TCPerspectiveCamera.h"


TCPerspectiveCamera::TCPerspectiveCamera(float fovy,float aspect,float near,float far)
    :TCCamera(),
    _fovy(fovy),
    _aspect(aspect),
    _zNear(near),
    _zFar(far)
{
    
}

TCPerspectiveCamera::~TCPerspectiveCamera()
{
    
}

TCPerspectiveCamera::TCPerspectiveCamera()
{
    
}

void TCPerspectiveCamera::initProjectionMatrix()
{
    _projectionMatrix.makePerspective(_fovy, _aspect, _zNear, _zFar);
}