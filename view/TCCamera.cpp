//
//  TCCamera.cpp
//  threecpp
//
//  Created by CastingJ on 16/9/28.
//
//

#include "TCCamera.h"


TCCamera::TCCamera()
{
    
}

TCCamera::~TCCamera()
{
    
}

TCCamera::TCCamera(const TCCamera& camera)
{
    _projectionMatrix = camera._projectionMatrix;
}