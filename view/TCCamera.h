//
//  TCCamera.h
//  threecpp
//
//  Created by CastingJ on 16/9/28.
//
//

#ifndef __threecpp__TCCamera__
#define __threecpp__TCCamera__

#include "TCMatrixd.h"

class TCCamera
{
public:
    TCCamera();
    virtual ~TCCamera();
 
    TCCamera(const TCCamera& camera);
    
protected:
    TCMatrixd _projectionMatrix;
    
};

#endif /* defined(__threecpp__TCCamera__) */
