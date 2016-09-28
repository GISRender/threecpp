//
//  TCView.h
//  threecpp
//
//  Created by CastingJ on 16/9/28.
//
//

#ifndef __threecpp__TCView__
#define __threecpp__TCView__

class TCCamera;

class TCView
{
public:
    TCView();
    virtual ~TCView();
    
    void update();
    
protected:
    TCCamera *_camera;
};

#endif /* defined(__threecpp__TCView__) */
