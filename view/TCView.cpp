//
//  TCView.cpp
//  threecpp
//
//  Created by CastingJ on 16/9/28.
//
//

#include "TCView.h"

#include <stdlib.h>

#include "TCCamera.h"

TCView::TCView()
    :_camera(NULL)
{
    
}

TCView::~TCView()
{
    if (_camera)
        delete _camera;
    
}

void update();
