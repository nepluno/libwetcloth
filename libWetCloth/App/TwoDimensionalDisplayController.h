//
// This file is part of the libWetCloth open source project
//
// The code is licensed under the same terms as a Clear BSD License but further
// restricted to academic and non-commercial use (commercial licenses may be
// obtained by contacting the faculty of the Columbia Computer Graphics Group
// or Columbia Technology Ventures).
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the disclaimer
// below) provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
// list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// * Neither the name of the copyright holder nor the names of its contributors may be used
// to endorse or promote products derived from this software without specific
// prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY THIS
// LICENSE. THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
// GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.

#ifndef __TWO_DIMENSIONAL_DISPLAY_CONTROLLER_H__
#define __TWO_DIMENSIONAL_DISPLAY_CONTROLLER_H__

#include <cmath>
#include <iostream>
#include <cassert>

#include "RenderingUtilities.h"
#include "Camera.h"

class TwoDimensionalDisplayController
{
public:
	TwoDimensionalDisplayController( int width, int height );
	
	void reshape( int w, int h );
	
	void keyboard( unsigned char key, int x, int y );
	
	void special( int key, int x, int y );
	
	void mouse( int button, int state, int x, int y );
	
	void motion( int x, int y );
	
	int getWindowWidth() const;
	int getWindowHeight() const;
	
	int setWindowWidth(int w);
	int setWindowHeight(int h);
	
	int getWorldWidth() const;
	int getWorldHeight() const;
	
	double getCenterX() const;
	double getCenterY() const;
	double getCenterZ() const;
	void setCenterX( double x );
	void setCenterY( double y );
	void setCenterZ( double z );
	
	void setScaleFactor( double scale );
	double getScaleFactor() const;
	void applyProjection() const;
	double getRatio() const;
	
	void initCamera(const renderingutils::Viewport& view);
	
	const Camera& getCamera() const;
	
	Camera& getCamera();
    
    int currentCameraIndex() const;
private:
	// Width of the window in pixels
	int m_window_width;
	// Height of the window in pixels
	int m_window_height;
	// Factor to 'zoom' in or out by
	double m_scale_factor;
	// Center of the display, x coord
	double m_center_x;
	// Center of the display, y coord
	double m_center_y;
	// Center of the display, z coord
	double m_center_z;
	// True if the user is dragging the display left
	bool m_left_drag;
	// True if the user is dragging the display right
	bool m_right_drag;
	// Last position of the cursor in a drag, x coord
	double m_last_x;
	// Last position of the cursor in a drag, y coord
	double m_last_y;
	
	bool m_right_part_click;
	
	unsigned m_modifiers;
	
	std::vector<Camera> m_cameras;
	int m_current_camera;
};

#endif
