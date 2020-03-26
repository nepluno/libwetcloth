//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
