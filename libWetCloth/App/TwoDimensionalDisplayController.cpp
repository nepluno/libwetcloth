//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


#include "TwoDimensionalDisplayController.h"
#include "ParticleSimulation.h"

TwoDimensionalDisplayController::TwoDimensionalDisplayController( int width, int height )
	: m_window_width(width)
	, m_window_height(height)
	, m_scale_factor(1.0)
	, m_center_x(0.0)
	, m_center_y(0.0)
	, m_center_z(0.0)
	, m_left_drag(false)
	, m_right_drag(false)
	, m_last_x(0)
	, m_last_y(0)
	, m_current_camera(-1)
	, m_modifiers(0)
	, m_right_part_click(false)
{
	m_cameras.push_back(Camera());
	m_current_camera = 0;
}

int TwoDimensionalDisplayController::setWindowWidth(int w)
{
	m_window_width = w;

	return w;
}

int TwoDimensionalDisplayController::setWindowHeight(int h)
{
	m_window_height = h;

	return h;
}

double TwoDimensionalDisplayController::getRatio() const
{
	return (double) m_window_width / (double) m_window_height;
}

void TwoDimensionalDisplayController::applyProjection() const
{
#ifdef RENDER_ENABLED
	// Reset the coordinate system before modifying
	if (g_rendering_enabled) {
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();

		double fov, zNear, zFar;
		m_cameras[m_current_camera].getPerspective( fov, zNear, zFar );
		gluPerspective( fov, getRatio(),  zNear, zFar );

		glMatrixMode( GL_MODELVIEW );
		glLoadIdentity();
		Eigen::Vector3d eye, center, up;
		m_cameras[m_current_camera].getLookAt( eye, center, up );
		gluLookAt( eye.x(), eye.y(), eye.z(), center.x(), center.y(), center.z(), up.x(), up.y(), up.z() );

		glutPostRedisplay();
	}
#endif
}

void TwoDimensionalDisplayController::reshape( int w, int h )
{

	assert( renderingutils::checkGLErrors() );
	// Record the new width and height
	m_window_width = w;
	m_window_height = h;
#ifdef RENDER_ENABLED
	if (g_rendering_enabled)
		glViewport(0, 0, getWindowWidth(), getWindowHeight());

	applyProjection();

	assert( renderingutils::checkGLErrors() );
#endif
}

int TwoDimensionalDisplayController::getWorldWidth() const
{
	double ratio = getRatio();
	return 2 * m_scale_factor / ratio;
}

int TwoDimensionalDisplayController::getWorldHeight() const
{
	return 2 * m_scale_factor;
}

int TwoDimensionalDisplayController::currentCameraIndex() const
{
	return m_current_camera;
}

void TwoDimensionalDisplayController::keyboard( unsigned char key, int x, int y )
{
#ifdef RENDER_ENABLED
	if ( key == '-' || key == '_' )
	{
		m_scale_factor += 0.1;
		reshape(m_window_width, m_window_height);
	}
	else if ( key == '=' || key == '+' )
	{
		m_scale_factor = std::max(0.1, m_scale_factor - 0.1);
		reshape(m_window_width, m_window_height);
	} else if (key == 'K' || key == 'k' )
	{
		m_cameras.push_back(m_cameras[m_current_camera]);
		m_current_camera = m_cameras.size() - 1;
		glutPostRedisplay();
	}
	else if ( key == 'L' || key == 'l' )
	{
		m_current_camera = (m_current_camera + 1) % m_cameras.size();
		applyProjection();
	}
	else if ( key == 'I' || key == 'i' )
	{
		std::cout << m_cameras[m_current_camera] << std::endl;
	}
#endif
}

void TwoDimensionalDisplayController::special( int key, int x, int y )
{
#ifdef RENDER_ENABLED
	if ( GLUT_KEY_UP == key )
	{
		m_center_y += 0.1;
		reshape(m_window_width, m_window_height);
	}
	else if ( GLUT_KEY_DOWN == key )
	{
		m_center_y -= 0.1;
		reshape(m_window_width, m_window_height);
	}
	else if ( GLUT_KEY_LEFT == key )
	{
		m_center_x -= 0.1;
		reshape(m_window_width, m_window_height);
	}
	else if ( GLUT_KEY_RIGHT == key )
	{
		m_center_x += 0.1;
		reshape(m_window_width, m_window_height);
	}
#endif
}

void TwoDimensionalDisplayController::mouse( int button, int state, int x, int y )
{
#ifdef RENDER_ENABLED
	if ( !m_right_drag && button == GLUT_LEFT_BUTTON && state == GLUT_DOWN )
	{
		m_right_part_click = (x > m_window_width / 2);

		m_left_drag = true;
		int iPart = 1;
		int mx = x % (m_window_width / iPart);
		int my = y;

		double r = (m_window_width / iPart) < m_window_height ? (m_window_width / iPart) : m_window_height;

		double nx = (2.0 * mx - m_window_width / iPart) / r;
		double ny = (m_window_height - 2.0 * my) / r;

		m_last_x = nx;
		m_last_y = ny;

		m_modifiers = glutGetModifiers();
	}
	if ( button == GLUT_LEFT_BUTTON && state == GLUT_UP )
	{
		m_right_part_click = false;
		m_left_drag = false;
		m_modifiers = 0;
	}

	if ( !m_left_drag && button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN )
	{
		m_right_part_click = (x > m_window_width / 2);

		int iPart = 1;
		int mx = x % (m_window_width / iPart);
		int my = y;

		double r = (m_window_width / iPart) < m_window_height ? (m_window_width / iPart) : m_window_height;

		double nx = (2.0 * mx - m_window_width / iPart) / r;
		double ny = (m_window_height - 2.0 * my) / r;

		m_last_x = nx;
		m_last_y = ny;
		m_right_drag = true;
		m_modifiers = glutGetModifiers();

	}
	if ( button == GLUT_RIGHT_BUTTON && state == GLUT_UP )
	{
		m_right_part_click = false;
		m_right_drag = false;
		m_modifiers = 0;
	}
#endif
}

void TwoDimensionalDisplayController::motion( int x, int y )
{
#ifdef RENDER_ENABLED
	if ( m_left_drag )
	{

		double ox = m_last_x;
		double oy = m_last_y;

		int iPart = 1;

		int mx = x % (m_window_width / iPart);
		int my = y;

		double r = (m_window_width / iPart) < m_window_height ? (m_window_width / iPart) : m_window_height;

		double nx = (2.0 * mx - m_window_width / iPart) / r;
		double ny = (m_window_height - 2.0 * my) / r;

		m_last_x = nx;
		m_last_y = ny;
		if (m_modifiers & GLUT_ACTIVE_SHIFT) {
			m_cameras[m_current_camera].pan(ox, oy, nx, ny);
		} else {
			m_cameras[m_current_camera].rotate(ox, oy, nx, ny);
		}
		applyProjection();


	}
	if ( m_right_drag )
	{
		double ox = m_last_x;
		double oy = m_last_y;

		int iPart = 1;

		int mx = x % (m_window_width / iPart);
		int my = y;

		double r = (m_window_width / iPart) < m_window_height ? (m_window_width / iPart) : m_window_height;

		double nx = (2.0 * mx - m_window_width / iPart) / r;
		double ny = (m_window_height - 2.0 * my) / r;

		m_last_x = nx;
		m_last_y = ny;

		m_cameras[m_current_camera].zoom(ox, oy, nx, ny);

		applyProjection();
	}
#endif
}

int TwoDimensionalDisplayController::getWindowWidth() const
{
	return m_window_width;
}

int TwoDimensionalDisplayController::getWindowHeight() const
{
	return m_window_height;
}

double TwoDimensionalDisplayController::getCenterX() const
{
	return m_center_x;
}

double TwoDimensionalDisplayController::getCenterY() const
{
	return m_center_y;
}

double TwoDimensionalDisplayController::getCenterZ() const
{
	return m_center_z;
}

void TwoDimensionalDisplayController::setCenterX( double x )
{
	m_center_x = x;
}

void TwoDimensionalDisplayController::setCenterY( double y )
{
	m_center_y = y;
}

void TwoDimensionalDisplayController::setCenterZ( double z )
{
	m_center_z = z;
}

void TwoDimensionalDisplayController::setScaleFactor( double scale )
{
	m_scale_factor = scale;
}

double TwoDimensionalDisplayController::getScaleFactor() const
{
	return m_scale_factor;
}

void TwoDimensionalDisplayController::initCamera(const renderingutils::Viewport& view)
{
	Vector3s bmin(view.cx - view.rx, view.cy - view.ry, view.cz - view.rz);
	Vector3s bmax(view.cx + view.rx, view.cy + view.ry, view.cz + view.rz);
	m_cameras[m_current_camera].init(bmin, bmax);
}

const Camera& TwoDimensionalDisplayController::getCamera() const
{
	return m_cameras[m_current_camera];
}

Camera& TwoDimensionalDisplayController::getCamera()
{
	return m_cameras[m_current_camera];
}
