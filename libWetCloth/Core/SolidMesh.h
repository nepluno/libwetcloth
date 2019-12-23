//
// This file is part of the libWetCloth open source project
//
// Copyright 2018 Yun (Raymond) Fei, Christopher Batty, Eitan Grinspun, and Changxi Zheng
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef SOLID_MESH_H__
#define SOLID_MESH_H__

#include "MathDefs.h"
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include "StringUtilities.h"

class SolidMesh
{
protected:
	std::vector<Vector3s> vertices;
	std::vector<Vector3i> indices;
public:
	const std::vector<Vector3s>& getVertices() const { return vertices; }
	const std::vector<Vector3i>& getIndices() const { return indices; }
    
    SolidMesh() {};
    
    SolidMesh(const std::string& szfn, const scalar scaling = 1.0) {
        ReadOBJ(szfn);
        for(Vector3s& v : vertices)
        {
            v *= scaling;
        }
    }
    
    void translate(const Vector3s& t)
    {
        for(Vector3s& v : vertices)
        {
            v += t;
        }
    }
    
    Vector3s getCenter() const
    {
        Vector3s ret = Vector3s::Zero();
        for(const Vector3s& v : vertices)
        {
            ret += v;
        }
        
        ret /= (scalar) vertices.size();
        return ret;
    }
    
    void boundingBox(Vector3s& bbx_min, Vector3s& bbx_max) const
    {
        if(!vertices.size()) {
            bbx_min = bbx_max = Vector3s::Zero();
            return;
        }
        
        bbx_min = Vector3s::Constant(std::numeric_limits<scalar>::max());
        bbx_max = Vector3s::Constant(std::numeric_limits<scalar>::min());
        
        for(const Vector3s& v : vertices)
        {
            bbx_min(0) = std::min(bbx_min(0), v(0));
            bbx_min(1) = std::min(bbx_min(1), v(1));
            bbx_min(2) = std::min(bbx_min(2), v(2));
            
            bbx_max(0) = std::max(bbx_max(0), v(0));
            bbx_max(1) = std::max(bbx_max(1), v(1));
            bbx_max(2) = std::max(bbx_max(2), v(2));
        }
    }
    
    void boundingBox(Vector3s& bbx_min, Vector3s& bbx_max, const Eigen::Quaternion<scalar>& rot) const
    {
        if(!vertices.size()) {
            bbx_min = bbx_max = Vector3s::Zero();
            return;
        }
        
        bbx_min = Vector3s::Constant(std::numeric_limits<scalar>::max());
        bbx_max = Vector3s::Constant(std::numeric_limits<scalar>::min());
        
        for(const Vector3s& v : vertices)
        {
            Vector3s rv = rot * v;
            
            bbx_min(0) = std::min(bbx_min(0), rv(0));
            bbx_min(1) = std::min(bbx_min(1), rv(1));
            bbx_min(2) = std::min(bbx_min(2), rv(2));
            
            bbx_max(0) = std::max(bbx_max(0), rv(0));
            bbx_max(1) = std::max(bbx_max(1), rv(1));
            bbx_max(2) = std::max(bbx_max(2), rv(2));
        }
    }
    
	void WriteOBJ(const std::string& szfn) const
	{
		std::ofstream fout(szfn);
		WriteOBJ(fout);
		fout.close();
	}
	void WriteOBJ(std::ostream& str) const
	{
		// vertices
		for(const Vector3s& v : vertices)
		{
			str << "v " << v.transpose() << std::endl;
		}
		
		for(const Vector3i& f : indices)
		{
			str << "f " << f(0)+1 << "//" << f(0)+1 << " " << f(1)+1 << "//" << f(1)+1 << " " << f(2)+1 << "//" << f(2)+1 << std::endl;
		}
	}
    
    void ReadOBJ(const std::string& szfn)
    {
        std::ifstream fin(szfn);
        ReadOBJ(fin);
        fin.close();
    }
    
    void ReadOBJ(std::istream& str)
    {
        std::string obj_command;
        
        int line_num = 0;
        while(!str.eof())
        {
            getline( str, obj_command );
            
            std::istringstream commandstream(obj_command);
            
            // First element of a command is the command's name
            std::string command;
            commandstream >> command;
            
            if( command == "v" ) {
                scalar x, y, z;
                if( !(commandstream >> x >> y >> z) )
                {
                    std::cout << "OBJPARSER: Invalid vertex command on line " << line_num << ", vertex must have three real coordinates." << std::endl;
                    continue;
                }
                
                vertices.push_back(Vector3s(x, y, z));
            } else if( command == "f" ) {
                int xidx,yidx,zidx;
                
                std::vector<std::string> vertstrngs;
                std::string vertcmmnd;
                while( commandstream >> vertcmmnd ) vertstrngs.push_back(vertcmmnd);
                
                if( vertstrngs.size() != 3 )
                {
                    std::cout << vertstrngs.size() << std::endl;
                    std::cout << "OBJPARSER: Invalid face command on line " << line_num << ", face must have three vertices." << std::endl;
                    continue;
                }
                
                std::string xstr = vertstrngs[0];
                std::string ystr = vertstrngs[1];
                std::string zstr = vertstrngs[2];
                
                // Attempt to extract vertex indices for x
                std::vector<std::string> vertex_spec;
                stringutils::tokenize( xstr, vertex_spec, "/" );
                if( vertex_spec.size() < 1 || vertex_spec.size() > 3 )
                {
                    std::cout << "OBJPARSER: Invalid face command on line " << line_num << ", x vertex specification must be of form v[/vt][/vn]." << std::endl;
                    continue;
                }
                std::istringstream xstream(vertex_spec[0]);
                if( !(xstream>>xidx) )
                {
                    std::cout << "OBJPARSER: Invalid face command on line " << line_num << ", x vertex index must be an integer." << std::endl;
                    continue;
                }
                if( xidx > (int) vertices.size() )
                {
                    std::cout << "OBJPARSER: Invalid face command on line " << line_num << ", x vertex index greater than number of vertices." << std::endl;
                    continue;
                }
                if( xidx < 0 ) if( xidx + ((int) vertices.size() ) < 0 )
                {
                    std::cout << "OBJPARSER: Invalid face command on line " << line_num << ", relative x vertex index out of bounds." << std::endl;
                    continue;
                }
                // Attempt to extract vertex indices for y
                vertex_spec.clear();
                stringutils::tokenize( ystr, vertex_spec, "/" );
                if( vertex_spec.size() < 1 || vertex_spec.size() > 3 )
                {
                    std::cout << "OBJPARSER: Invalid face command on line " << line_num << ", y vertex specification must be of form v[/vt][/vn]." << std::endl;
                    continue;
                }
                std::istringstream ystream(vertex_spec[0]);
                if( !(ystream>>yidx) )
                {
                    std::cout << "OBJPARSER: Invalid face command on line " << line_num << ", y vertex index must be an integer." << std::endl;
                    continue;
                }
                if( yidx > (int) vertices.size() )
                {
                    std::cout << "OBJPARSER: Invalid face command on line " << line_num << ", y vertex index greater than number of vertices." << std::endl;
                    continue;
                }
                if( yidx < 0 ) if(  yidx + ( (int) vertices.size() )< 0 )
                {
                    std::cout << "OBJPARSER: Invalid face command on line " << line_num << ", relative y vertex index out of bounds." << std::endl;
                    continue;
                }
                // Attempt to extract vertex indices for z
                vertex_spec.clear();
                stringutils::tokenize( zstr, vertex_spec, "/" );
                if( vertex_spec.size() < 1 || vertex_spec.size() > 3 )
                {
                    std::cout << "OBJPARSER: Invalid face command on line " << line_num << ", z vertex specification must be of form v[/vt][/vn]." << std::endl;
                    continue;
                }
                std::istringstream zstream(vertex_spec[0]);
                if( !(zstream>>zidx) )
                {
                    std::cout << "OBJPARSER: Invalid face command on line " << line_num << ", z vertex index must be an integer." << std::endl;
                    continue;
                }
                if( zidx > (int) vertices.size() )
                {
                    std::cout << "OBJPARSER: Invalid face command on line " << line_num << ", z vertex index greater than number of vertices." << std::endl;
                    continue;
                }
                if( zidx < 0 ) if( zidx + ( (int) vertices.size() ) < 0 )
                {
                    std::cout << "OBJPARSER: Invalid face command on line " << line_num << ", relative z vertex index out of bounds." << std::endl;
                    continue;
                }
                
                // Generate a triangular face
                indices.push_back(Vector3i(xidx-1,yidx-1,zidx-1));
            }
            
            ++line_num;
        }
    }
};

#endif
