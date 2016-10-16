
#ifndef  Mesh_h
#define  Mesh_h

#include <vector>

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
//#include "quaternion.h"

// ====== implementation only ?
#include "arrayAlgs.h"
#include "raytrace.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


class Disk3D{  // used for acceleration of raytracing
    public:
    Vec3d  center;
    double Rbound;
    Vec3d  normal;
    double Cplane;
};

class Polygon{
    public:
    //Disk3D * disk;
    std::vector<int> ipoints;
    std::vector<int> inormals;
    std::vector<int> iedges;
    //bool convex=false;
};

class Mesh{
    public:
    int    rendered_shape;
    Vec3d  center; // center of spherical bounding box
    double Rbound; // radius of spherical bounding box
    /*
    int npoints;
    int ntris;
    int npolys;
    int nedges;

    Vec3d  * points;
    Vec3i  * tris;
    int    * polyNs;
    double ** polys;
    Vec2i  * edges;
    */

    std::vector<Vec3d>      points;
    std::vector<Vec3d>      normals;
    std::vector<Vec3i>      triangles;
    std::vector<Polygon*>    polygons;
    //std::vector<double[4]>  planes;
    std::vector<Vec2i>      edges;

    Disk3D * disks;  // used for acceleration of raytracing; for each polygon there is one disk

    int fromFileOBJ( std::string fname ){
        //std::ofstream infile;
        std::ifstream infile ( fname );
        //myfile << "Writing this to a file.\n";
        std::string line,word;
        //std::stringstream ss;

        int npoints0   = points  .size();
        int npolygons0 = polygons.size();
        int i=0;
        while ( std::getline(infile,line) ){
            std::cout << line << '\n';
            //ss.str(line);
            std::stringstream ss(line);
            getline(ss, word, ' ');
            if ( word.compare("v")==0 ){
                Vec3d p;
                getline(ss, word, ' '); p.x = stof(word);
                getline(ss, word, ' '); p.y = stof(word);
                getline(ss, word, ' '); p.z = stof(word);
                points.push_back(p);
                i++;
            }else if ( word.compare("f")==0  ){
                Polygon * p = new Polygon();
                polygons.push_back(p);
                while ( getline(ss, word, ' ')) {
                    size_t of;
                    p->ipoints .push_back( stoi(word,&of)-1 );
                    //std::cout << word.substr(of) << '\n';
                    p->inormals.push_back( stoi(word.substr(of+2))-1 );
                };
                i++;
            }else if ( word.compare("vn")==0 ){
                Vec3d p;
                getline(ss, word, ' '); p.x = stof(word);
                getline(ss, word, ' '); p.y = stof(word);
                getline(ss, word, ' '); p.z = stof(word);
                normals.push_back(p);
                i++;
            };
        }
        infile.close();
        return i;
    }

    int pickVertex( const Vec3d &ray0, const Vec3d &hRay ){
        double r2min=1e+300;
        int imin=0;
        for(int i=0; i<points.size(); i++){
            double t;
            double r2 = rayPointDistance2( ray0, hRay, points[i], t );
            if(r2<r2min){ imin=i; r2min=r2; }
        }
        return imin;
    }

    double ray( const Vec3d &ray0, const Vec3d &hRay, Vec3d& normal ){
        //int    imin  = 0;

        Vec3d hX,hY;
        hRay.getSomeOrtho( hX, hY );

        double t_min = 1e+300;
        //Vec3d hitpos_min,normal_min;
        for(int i=0; i<triangles.size(); i++ ){
            Vec3i itri = triangles[i];
            Vec3d A = points[itri.x];
            Vec3d B = points[itri.y];
            Vec3d C = points[itri.z];
            Vec3d normal_;
            bool inside_;
            //double t = rayTriangle( ray0, hRay, A, B, C, inside_, hitpos_, normal_ );
            double t = rayTriangle2( ray0, hRay, hX, hY, A, B, C, normal_ );
            //printf( "t=%f\n", t );
            inside_ = (t<0.9e+300 )&&(t>0);
            if( inside_ && ( t<t_min ) ){
                t_min = t;
                normal = normal_;
            }
        };

    };

};

#endif



