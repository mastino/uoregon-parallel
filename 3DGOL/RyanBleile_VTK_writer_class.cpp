// Ryan Bleile 
// My code is written using c++ and compiled with icpc or g++
//
// need to include :

#include <sstream>
#include <fstream>
#include <cstring>

/*
nx = number of cells in the x dimension
ny = number of cells in the y dimension
nz = number of cells in the z dimension
*/
void Mesh::PrintVtk( int iteration )
{

    stringstream smesh;
    smesh << "# vtk DataFile Version 3.0" << endl;
    smesh << "vtk output" << endl;
    smesh << "ASCII" << endl;
    smesh << "DATASET RECTILINEAR_GRID" << endl;
    smesh << "DIMENSIONS " << nx+1 << " " << ny+1 << " " << nz+1 << endl;
    smesh << "X_COORDINATES " << nx+1 << " float" << endl;
    for( int i = 0; i <= nx; i++ )    smesh << i << " ";
    smesh << endl;

    smesh << "Y_COORDINATES " << ny+1 << " float" << endl;
    for( int i = 0; i <= ny; i++ )    smesh << i << " ";
    smesh << endl;

    smesh << "Z_COORDINATES " << nz+1 << " float" << endl;
    for( int i = 0; i <= nz; i++ )    smesh << i << " ";
    smesh << endl;

    smesh << "CELL_DATA " << nx*ny*nz << endl;
    smesh << "SCALARS isAlive int" << endl;
    smesh << "LOOKUP_TABLE default" << endl;

//1D loop over cells. Can do this with a 3D loop and then comvert to a 1D index using : ci = x + y*nx + z*nx*ny
    for( unsigned int ci = 0; ci < nx*ny*nz; ci++ )
    {
        if( ci != 0 && ci%9 == 0)
            smesh << endl;

        smesh << cell[ci].isAlive() << " ";
    }

    stringstream filename;

    filename << "VTK/GameOfLife" << std::setw(6) << std::setfill('0') << iteration << ".vtk";
    string fname = filename.str();
    ofstream out( fname.c_str() );

    out << smesh.rdbuf();

    out.close();

} 
