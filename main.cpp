#include "TriangulateWires.hpp"
#include <unistd.h>
#include <omp.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////

void Usage( const string &executable)
{
    cout << "Usage : " << executable << " Arguments " << endl;
    cout << "-i  : An input file containing mesh edges  (.off, ,mtl) " << endl;
    cout << "    : with off file, specify beam radius" << endl;
    cout << "-n  : Number of sides on a beam (>=3) " << endl;
    cout << "-o  : An output off file " << endl;
    cout << "-r  : Radius of the wire" << endl;
    cout << "-R  : Read beam radius from the file" << endl;
    cout << "     (This option overrides previously set beam radius)" << endl;
    cout << "-v  : Reduce node volume (0 No: 1 Yes" << endl;
    cout << "-x  : Allow intersections at balls and beams " << endl;
    cout << "-t  : Number of openmp threads  " << endl;
    cout << "-h  : Help" << endl;
    cout << "Examples: " << endl;
    cout << "  " << executable << " -i model.off -n 6 -o model.off -r 0.01 "<< endl;
}

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{
    string infile, outfile = "lattice.off";
    string radiusFile;
    int    polySides = 6;
    int    minArgs   = 0;
    bool   reduceVolume = 1;
    bool   manifold_surface = 1;
    int    numThreads = 1;
    double beamRadius = 0;

    int c;
    while ((c = getopt (argc, argv, "n:i:o:r:R:t:xv:h")) != -1) {
        switch (c)
        {
        case 'n':
            polySides = atoi(optarg);
            break;
        case 'i':
            infile = optarg;
            break;
        case 'o':
            outfile = optarg;
            break;
        case 'r':
            beamRadius  = stod(optarg);
            break;
        case 'R':
            radiusFile  = optarg;
            break;
        case 'v':
            reduceVolume  = atoi(optarg);
            break;
        case 'x':
            manifold_surface = 0;
            break;
        case 't':
            numThreads = atoi(optarg);
            break;
        case 'h':
            Usage( argv[0] );
            return 1;
        }
    }

    if( infile.empty() ) {
        Usage( argv[0]);
        return 1;
    }

    omp_set_num_threads(numThreads);

    TriangulateWires triwires;
    triwires.readGraph( infile );
    if( !radiusFile.empty() ) triwires.setBeamRadius( radiusFile);
    triwires.setPolySides(polySides);
    triwires.setBeamRadius(beamRadius);
    triwires.setUnion(manifold_surface);
    triwires.setReduceNodeVolume(reduceVolume);

    triwires.triangulate();
    triwires.saveAs( outfile );
    cout << "Done " << endl;

    return 0;
}
