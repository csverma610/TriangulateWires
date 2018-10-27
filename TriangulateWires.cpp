#include "TriangulateWires.hpp"
#include <iomanip>
#include <omp.h>

#ifdef BOOST_TIMER
#include <boost/timer/timer.hpp>
#endif

using namespace std;

using Vec3D = quickhull::Vector3<double>;

////////////////////////////////////////////////////////////////////////////////////

void TriangulateWires::alignAlong(vector<Array3D> &points, const Array3D &srcVec, const Array3D &dstVec)
{
    Array3D  perpAxis;
    cross_product( srcVec, dstVec, perpAxis);

    double dl = magnitude(perpAxis);
    if( fabs(dl) < 1.0E-15) return;
    perpAxis[0] /= dl;
    perpAxis[1] /= dl;
    perpAxis[2] /= dl;

    double t = getVecAngle(srcVec, dstVec);
    if( fabs(t) < 1.0E-15) return;

    double qcos = cos(0.5*t);
    double qsin = sin(0.5*t);

    Quaternion<double> q(qcos, qsin*perpAxis[0], qsin*perpAxis[1], qsin*perpAxis[2]);
    Quaternion<double> q1 = q.conjugate();
    Quaternion<double> result;

    for( int i = 0; i < points.size(); i++) {
        Array3D p3d = points[i];
        Quaternion<double> v(0.0, p3d[0], p3d[1], p3d[2]);
        result = q*v*q1;
        p3d[0] = result.get(1);
        p3d[1] = result.get(2);
        p3d[2] = result.get(3);
        points[i] = p3d;
    }
}

////////////////////////////////////////////////////////////////////////////////////

void TriangulateWires::readMTL( const string &filename)
{
    ifstream infile(filename.c_str(), ios::in);
    if( infile.fail() ) return;

    size_t numNodes;
    infile >> numNodes;
    junctions.resize(numNodes);
    string str;
    double x, y, z, r;
    for(size_t i = 0; i < numNodes; i++) {
        infile >> x >> y >> z >> r;
        getline(infile, str);
        junctions[i].xyz = {x,y,z};
        junctions[i].radius = r;
    }

    size_t numEdges;
    infile >> numEdges;

    beams.resize(numEdges);
    size_t n0, n1;
    for(size_t i = 0; i < numEdges; i++) {
        infile >> n0 >> n1;
        beams[i].nodes = {n0,n1};
    }
    double r0, r1;

    for( size_t i = 0; i < numEdges; i++) {
        int v0 = beams[i].nodes[0];
        int v1 = beams[i].nodes[1];
        if( min(v0,v1) == v0) {
            r0 = junctions[v0].radius;
            r1 = junctions[v1].radius;
        }

        if( min(v0,v1) == v1) {
            r0 = junctions[v1].radius;
            r1 = junctions[v0].radius;
        }
        beams[i].profile[0].radius = r0;
        beams[i].profile[1].radius = r1;
    }
    nodeRadiusAssigned = 1;
    edgeRadiusAssigned = 1;
}
////////////////////////////////////////////////////////////////////////////////////

void TriangulateWires::readOFF(const string &filename)
{
    junctions.clear();
    beams.clear();

    ifstream ifile(filename.c_str(), ios::in);
    if( ifile.fail() ) {
        cout << "Error: Can't open input file " << endl;
        return;
    }

    string str;
    ifile >> str;
    if( str != "OFF") {
        cout << "Error: input file not in the off format" << endl;
        return;
    }

    int numPoints, numFaces, numEdges;
    ifile >> numPoints >> numFaces >> numEdges;

    if( numPoints < 1) {
        cout << "Warning: Input file has no points " << endl;
        return;
    }

    if( numPoints ) {
        double x, y, z;
        junctions.resize(numPoints);
        for( size_t i = 0; i < numPoints; i++) {
            ifile >> x >> y >> z;
            junctions[i].xyz = {x,y,z};
        }
    }

    size_t n0, n1;
    if( numFaces) {
        set<pair<size_t,size_t>> eSet;
        int nn;
        vector<int> fnodes;
        for( size_t i = 0; i < numFaces; i++) {
            ifile >> nn;
            fnodes.resize(nn);
            for( int j = 0; j < nn; j++) ifile >> fnodes[j];
            for( int j = 0; j < nn; j++) {
                const auto v0 = fnodes[(j+1)%nn];
                const auto v1 = fnodes[(j+2)%nn];
                const auto vmin = min(v0,v1);
                const auto vmax = max(v0,v1);
                eSet.insert( make_pair(vmin,vmax));
            }
        }
        beams.resize(eSet.size());
        int index = 0;
        for(auto edge: eSet) {
            n0 = edge.first;
            n1 = edge.second;
            assert( n0 != n1);
            beams[index].nodes = {n0,n1};
            index++;
        }
    }

    if( numEdges ) {
        beams.resize(numEdges);
        for(size_t i = 0; i < numEdges; i++) {
            ifile >> n0 >> n1;
            assert(n0 != n1);
            beams[i].nodes = {n0,n1};
        }
    }

    edgeRadiusAssigned = 0;
    nodeRadiusAssigned = 0;
}

////////////////////////////////////////////////////////////////////////////////////

void TriangulateWires::readGraph( const string &filename)
{
    if( filename.rfind(".off") != std::string::npos) readOFF( filename);
    if( filename.rfind(".mtl") != std::string::npos) readMTL( filename);

    size_t numEdges = beams.size();

    minEdgeLength =  numeric_limits<double>::max();
    maxEdgeLength = -minEdgeLength;

    for( size_t i  = 0; i < numEdges; i++) {
        auto v0 =  beams[i].nodes[0];
        auto v1 =  beams[i].nodes[1];
        const Array3D &p0 = junctions[v0].xyz;
        const Array3D &p1 = junctions[v1].xyz;

        double dx = p1[0] - p0[0];
        double dy = p1[1] - p0[1];
        double dz = p1[2] - p0[2];
        double dl = sqrt(dx*dx + dy*dy + dz*dz);
        beams[i].length = dl;
        minEdgeLength = min(minEdgeLength, dl);
        maxEdgeLength = max(maxEdgeLength, dl);
    }

    assert( minEdgeLength > 1.0E-15);

    // Create node-edge relation ...
    for( size_t i  = 0; i < numEdges; i++) {
        auto v0 =  beams[i].nodes[0];
        auto v1 =  beams[i].nodes[1];
        junctions[v0].adjEdges.insert(i);
        junctions[v1].adjEdges.insert(i);
    }
}

////////////////////////////////////////////////////////////////////////////////
void TriangulateWires:: checkAndCorrectBeamRadius(void)
{
    // once the beam radii have been assigned, ensure that they are no more than 0.25 edge length
    double beamLength,maxRadius;
    int  numEdges = beams.size();
    int v0,v1;

    for( int edgeID = 0; edgeID < numEdges; edgeID++) {
        v0 = beams[edgeID].nodes[0];
        v1 =  beams[edgeID].nodes[1];

        const Array3D &p0 = junctions[v0].xyz;
        const Array3D &p1 = junctions[v1].xyz;
        beamLength = length(p0, p1);
        maxRadius = 0.25*beamLength;
        if (beams[edgeID].profile[0].radius > maxRadius){
            beams[edgeID].profile[0].radius = maxRadius;
        }
        if (beams[edgeID].profile[1].radius > maxRadius){
            beams[edgeID].profile[1].radius = maxRadius;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void TriangulateWires:: setBeamRadius(const string &filename)
{
    ifstream ifile(filename.c_str(), ios::in);
    if( ifile.fail() ) {
        cout << "Error: Can't open input file " << endl;
        return;
    }
    int id;
    double r;
    for( int i = 0; i < beams.size(); i++) {
        ifile >> id >> r;
        setBeamRadius(id,r);
    }
    edgeRadiusAssigned = 1;
}
////////////////////////////////////////////////////////////////////////////////

void TriangulateWires:: setBeamRadius(double r)
{
    if( r == 0.0) return;

    uniformBeamRadius = r;
    int  numEdges = beams.size();
    #pragma omp parallel for
    for( int i = 0; i < numEdges; i++) {
        beams[i].profile[0].radius = r;
        beams[i].profile[1].radius = r;
    }
    edgeRadiusAssigned = 1;
    nodeRadiusAssigned = 0;
}
////////////////////////////////////////////////////////////////////////////////
void TriangulateWires:: setBeamRadius( const vector<double> &r)
{
    int numEdges = beams.size();
    #pragma omp parallel for
    for( int i = 0; i < numEdges; i++) {
        beams[i].profile[0].radius = r[i];
        beams[i].profile[1].radius = r[i];
    }
    edgeRadiusAssigned = 1;
    nodeRadiusAssigned = 0;
}

////////////////////////////////////////////////////////////////////////////////

bool TriangulateWires :: isNormalProfile( const vector<Array3D> &points, const Array3D &axis)
{
    cout << "Check Norm" << endl;
    double xc = 0.0;
    double yc = 0.0;
    double zc = 0.0;

    int nPoints = points.size();

    for( int i = 0; i < nPoints; i++) {
        xc += points[i][0];
        yc += points[i][1];
        zc += points[i][2];
    }
    xc /= (double)nPoints;
    yc /= (double)nPoints;
    zc /= (double)nPoints;

    Array3D vec;
    vec[0] = points[0][0] - xc;
    vec[1] = points[0][1] - yc;
    vec[2] = points[0][2] - zc;

    double dproduct = fabs(vec[0]*axis[0] + vec[1]*axis[1] + vec[2]*axis[2]);

    if( dproduct < 1.0E-05) return 1;

    cout << "Dot Product " << dproduct << endl;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

void TriangulateWires :: setProfileCenter(size_t edgeID, size_t nodeID)
{
    int v0 = min(beams[edgeID].nodes[0], beams[edgeID].nodes[1] );
    int v1 = max(beams[edgeID].nodes[0], beams[edgeID].nodes[1] );

    int id = 0;
    if( v1 == nodeID ) id = 1;

    const Array3D &p0 = junctions[v0].xyz;
    const Array3D &p1 = junctions[v1].xyz;

    double dx = p1[0] - p0[0];
    double dy = p1[1] - p0[1];
    double dz = p1[2] - p0[2];
    double dl = sqrt(dx*dx + dy*dy + dz*dz);

    if( dl < 1.0E-10) {
        cout << "P0 " << p0[0] << " " << p0[1] << " " << p0[2] << endl;
        cout << "P1 " << p1[0] << " " << p1[1] << " " << p1[2] << endl;
        cout << "Warning: Beam not triangulated " << edgeID << " " << v0 << " " << v1 << " Length " << dl << endl;
        exit(0);
    }

    double vx  = dx/dl;
    double vy  = dy/dl;
    double vz  = dz/dl;

    Array3D pos;
    double offset;

    if( id == 0) {
        offset = junctions[v0].radius;
        pos[0] = junctions[v0].xyz[0] + offset*vx;
        pos[1] = junctions[v0].xyz[1] + offset*vy;
        pos[2] = junctions[v0].xyz[2] + offset*vz;
    }

    if( id == 1) {
        offset = junctions[v1].radius;
        pos[0] = junctions[v1].xyz[0] - offset*vx;
        pos[1] = junctions[v1].xyz[1] - offset*vy;
        pos[2] = junctions[v1].xyz[2] - offset*vz;
    }

    beams[edgeID].profile[id].center = pos;
}

////////////////////////////////////////////////////////////////////////////////

void TriangulateWires :: setProfile(size_t edgeID, size_t nodeID)
{
    auto v0 = min(beams[edgeID].nodes[0], beams[edgeID].nodes[1] );
    auto v1 = max(beams[edgeID].nodes[0], beams[edgeID].nodes[1] );

    int id = 0;
    if( v1 == nodeID ) id = 1;

    const Array3D &p0 = junctions[v0].xyz;
    const Array3D &p1 = junctions[v1].xyz;

    double dx = p1[0] - p0[0];
    double dy = p1[1] - p0[1];
    double dz = p1[2] - p0[2];
    double dl = sqrt(dx*dx + dy*dy + dz*dz);

    if( dl < 1.0E-10) {
        cout << "P0 " << p0[0] << " " << p0[1] << " " << p0[2] << endl;
        cout << "P1 " << p1[0] << " " << p1[1] << " " << p1[2] << endl;
        cout << "Warning: Beam not triangulated " << edgeID << " " << v0 << " " << v1 << " Length " << dl << endl;
        exit(0);
    }

    const Array3D axis = {dx/dl, dy/dl, dz/dl};
    const Array3D srcVec = {0.0, 0.0, 1.0};
    double dtheta = 2.0*M_PI/(double)polySides;

    vector<Array3D> points(polySides);
    double radius;
    // Starting prodile ../ Perturb the coordinates to avoid collinearity  for convex hull
    radius = beams[edgeID].profile[id].radius;
    assert( radius > 1.0E-15);
    for( int i = 0; i < polySides; i++) {
        double t  = i*dtheta + (M_PI/180.0)*rand()/RAND_MAX;
        points[i] = { radius*cos(t), radius*sin(t), 0.0};
    }

    const Array3D &pos = beams[edgeID].profile[id].center;

    vector<Array3D> testPoints(3);
    testPoints[0] = {0.0, 0.0, 0.0};
    testPoints[1] = {1.0, 0.0, 0.0};
    testPoints[2] = {0.0, 1.0, 0.0};
    alignAlong(testPoints, srcVec, axis);
    for( int i = 0; i < 3; i++) {
        testPoints[i][0] += pos[0];
        testPoints[i][1] += pos[1];
        testPoints[i][2] += pos[2];
    }

    auto normal = getTriangleNormal(testPoints[0], testPoints[1], testPoints[2]);
    double t = getVecAngle(normal, axis);

    if( t < 0.5*M_PI) std::reverse(points.begin(), points.end());

    alignAlong(points, srcVec, axis);
    for( int i = 0; i < polySides; i++) {
        points[i][0] += pos[0];
        points[i][1] += pos[1];
        points[i][2] += pos[2];
    }

    assert( isNormalProfile( points, axis) );

    beams[edgeID].profile[id].nodes = points;
}

////////////////////////////////////////////////////////////////////////////////
bool TriangulateWires :: checkIntersections(int nodeID)
{
    auto getLength = [] ( const Array3D &pi, const Array3D &pj) {
        double dx = pj[0] - pi[0];
        double dy = pj[1] - pi[1];
        double dz = pj[2] - pi[2];
        return sqrt(dx*dx + dy*dy + dz*dz);
    };

    int nEdges = junctions[nodeID].adjEdges.size();

    vector<int> edgeid;
    for( auto eid : junctions[nodeID].adjEdges) edgeid.push_back(eid);

    Array3D pi, pj;
    double  ri, rj;
    for( int i = 0; i < nEdges; i++)  {
        int id = edgeid[i];
        int v0 = beams[id].nodes[0];
        int v1 = beams[id].nodes[1];
        if( min(v0,v1) == nodeID) {
            pi = beams[id].profile[0].center;
            ri = beams[id].profile[0].radius;
        }
        if( max(v0,v1) == nodeID) {
            pi = beams[id].profile[1].center;
            ri = beams[id].profile[1].radius;
        }
        for( int j = i+1; j < nEdges; j++) {
            int id = edgeid[j];
            int v0 = beams[id].nodes[0];
            int v1 = beams[id].nodes[1];
            if( min(v0,v1) == nodeID) {
                pj = beams[id].profile[0].center;
                rj = beams[id].profile[0].radius;
            }
            if( max(v0,v1) == nodeID) {
                pj = beams[id].profile[1].center;
                rj = beams[id].profile[1].radius;
            }
            double rij = 1.01*(ri+rj);
            if( getLength(pi,pj) < rij ) {
                return 1;
            }
        }
    }

    return 0;
}

/////////////////////////////////////////////////////////////////////////////

void TriangulateWires :: processNode( size_t nodeID)
{
    const int nEdges = junctions[nodeID].adjEdges.size();

    if( nEdges > 1) {
        setConvexHullRadius(nodeID);
        setConvexHull(nodeID);
        if( reduceNodeVolume ) shrinkNode(nodeID);
    }

    // If the end node of a beam has no neighbours, put cap triangles ...
    if( nEdges == 1) {
        auto iedge = *junctions[nodeID].adjEdges.begin();
        setProfileCenter(iedge, nodeID);
        setProfile(iedge, nodeID);
    }
}

/////////////////////////////////////////////////////////////////////////////

void TriangulateWires :: setConvexHullRadius(size_t nodeID)
{
    const int nEdges = junctions[nodeID].adjEdges.size();
    if( nEdges < 2) return;

    double minlen = std::numeric_limits<double>::max();
    int    pid;

    Array3D ab, cd;

    double maxdist = 0.0;
    double ri, rj;
    for( auto iedge : junctions[nodeID].adjEdges) {
        int va = nodeID;
        int vb = beams[iedge].getOtherNode(nodeID);
        const Array3D &a = junctions[va].xyz;
        const Array3D &b = junctions[vb].xyz;
        ab[0]  = b[0] - a[0];
        ab[1]  = b[1] - a[1];
        ab[2]  = b[2] - a[2];
        minlen = min(minlen, beams[iedge].length);
        pid     = beams[iedge].getProfileID(nodeID);
        ri      = beams[iedge].profile[pid].radius;

        double maxdij = ri;
        for( auto jedge : junctions[nodeID].adjEdges) {
            if( iedge != jedge ) {
                int vc = nodeID;
                int vd = beams[jedge].getOtherNode(nodeID);
                const Array3D &c = junctions[vc].xyz;
                const Array3D &d = junctions[vd].xyz;
                cd[0]  = d[0] - c[0];
                cd[1]  = d[1] - c[1];
                cd[2]  = d[2] - c[2];
                double t  = getVecAngle(ab,cd);
                pid       = beams[jedge].getProfileID(nodeID);
                double rj = beams[jedge].profile[pid].radius;

                if( fabs(t) < 0.5*M_PI) {
                    double dij = (ri + rj*cos(t))/sin(t);
                    maxdij    = max(maxdij, dij);
                }
            }
        }

        pid = beams[iedge].getProfileID(nodeID);
        beams[iedge].profile[pid].minOffset = maxdij;
        maxdist = max(maxdist, maxdij);
    }

    junctions[nodeID].maxRadius  = 0.49*minlen;
    junctions[nodeID].minRadius  = ri;
    junctions[nodeID].radius     = max(ri,maxdist);

    const Array3D &pcenter = junctions[nodeID].xyz;
    double  radius  = junctions[nodeID].radius;
    assert( radius > 0.0);

    while(1) {
        for( auto eid : junctions[nodeID].adjEdges) setProfileCenter(eid, nodeID);

        const int intersects = checkIntersections(nodeID);
        if( intersects == 0) break;
        if( junctions[nodeID].radius < 1.0E-10) {
            cout << "Warning: intersections may be at the node " << nodeID << endl;
            break;
        }
        if( junctions[nodeID].radius > junctions[nodeID].maxRadius ) {
            if( allowBeamThinning) {
                //       cout << "Warning: Beams at the node " <<  nodeID << " intersect" << endl;
                for( auto eid : junctions[nodeID].adjEdges) {
                    const auto v0 = beams[eid].nodes[0];
                    const auto v1 = beams[eid].nodes[1];
                    if( min(v0,v1) == nodeID) beams[eid].profile[0].radius *= 0.90;
                    if( max(v0,v1) == nodeID) beams[eid].profile[1].radius *= 0.90;
                }
            }
            break;
        } else {
            junctions[nodeID].radius *= 1.01;
        }
    }

    for( auto eid : junctions[nodeID].adjEdges) setProfile(eid, nodeID);
}

/////////////////////////////////////////////////////////////////////////////

void TriangulateWires :: refineMesh( TriMesh &nodemesh, size_t nodeID)
{
    int numTris = nodemesh.faces.size();

    vector<Array3D> normals(numTris);

    for( int i = 0; i < numTris; i++) {
        normals[i] = {0.0, 0.0, 0.0};
        for( int j = 0; j < 3; j++) {
            int id = nodemesh.faces[i].nodes[j];
            const Array3D &p = nodemesh.vCoords[id];
            normals[i][0] = p[0];
            normals[i][1] = p[1];
            normals[i][2] = p[2];
        }
        normals[i][0] /= 3.0;
        normals[i][1] /= 3.0;
        normals[i][2] /= 3.0;
    }

    Array3D pnew = nodemesh.vCoords[nodeID];
    double mindist = numeric_limits<double>::max();
    int triid;
    for( int i = 0; i < numTris; i++) {
        double dx = normals[i][0] - pnew[0];
        double dy = normals[i][1] - pnew[1];
        double dz = normals[i][2] - pnew[2];
        double l2 = dx*dx + dy*dy + dz*dz;
        if( l2 < mindist) {
            triid = i;
            mindist = l2;
        }
    }

    nodemesh.faces.erase(nodemesh.faces.begin() + triid);

    size_t v0 = nodemesh.faces[triid].nodes[0];
    size_t v1 = nodemesh.faces[triid].nodes[1];
    size_t v2 = nodemesh.faces[triid].nodes[2];

    Array3I tnodes;
    tnodes = {v0,v1,nodeID};
    nodemesh.faces.emplace_back( tnodes, nodeID, 0);
    tnodes = {v1,v2,nodeID};
    nodemesh.faces.emplace_back( tnodes, nodeID, 0);
    tnodes = {v2,v0,nodeID};
    nodemesh.faces.emplace_back( tnodes, nodeID, 0);
}
/////////////////////////////////////////////////////////////////////////////

int TriangulateWires :: setConvexHull(size_t nodeID)
{
    size_t groupID = nodeID;

    std::vector<int>   nodeMark;
    std::vector<int>   globalNodeID;
    TriMesh nodemesh;

    const Array3D &center = junctions[nodeID].xyz;

    int numEdges = junctions[nodeID].adjEdges.size();

    for( auto eid : junctions[nodeID].adjEdges) {
        const auto v0 = beams[eid].nodes[0];
        const auto v1 = beams[eid].nodes[1];
        if( min(v0,v1) == nodeID) {
            size_t offset = 2*polySides*eid;
            for( int i = 0; i < polySides; i++) {
                const Array3D &p = beams[eid].profile[0].nodes[i];
                nodemesh.vCoords.push_back( {p[0]-center[0], p[1]-center[1], p[2]-center[2]});
                nodeMark.push_back(eid);
                globalNodeID.push_back(offset + i);
            }
        }
        if( max(v0,v1)== nodeID) {
            size_t offset = 2*polySides*eid + polySides;
            for( int i = 0; i < polySides; i++) {
                const Array3D &p = beams[eid].profile[1].nodes[i];
                nodemesh.vCoords.push_back( {p[0]-center[0], p[1]-center[1], p[2]-center[2]} );
                nodeMark.push_back(eid);
                globalNodeID.push_back(offset + i);
            }
        }
    }

    junctions[nodeID].numBoundNodes = nodemesh.vCoords.size();

    size_t numNodes = nodemesh.vCoords.size();
    for( size_t i = 0; i < numNodes; i++) {
        double x = nodemesh.vCoords[i][0];
        double y = nodemesh.vCoords[i][1];
        double z = nodemesh.vCoords[i][2];
        double r = sqrt(x*x + y*y + z*z);
        nodemesh.vCoords[i][0] = x/r;
        nodemesh.vCoords[i][1] = y/r;
        nodemesh.vCoords[i][2] = z/r;
    }

    std::vector<Vec3D> pc(numNodes);
    int index = 0;
    for( auto v : nodemesh.vCoords) pc[index++] = {v[0], v[1], v[2]};

    quickhull::QuickHull<double> qh;
    auto hull = qh.getConvexHull(pc,true,true);
    auto tris = hull.getIndexBuffer();

    Array3I tnodes;
    std::set<int> vSet;

    int numTriangles = tris.size()/3;

    nodemesh.faces.reserve(numEdges*polySides);

    for( size_t i = 0; i < numTriangles; i++) {
        size_t v0 = tris[3*i+0];
        size_t v1 = tris[3*i+1];
        size_t v2 = tris[3*i+2];
        if( (nodeMark[v0] == nodeMark[v1]) && (nodeMark[v1] ==  nodeMark[v2]) ) {
            continue;
        }
        if( v0 == v1 || v1 == v2 || v2 == v0) {
            cout << v0 << " " << v1 << " " << v2 << endl;
            cout << "Warning: A degenerate triangle created at node:  " << nodeID << endl;
            continue;
        }
        vSet.insert(v0);
        vSet.insert(v1);
        vSet.insert(v2);
        tnodes = {v0,v1,v2};
        nodemesh.faces.emplace_back(tnodes, nodeID, 0);
    }

    // Insert the missing nodes ....
    if( vSet.size() != nodemesh.vCoords.size() ) {
        cout << "Warning: Some nodes missed: two more triangles will be added " << endl;
        for( int i = 0; i < numEdges*polySides; i++) {
            if( vSet.find(i) == vSet.end() ) refineMesh( nodemesh, i);
        }
    }

    numTriangles = nodemesh.faces.size();
    junctions[nodeID].faces.resize(numTriangles);

    index = 0;
    for( auto f : nodemesh.faces) {
        auto v0 = f.nodes[0];
        auto v1 = f.nodes[1];
        auto v2 = f.nodes[2];
        tnodes[0] = globalNodeID[ v0 ];
        tnodes[1] = globalNodeID[ v2 ];
        tnodes[2] = globalNodeID[ v1 ];
        junctions[nodeID].faces[index++] = {tnodes, groupID, 0};
    }

    return 0;
}

/////////////////////////////////////////////////////////////////////////////
void TriangulateWires :: triangulate()
{
#ifdef BOOST_TIMER
    boost::timer::auto_cpu_timer watch;
#endif

    if( !edgeRadiusAssigned )  {
        cout << "Warning: Beam radius not specified, estimating using min edgelength " << endl;
        setBeamRadius( 0.10*minEdgeLength );
    }

    int numBeams = beams.size();

    cout << "Building lattice " << endl;
    cout << "#Beams                    : " << numBeams << endl;
    cout << "#MeshNodes  predicted     : " << 2*numBeams*polySides << endl;
    cout << "#MeshTriangles predicted  : " << 4*numBeams*polySides << endl;

    trimesh.vCoords.clear();
    trimesh.faces.clear();

    #pragma omp parallel for
    for( int i = 0; i < numBeams; i++) {
        auto v0 = min(beams[i].nodes[0], beams[i].nodes[1]);
        auto v1 = max(beams[i].nodes[0], beams[i].nodes[1]);
        beams[i].profile[0].center = junctions[v0].xyz;
        beams[i].profile[1].center = junctions[v1].xyz;
    }

    if( make_union ) {
        int numNodes = junctions.size();
        #pragma omp parallel for
        for( int inode = 0; inode < numNodes; inode++) {
            processNode(inode);
        }

    } else {
        #pragma omp parallel for
        for( int i = 0; i < numBeams; i++) {
            setProfile( i, beams[i].nodes[0]);
            setProfile( i, beams[i].nodes[1]);
        }
    }

    assert( trimesh.vCoords.empty() );
    trimesh.vCoords.resize( 2*numBeams*polySides);

    int numEdges = beams.size();
    #pragma omp parallel for
    for( int i = 0; i < numEdges; i++) {
        size_t offset = 2*i*polySides;
        size_t index = 0;
        for( int j = 0; j < 2; j++) {
            for( int k = 0; k < polySides; k++) {
                size_t id = offset + index;
                trimesh.vCoords[id] = beams[i].profile[j].nodes[k];
                index++;
            }
        }
    }

    assert( trimesh.faces.empty() );
    trimesh.faces.reserve(4*numBeams*polySides);
    size_t numNodes = junctions.size();
    for( size_t i = 0; i < numNodes; i++) {
        for( auto f : junctions[i].faces) trimesh.faces.push_back(f);
    }

    buildSideTriangles();

    if( make_union ) {
        #pragma omp parallel for
        for( int inode = 0; inode < junctions.size(); inode++) {
            int nEdges = junctions[inode].adjEdges.size();
            if( nEdges  == 1 ) {
                auto iedge = *junctions[inode].adjEdges.begin();
                buildCapTriangles(iedge, inode);
            }
        }
    } else {
        buildCapTriangles();
    }
}

/////////////////////////////////////////////////////////////////////////////////

void TriangulateWires :: buildSideTriangles()
{
    cout << "Building beam triangles " << endl;
    int numBeams = beams.size();

    size_t foffset = trimesh.faces.size();
    trimesh.faces.resize( foffset + 2*numBeams*polySides);

#ifdef BOOST_TIMER
    boost::timer::auto_cpu_timer watch;
#endif

    Array3I tri;
    #pragma omp parallel for private(tri)
    for( int i = 0; i < numBeams; i++) {
        size_t groupID = junctions.size() + i;
        size_t voffset = 2*i*polySides;
        for( int k = 0; k < polySides; k++) {
            tri[0] = voffset + k;
            tri[1] = voffset + k + polySides;
            tri[2] = voffset + (k+1)%polySides;
            trimesh.faces[foffset + voffset + 2*k] = {tri, groupID, 1};

            tri[0] = voffset + (k+1)%polySides;
            tri[1] = voffset +  k+ polySides;
            tri[2] = voffset + (k+1)%polySides + polySides;
            trimesh.faces[foffset + voffset + 2*k+1] = {tri, groupID, 1};
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void TriangulateWires :: buildCapTriangles(size_t edgeID, size_t nodeID)
{
    size_t offset1 = 2*polySides*edgeID;

    int id = beams[edgeID].getProfileID(nodeID);

    auto groupID = nodeID;

    vector<Face> newfaces(polySides-2);

    Array3I tri;
    if( id == 0) {
        for( int k = 0; k < polySides-2; k++) {
            tri[0] = offset1;
            tri[1] = offset1 + (k+1)%polySides;
            tri[2] = offset1 + (k+2)%polySides;
            newfaces[k] = {tri, groupID, 1};
        }
    }

    if( id == 1) {
        offset1 += polySides;
        for( int k = 0; k < polySides-2; k++) {
            tri[0] = offset1;
            tri[1] = offset1 + (k+2)%polySides;
            tri[2] = offset1 + (k+1)%polySides;
            newfaces[k] = {tri, groupID, 1};
        }
    }

    #pragma omp critical
    for( int k = 0; k < polySides-2; k++)
        trimesh.faces.push_back( newfaces[k] );
}
///////////////////////////////////////////////////////////////////////////////

void TriangulateWires :: buildCapTriangles()
{
    int numEdges = beams.size();
    #pragma omp parallel for
    for( int iedge = 0; iedge < numEdges; iedge++) {
        buildCapTriangles(iedge, beams[iedge].nodes[0]);
        buildCapTriangles(iedge, beams[iedge].nodes[1]);
    }
}

/////////////////////////////////////////////////////////////////////////////

void TriangulateWires :: shrinkNode( size_t nodeID)
{
    double minRadius = junctions[nodeID].minRadius;
    Array3D target;
    for( auto iedge : junctions[nodeID].adjEdges) {
        int  pid    = beams[iedge].getProfileID(nodeID);
        double offset = beams[iedge].profile[pid].minOffset;
        if( offset > minRadius) {
            auto v0 = min(beams[iedge].nodes[0], beams[iedge].nodes[1] );
            auto v1 = max(beams[iedge].nodes[0], beams[iedge].nodes[1] );

            const auto &p0 = junctions[v0].xyz;
            const auto &p1 = junctions[v1].xyz;

            double dx = p1[0] - p0[0];
            double dy = p1[1] - p0[1];
            double dz = p1[2] - p0[2];
            double dl = sqrt(dx*dx + dy*dy + dz*dz);

            double vx  = dx/dl;
            double vy  = dy/dl;
            double vz  = dz/dl;


            if( pid == 0) {
                target[0] = junctions[v0].xyz[0] + offset*vx;
                target[1] = junctions[v0].xyz[1] + offset*vy;
                target[2] = junctions[v0].xyz[2] + offset*vz;
            } else {
                target[0] = junctions[v1].xyz[0] - offset*vx;
                target[1] = junctions[v1].xyz[1] - offset*vy;
                target[2] = junctions[v1].xyz[2] - offset*vz;
            }

            auto center = beams[iedge].profile[pid].center;
            dx = target[0] - center[0];
            dy = target[1] - center[1];
            dz = target[2] - center[2];
            for( size_t i = 0; i < beams[iedge].profile[pid].nodes.size(); i++) {
                beams[iedge].profile[pid].nodes[i][0] += dx;
                beams[iedge].profile[pid].nodes[i][1] += dy;
                beams[iedge].profile[pid].nodes[i][2] += dz;
            }
        }
    }

}
/////////////////////////////////////////////////////////////////////////////


void TriangulateWires :: saveAs( const string &filename)
{

    if( filename.rfind(".off") != string::npos) {

        cout << "Saving triangular mesh in file : " << filename << endl;
        ofstream ofile(filename.c_str(), ios::out);
        ofile << "OFF" << endl;

        ofile << trimesh.vCoords.size() << " " << trimesh.faces.size() << " 0 " << endl;

        for( auto v : trimesh.vCoords) {
            ofile << v[0] << " " << v[1] << "  " << v[2] << endl;
        }

        for( auto f : trimesh.faces) {
            ofile << "3 " << f.nodes[0] << " " << f.nodes[1] << "  " << f.nodes[2] << endl;
        }
    }

    if( filename.rfind(".stl") != string::npos) {
        cout << "Saving triangular mesh in file : " << filename << endl;

        int ntri = trimesh.faces.size();

        ofstream ofile(filename.c_str(), ios::out| ios::binary);
        std::string header = "solid " + filename + "-output";
        char head[80];
        strncpy(head,header.c_str(),sizeof(head)-1);

        ofile.write( head, 80*sizeof(char));
        ofile.write( (char*)&ntri, 4*sizeof(char));

        char attrib[2] = "0";

        float val;
        for( int i = 0; i < ntri; i++) {
            size_t i0 = trimesh.faces[i].nodes[0];
            size_t i1 = trimesh.faces[i].nodes[1];
            size_t i2 = trimesh.faces[i].nodes[2];
            const auto &p0 = trimesh.vCoords[i0];
            const auto &p1 = trimesh.vCoords[i1];
            const auto &p2 = trimesh.vCoords[i2];
            Array3D normal = getTriangleNormal(p0, p1, p2);

            for( int j = 0; j < 3; j++) {
                val = normal[j];
                ofile.write( reinterpret_cast<const char *>(&val), sizeof(float));
            }

            for( int j = 0; j < 3; j++) {
                val = p0[j];
                ofile.write( reinterpret_cast<const char *>(&val), sizeof(float));
            }

            for( int j = 0; j < 3; j++) {
                val = p1[j];
                ofile.write( reinterpret_cast<const char *>(&val), sizeof(float));
            }

            for( int j = 0; j < 3; j++) {
                val = p2[j];
                ofile.write( reinterpret_cast<const char *>(&val), sizeof(float));
            }

            ofile.write( (char*)&attrib, 2);
        }
    }

    /*
       ofstream ofile1("facetag.dat", ios::out);
       ofile1 << trimesh.faceGroup.size() << endl;
       for( size_t i = 0; i < trimesh.faceGroup.size(); i++)
            ofile1 << trimesh.faceEntity[i] << " " << trimesh.faceGroup[i] << endl;
        ofile1.close();

       size_t nedges = 2*beams.size()*polySides;
       ofstream ofile2("profiles.dat", ios::out);
       ofile2 << nedges << endl;
       for( size_t i = 0; i < beams.size(); i++) {
        for( size_t j = 0; j < polySides; j++) {
            auto v0 = beams[i].profile[0].globalNodes[j];
            auto v1 = beams[i].profile[0].globalNodes[(j+1)%polySides];
            ofile2 << v0 << " " << v1 << endl;
        }

        for( size_t j = 0; j < polySides; j++) {
            auto v0 = beams[i].profile[1].globalNodes[j];
            auto v1 = beams[i].profile[1].globalNodes[(j+1)%polySides];
            ofile2 << v0 << " " << v1 << endl;
        }
        }
        ofile2.close();
    */

    ofstream ofile3("facenormal.dat", ios::out);
    ofile3 << trimesh.faces.size() << endl;
    for( int i = 0; i < trimesh.faces.size(); i++) {
        int i0 = trimesh.faces[i].nodes[0];
        int i1 = trimesh.faces[i].nodes[1];
        int i2 = trimesh.faces[i].nodes[2];
        const auto &p0 = trimesh.vCoords[i0];
        const auto &p1 = trimesh.vCoords[i1];
        const auto &p2 = trimesh.vCoords[i2];
        Array3D normal = getTriangleNormal(p0, p1, p2);
        ofile3 << normal[0] << " " << normal[1] << " " << normal[2] << endl;
    }
    ofile3.close();
}

/////////////////////////////////////////////////////////////////////////////
