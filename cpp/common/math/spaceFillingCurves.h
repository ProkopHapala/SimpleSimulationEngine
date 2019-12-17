#ifndef  SpaceFillingCurves_h
#define  SpaceFillingCurves_h

// https://en.wikipedia.org/wiki/Space-filling_tree
// https://en.wikipedia.org/wiki/H_tree
// https://en.wikipedia.org/wiki/Hilbert_curve

// https://commons.wikimedia.org/wiki/File:Mandeltree.svg

// http://people.csail.mit.edu/jaffer/Geometry/HSFC
// http://and-what-happened.blogspot.cz/2011/08/fast-2d-and-3d-hilbert-curves-and.html
// http://blog.notdot.net/2009/11/Damn-Cool-Algorithms-Spatial-indexing-with-Quadtrees-and-Hilbert-Curves

// https://en.wikipedia.org/wiki/Rapidly-exploring_random_tree

namespace HiblertCurve{

inline void rot(int n, int& x, int& y, int rx, int ry) {
    if (ry == 0) {
        if (rx == 1) { x = n-1-x; y = n-1-y; }
        //if (rx == 1) { x = n-1-y; y = n-1-x; }
        int t  = x; x = y; y = t;
    }
}

inline int xy2d (int n, int x, int y) {
    int rx, ry, s, d=0;
    for (s=n>>1; s>0; s>>=1) {
        rx = (x & s) > 0;
        ry = (y & s) > 0;
        d += s * s * ((3 * rx) ^ ry);
        rot(s, x, y, rx, ry);
    }
    return d;
}

inline void d2xy(int n, int d, int& x, int& y) {
    int rx, ry, s, t=d;
    x=y=0;
    for (s=1; s<n; s<<=1) {
        rx = 1 & (t>>1);
        ry = 1 & (t ^ rx);
        rot(s, x, y, rx, ry);
        x += s * rx;
        y += s * ry;
        t >>= 2;
    }
}


/*

//rotate/flip a quadrant appropriately
inline void rot(int n, int *x, int *y, int rx, int ry) {
    if (ry == 0) {
        if (rx == 1) {
            *x = n-1 - *x;
            *y = n-1 - *y;
        }

        //Swap x and y
        int t  = *x;
        *x = *y;
        *y = t;
    }
}

int xy2d (int n, int x, int y) {
    int rx, ry, s, d=0;
    for (s=n/2; s>0; s/=2) {
        rx = (x & s) > 0;
        ry = (y & s) > 0;
        d += s * s * ((3 * rx) ^ ry);
        rot(s, &x, &y, rx, ry);
    }
    return d;
}

//convert d to (x,y)
inline void d2xy(int n, int d, int *x, int *y) {
    int rx, ry, s, t=d;
    *x = *y = 0;
    for (s=1; s<n; s*=2) {
        rx = 1 & (t/2);
        ry = 1 & (t ^ rx);
        rot(s, x, y, rx, ry);
        *x += s * rx;
        *y += s * ry;
        t /= 4;
    }
}

*/

}

#endif
