#ifndef Rival_h
#define Rival_h

struct Weather{
    float cloudDensity;
    float cloudHeight; // if 0 it is fog
    float Temperature;
    float snow;
    float mud;
}

struct Terrain{
    float coverGroups[3]; // cover of different size 1-soldier, 2-gun, 3-tank
}

struct Battlefield{
    Weather weather;
    Terrain terrain;
};

#endif
