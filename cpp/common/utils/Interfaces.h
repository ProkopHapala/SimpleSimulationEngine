#ifndef  Interfaces_h
#define  Interfaces_h

#include "Vec2.h"
#include "Vec3.h"

class AnyControler{ public:
    //virtual void control(void* obj, double dt)=0;
    virtual void update(double dt)=0;
};

template<typename T>
class Solver{ public:
    //virtual void control(void* obj, double dt)=0;
    virtual void solve(int n, T* x, T* b)=0;
};

// Interf
/**
 * @brief The BindLoader class is an abstract base class that defines the interface for objects that can bind another object and load data from it
 */
class BindLoader {
public:
    /**
     * @brief Binds and loads data from specified object.
     * @param o A pointer to the object to bind and load data from
     * @return An integer indicating the result of the operation.
     */
    virtual int bindLoad(void* o) = 0;
};

/**
 * @brief The Picker class is an abstract base class that defines the interface for picking operations.
 */
class Picker {
public:
    /**
     * @brief Picks the nearest object along the given ray.
     * 
     * @param ray0 The starting point of the ray.
     * @param hray The direction and length of the ray.
     * @param ipick The index of the picked object.
     * @param mask The mask specifying which objects to consider for picking.
     * @return The type of the picked object.
     */
    virtual int pick_nearest(Vec3d ray0, Vec3d hray, int& ipick, int mask = 0xFFFFFFFF, double Rmax=0.0 ) = 0;

    /**
     * @brief Picks all objects intersected by the given ray.
     * 
     * @param ray0 The starting point of the ray.
     * @param hray The direction and length of the ray.
     * @param out  The output array of picked objects. 
     * @param mask The mask specifying which objects to consider for picking.
     * @return The number of objects added to the picking list.
     */
    virtual int pick_all(Vec3d ray0, Vec3d hray, int* out, int mask = 0xFFFFFFFF, double Rmax=0.0 ) = 0;

    /**
     * @brief Gets the object associated with the specified index in the picking list.
     * 
     * @param picked The index of the picked object.
     * @param mask   The mask specifying which objects to consider for picking.
     * @return       A pointer to the picked object.
     */
    virtual void* getPickedObject(int picked, int mask = 0xFFFFFFFF ) = 0;

};

#endif
