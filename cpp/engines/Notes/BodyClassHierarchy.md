
### Problem:

We need hierarchy of following objects:

* ``StaticPointBody      {pos}``
* ``DynamicPointBody {pos,vpos}`` 
* ``StaticRigidBody      {pos,rot}`` 
* ``DynamicRigidBody {pos,rot,vpos,vrot}`` 
* ``CollisionPointBody  {pos,R}``
* ``CollisionRigidBody  {pos,dims}``


```C++
class DynamicPointBody{
	Vec3d pos;
	Vec3d vel;
};

class DynamicRigidBody{
	//Vec3d pos;   // from DynamicPointBody
	//Vec3d vel;   // from  DynamicPointBody
	Mat3d rot;
	Vec3d vrot;
};

class StaticRigidBody{
	Vec3d pos;
	Mat3d rot;
};
```

Solution leads to diamond problem:

```C++
class StaticPointBody { // dummy class
	Vec3d pos;
};

class StaticRigidBody : public StaticPointBody {
	//Vec3d pos;  // from StaticPointBody
	Mat3d rot;
};

class DynamicPointBody : public StaticPointBody {
	//Vec3d pos;  // from StaticPointBody
	Vec3d vel;
};

class DynamicRigidBody : public DynamicPointBody {
	//Vec3d pos;  // from StaticPointBody
	//Vec3d vel;   // from  DynamicPointBody
	//Mat3d rot;   // from  DynamicPointBody
	Vec3d vrot;
};
```

Solution with contain:

```C++
class StaticPointBody { // dummy class
	double R;
	Vec3d pos;
};

class StaticRigidBody : public StaticPointBody {
	//double R;   // from StaticPointBody
	//Vec3d pos;  // from StaticPointBody
	Mat3d rot;
	Vec3d dims;
};

class DynamicPointBody : public StaticPointBody {
	//double R;   // from StaticPointBody
	//Vec3d pos;  // from StaticPointBody
	double mass;
	Vec3d vel;
};

class DynamicRigidBody : public DynamicPointBody, public StaticRigidBody {
	//double R;    // from StaticPointBody
	//Vec3d pos;   // from StaticPointBody
	//Mat3d rot;   // from StaticRigidBody
	//Vec3d dims;  // from StaticRigidBody
	//double mass; // from DynamicPointBody
	//Vec3d vel;   // from DynamicPointBody
	Vec3d vrot;
	Vec3d I;
};
```

without diamond

```C++


class _Rotation {
	Mat3d rot;
};

class DynamicPointBody : public StaticPointBody {
	Vec3d pos;  // from PointBody
	double mass;
	Vec3d vel;
};

class DynamicRigidBody : public DynamicPointBody, public _Rotation {
	//Vec3d pos;   // from PointBody
	//Mat3d rot;   // from _Rotation 
	//Vec3d vel;   // from PointBody
	Vec3d vrot;
	Vec3d I;
};
```



```C++
class ObjectID{
    int id,kind,shape;
};

class _Position{
    double R;
    Vec3d pos;
};

class _Rotation{
    Mat3d rot;
    Vec3d dims;
};

// ============= Bodies

class PointBody : public _Position{
    Vec3d  pos;
    double mass;
};

class RigidBody : public PointBody, public _Rotation {
    Vec3d vrot;
    Vec3d rotMass;
};

class StaticObject : public _Position, public _Rotation {
    _Position* controlPos;
    _Rotation* controlRot;
};

class RigidBodyID : public RigidBody, public ObjectID {};

class StaticObjectID : public StaticObject, public ObjectID {};
```
