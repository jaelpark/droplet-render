#ifndef SCENE_SURFACE_H
#define SCENE_SURFACE_H

namespace Node{

//additional layer of abstraction - the openvdb compile times are ridiculous
class BaseSurfaceNode1 : public virtual BaseSurfaceNode{
public:
    BaseSurfaceNode1(uint, NodeTree *);
    ~BaseSurfaceNode1();
	openvdb::FloatGrid::Ptr ComputeLevelSet(openvdb::math::Transform::Ptr, float) const;
    openvdb::FloatGrid::Ptr pbgrid; //billowing grid
    std::vector<openvdb::Vec3s> vl;
    std::vector<openvdb::Vec3I> tl;
    std::vector<openvdb::Vec4I> ql;
};

class SurfaceInput : public BaseSurfaceNode1, public ISurfaceInput{
public:
    SurfaceInput(uint, NodeTree *);
    ~SurfaceInput();
    void Evaluate(const void *);
};

class Displacement : public BaseSurfaceNode1, public IDisplacement{
public:
    Displacement(uint _level, NodeTree *pnt, float);
    ~Displacement();
    void Evaluate(const void *);
	float resf;
};

}

#endif
