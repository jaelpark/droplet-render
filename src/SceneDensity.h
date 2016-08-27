#ifndef SCENE_DENSITY_H
#define SCENE_DENSITY_H

namespace SceneData{

class PostFog : public BaseObject{
public:
	PostFog(Node::NodeTree *, openvdb::FloatGrid::Ptr);
	~PostFog();
	openvdb::FloatGrid::Ptr pdgrid; //input grid
};

}

namespace Node{

//additional layer of abstraction - the openvdb compile times are ridiculous
class BaseFogNode1 : public virtual BaseFogNode{
public:
    BaseFogNode1(uint, NodeTree *);
    ~BaseFogNode1();
    openvdb::FloatGrid::Ptr pdgrid;
	//openvdb::Vec3SGrid::Ptr pvgrid; //sadly, there's no 4d vector grid
};

class BaseVectorFieldNode1 : public virtual BaseVectorFieldNode{
public:
	BaseVectorFieldNode1(uint, NodeTree *);
	~BaseVectorFieldNode1();
	openvdb::Vec3SGrid::Ptr pvgrid;
};

class ParticleInput : public BaseFogNode1, public IParticleInput{
public:
	ParticleInput(uint, NodeTree *);
    ~ParticleInput();
    void Evaluate(const void *);
};

class SmokeCache : public BaseFogNode1, public ISmokeCache{
public:
	SmokeCache(uint, NodeTree *);
	~SmokeCache();
	void Evaluate(const void *);
};

class FogPostInput : public BaseFogNode1, public IFogPostInput{
public:
	FogPostInput(uint, NodeTree *);
	~FogPostInput();
	void Evaluate(const void *);
};

class Composite : public BaseFogNode1, public IComposite{
public:
	Composite(uint, NodeTree *);
	~Composite();
	void Evaluate(const void *);
};

class Advection : public BaseFogNode1, public IAdvection{
public:
	Advection(uint, NodeTree *);
	~Advection();
	void Evaluate(const void *);
};

class VectorFieldSampler : public IVectorFieldSampler{
public:
	VectorFieldSampler(uint, NodeTree *);
	~VectorFieldSampler();
	void Evaluate(const void *);
};

}

#endif
