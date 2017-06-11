#ifndef SCENE_DENSITY_H
#define SCENE_DENSITY_H

namespace SceneData{

class PostFog : public BaseObject{
public:
	PostFog(Node::NodeTree *, openvdb::FloatGrid::Ptr, const dfloat3 *, uint);
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
	virtual void Clear();
	openvdb::FloatGrid::Ptr pdgrid;
};

class BaseVectorFieldNode1 : public virtual BaseVectorFieldNode{
public:
	BaseVectorFieldNode1(uint, NodeTree *);
	~BaseVectorFieldNode1();
	virtual void Clear();
	openvdb::Vec3SGrid::Ptr pvgrid;
};

class ParticleInput : public BaseFogNode1, public IParticleInput{
public:
	ParticleInput(uint, NodeTree *);
	~ParticleInput();
	void Evaluate(const void *);
};

class FieldInput : public BaseFogNode1, public BaseVectorFieldNode1, public IFieldInput{
public:
	FieldInput(uint, NodeTree *);
	~FieldInput();
	void Evaluate(const void *);
	void Clear();
};

class SmokeCache : public BaseFogNode1, public BaseVectorFieldNode1, public ISmokeCache{
public:
	SmokeCache(uint, NodeTree *);
	~SmokeCache();
	void Evaluate(const void *);
	void Clear();
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

class Combine : public BaseFogNode1, public ICombine{
public:
	Combine(uint, NodeTree *, char);
	~Combine();
	void Evaluate(const void *);
	char opch;
};

class Advection : public BaseFogNode1, public IAdvection{
public:
	Advection(uint, NodeTree *, uint);
	~Advection();
	void Evaluate(const void *);
	uint flags;
};

}

#endif
