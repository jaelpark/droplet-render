#ifndef SCENE_DENSITY_H
#define SCENE_DENSITY_H

namespace Node{

//additional layer of abstraction - the openvdb compile times are ridiculous
class BaseFogNode1 : public virtual BaseFogNode{
public:
    BaseFogNode1(uint, NodeTree *);
    ~BaseFogNode1();
    openvdb::FloatGrid::Ptr pdgrid;
	openvdb::Vec3SGrid::Ptr pvgrid; //sadly, there's no 4d vector grid
};

class ParticleInput : public BaseFogNode1, public IParticleInput{
public:
	ParticleInput(uint, NodeTree *);
    ~ParticleInput();
    void Evaluate(const void *);
};

class Advection : public BaseFogNode1, public IAdvection{
public:
	Advection(uint, NodeTree *);
	~Advection();
	void Evaluate(const void *);
};

}

#endif
