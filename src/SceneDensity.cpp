#include "main.h"
#include "node.h"
#include "scene.h"
#include "noise.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/GridTransformer.h>
#include <openvdb/tools/Composite.h>

#include "SceneDensity.h"

namespace Node{

template<class T>
using InputNodeParams = std::tuple<T *, openvdb::math::Transform::Ptr>;
enum INP{
	INP_OBJECT,
	INP_TRANSFORM
};

BaseFogNode1::BaseFogNode1(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt){
    pdgrid = openvdb::FloatGrid::create();
    pdgrid->setGridClass(openvdb::GRID_FOG_VOLUME);
}

BaseFogNode1::~BaseFogNode1(){
    //
}

BaseFogNode * BaseFogNode::Create(uint level, NodeTree *pnt){
    return new BaseFogNode1(level,pnt);
}

ParticleInput::ParticleInput(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseFogNode1(_level,pnt), IParticleInput(_level,pnt){
    //
    //DebugPrintf(">> ParticleInput()\n");
}

ParticleInput::~ParticleInput(){
    //
}

void ParticleInput::Evaluate(const void *pp){
	InputNodeParams<ParticleSystem> *pd = (InputNodeParams<ParticleSystem>*)pp;
	ParticleSystem *pps = std::get<INP_OBJECT>(*pd);

	BaseValueNode<float> *prasres = dynamic_cast<BaseValueNode<float>*>(pnodes[IParticleInput::INPUT_RASTERIZATIONRES]);
	BaseValueNode<float> *pweight = dynamic_cast<BaseValueNode<float>*>(pnodes[IParticleInput::INPUT_WEIGHT]);

	openvdb::math::Transform::Ptr pgridtr = std::get<INP_TRANSFORM>(*pd);
	pdgrid->setTransform(pgridtr);

	openvdb::FloatGrid::Ptr ptgrid;
	if(pgridtr->voxelSize().x() < prasres->result){
		openvdb::math::Transform::Ptr pgridtr1 = openvdb::math::Transform::createLinearTransform(prasres->result);
		ptgrid = openvdb::FloatGrid::create();
        ptgrid->setGridClass(openvdb::GRID_FOG_VOLUME);
		ptgrid->setTransform(pgridtr1);
	}else ptgrid = pdgrid;

	DebugPrintf("> Rasterizing particles...\n");

	openvdb::FloatGrid::Accessor grida = ptgrid->getAccessor();
	for(uint i = 0; i < pps->vl.size(); ++i){
		openvdb::Vec3d t(pps->vl[i].x,pps->vl[i].y,pps->vl[i].z);
		openvdb::Vec3d c = ptgrid->transform().worldToIndex(t); //assume cell-centered indices
		openvdb::Vec3f f = openvdb::Vec3f(floorf(c.x()-0.5f),floorf(c.y()-0.5f),floorf(c.z()-0.5f));
		openvdb::Vec3f b = c-f;
		openvdb::Vec3f B = openvdb::Vec3f(1.0f)-b;

		//TODO: particle weight factor
		openvdb::Coord q((int)f.x(),(int)f.y(),(int)f.z());
		grida.modifyValue(q.offsetBy(0,0,0),[&](float &v){v += pweight->result*B.x()*B.y()*B.z();});
		grida.modifyValue(q.offsetBy(1,0,0),[&](float &v){v += pweight->result*b.x()*B.y()*B.z();});
		grida.modifyValue(q.offsetBy(0,1,0),[&](float &v){v += pweight->result*B.x()*b.y()*B.z();});
		grida.modifyValue(q.offsetBy(1,1,0),[&](float &v){v += pweight->result*b.x()*b.y()*B.z();});
		grida.modifyValue(q.offsetBy(0,0,1),[&](float &v){v += pweight->result*B.x()*B.y()*b.z();});
		grida.modifyValue(q.offsetBy(1,0,1),[&](float &v){v += pweight->result*b.x()*B.y()*b.z();});
		grida.modifyValue(q.offsetBy(0,1,1),[&](float &v){v += pweight->result*B.x()*b.y()*b.z();});
		grida.modifyValue(q.offsetBy(1,1,1),[&](float &v){v += pweight->result*b.x()*b.y()*b.z();});
	}

	if(pgridtr->voxelSize().x() < prasres->result){
		DebugPrintf("> Upsampling particle fog...\n");
		openvdb::tools::resampleToMatch<openvdb::tools::BoxSampler>(*ptgrid,*pdgrid);
	}else DebugPrintf("Used native grid resolution for particle rasterization.\n");

	//advection: trace voxel if density < threshold

	pdgrid->tree().prune();
}

Node::IParticleInput * IParticleInput::Create(uint level, NodeTree *pnt){
	return new ParticleInput(level,pnt);
}

Advection::Advection(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseFogNode1(_level,pnt), IAdvection(_level,pnt){
	//
}

Advection::~Advection(){
	//
}

void Advection::Evaluate(const void *pp){
	//
}

Node::IAdvection * IAdvection::Create(uint level, NodeTree *pnt){
	return new Advection(level,pnt);
}

}
