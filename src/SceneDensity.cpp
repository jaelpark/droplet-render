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

BaseVectorFieldNode1::BaseVectorFieldNode1(uint _level, NodeTree *pnt) : BaseVectorFieldNode(_level,pnt){
	pvgrid = openvdb::Vec3SGrid::create();
	pvgrid->setGridClass(openvdb::GRID_FOG_VOLUME); //this probably doesn't matter here
}

BaseVectorFieldNode1::~BaseVectorFieldNode1(){
	//
}

BaseVectorFieldNode * BaseVectorFieldNode::Create(uint level, NodeTree *pnt){
	return new BaseVectorFieldNode1(level,pnt);
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

	pntree->EvaluateNodes0(0,level+1,emask);

	openvdb::math::Transform::Ptr pgridtr = std::get<INP_TRANSFORM>(*pd);
	pdgrid->setTransform(pgridtr);

	openvdb::FloatGrid::Ptr ptgrid;
	if(pgridtr->voxelSize().x() < prasres->result.local()){
		openvdb::math::Transform::Ptr pgridtr1 = openvdb::math::Transform::createLinearTransform(prasres->result.local());
		ptgrid = openvdb::FloatGrid::create();
        ptgrid->setGridClass(openvdb::GRID_FOG_VOLUME);
		ptgrid->setTransform(pgridtr1);
	}else ptgrid = pdgrid;

	DebugPrintf("> Rasterizing particles...\n");

	openvdb::FloatGrid::Accessor grida = ptgrid->getAccessor();
	for(uint i = 0; i < pps->vl.size(); ++i){
		//pntree->EvaluateNodes0(0,level+1,emask);
		openvdb::Vec3d t(pps->vl[i].x,pps->vl[i].y,pps->vl[i].z);
		openvdb::Vec3d c = ptgrid->transform().worldToIndex(t); //assume cell-centered indices
		openvdb::Vec3f f = openvdb::Vec3f(floorf(c.x()-0.5f),floorf(c.y()-0.5f),floorf(c.z()-0.5f));
		openvdb::Vec3f b = c-f;
		openvdb::Vec3f B = openvdb::Vec3f(1.0f)-b;

		//TODO: particle weight factor
		openvdb::Coord q((int)f.x(),(int)f.y(),(int)f.z());
		grida.modifyValue(q.offsetBy(0,0,0),[&](float &v){v += pweight->result.local()*B.x()*B.y()*B.z();});
		grida.modifyValue(q.offsetBy(1,0,0),[&](float &v){v += pweight->result.local()*b.x()*B.y()*B.z();});
		grida.modifyValue(q.offsetBy(0,1,0),[&](float &v){v += pweight->result.local()*B.x()*b.y()*B.z();});
		grida.modifyValue(q.offsetBy(1,1,0),[&](float &v){v += pweight->result.local()*b.x()*b.y()*B.z();});
		grida.modifyValue(q.offsetBy(0,0,1),[&](float &v){v += pweight->result.local()*B.x()*B.y()*b.z();});
		grida.modifyValue(q.offsetBy(1,0,1),[&](float &v){v += pweight->result.local()*b.x()*B.y()*b.z();});
		grida.modifyValue(q.offsetBy(0,1,1),[&](float &v){v += pweight->result.local()*B.x()*b.y()*b.z();});
		grida.modifyValue(q.offsetBy(1,1,1),[&](float &v){v += pweight->result.local()*b.x()*b.y()*b.z();});
	}

	if(pgridtr->voxelSize().x() < prasres->result.local()){
		DebugPrintf("> Upsampling particle fog...\n");
		openvdb::tools::resampleToMatch<openvdb::tools::BoxSampler>(*ptgrid,*pdgrid);
	}else DebugPrintf("Used native grid resolution for particle rasterization.\n");

	//advection: trace voxel if density < threshold

	/*
	options:
	-include the velocity with the fog:
		-need bool input to enable computation
		-hard to separate fog/velocity - e.g. for need to use different velocity field than the fog provides
	-separate ParticleInput (ParticleVelocity) node
		-use the VectorField base class
	*/

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
	InputNodeParams<SceneObject> *pd = (InputNodeParams<SceneObject>*)pp;

	BaseValueNode<float> *pthrs = dynamic_cast<BaseValueNode<float>*>(pnodes[IAdvection::INPUT_THRESHOLD]);
	BaseValueNode<float> *pdist = dynamic_cast<BaseValueNode<float>*>(pnodes[IAdvection::INPUT_DISTANCE]);
	BaseValueNode<int> *piters = dynamic_cast<BaseValueNode<int>*>(pnodes[IAdvection::INPUT_ITERATIONS]);
	BaseVectorFieldNode1 *pvfn = dynamic_cast<BaseVectorFieldNode1*>(pnodes[IAdvection::INPUT_VELOCITY]);
	BaseFogNode1 *pnode = dynamic_cast<BaseFogNode1*>(pnodes[IAdvection::INPUT_FOG]);

	pntree->EvaluateNodes0(0,level+1,emask);

	openvdb::math::Transform::Ptr pgridtr = std::get<INP_TRANSFORM>(*pd);

	openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> samplerd(*pnode->pdgrid);
	openvdb::tools::GridSampler<openvdb::Vec3SGrid, openvdb::tools::StaggeredBoxSampler> samplerv(*pvfn->pvgrid);

	DebugPrintf("> Advecting fog volume...\n");

	typedef std::tuple<openvdb::FloatGrid::Ptr, openvdb::FloatGrid::Accessor> FloatGridT;
    tbb::enumerable_thread_specific<FloatGridT> tgrida([&]()->FloatGridT{
        openvdb::FloatGrid::Ptr ptgrid = openvdb::FloatGrid::create();
        ptgrid->setTransform(pgridtr);
        ptgrid->setGridClass(openvdb::GRID_FOG_VOLUME);
        return FloatGridT(ptgrid,ptgrid->getAccessor());
    });
    tbb::parallel_for(openvdb::tree::IteratorRange<openvdb::FloatGrid::ValueOnIter>(pnode->pdgrid->beginValueOn()),[&](openvdb::tree::IteratorRange<openvdb::FloatGrid::ValueOnIter> &r){
        FloatGridT &fgt = tgrida.local();
        for(; r; ++r){
			pntree->EvaluateNodes0(0,level+1,emask);
			
            const openvdb::FloatGrid::ValueOnIter &m = r.iterator();

			float f = m.getValue();
			if(f > pthrs->result.local())
				continue;

			openvdb::Coord c = m.getCoord();
			openvdb::math::Vec3s posw = pnode->pdgrid->transform().indexToWorld(c);

			float s = pdist->result.local()/(float)piters->result.local();

			float4 rc = float4::load(posw.asPointer());
			for(uint i = 0; i < piters->result.local(); ++i){
				openvdb::math::Vec3s v = samplerv.wsSample(posw);
				rc += s*float4::load((dfloat3*)v.asPointer());
				//
				//float4::store((dfloat3*)posw.asPointer(),rc);
				//float p = samplerd.wsSample(posw); //TODO: do actual integration instead of sampling the last value?

				//-when sampling noise velocity, do it directly somehow, not by first rendering it to a grid
				//-use virtual Sample() for the BaseVectorField nodes, to either sample grid or return noise
				//-just need to figure out how the actual sampler is created and stored
				//Calling Sample() will evaluate the nodes. Need multithreading support.
			}

			float4::store((dfloat3*)posw.asPointer(),rc);
			float p = samplerd.wsSample(posw);

			std::get<1>(fgt).setValue(c,f+p);
        }
    });

	//TODO: try compSum
	openvdb::FloatGrid::Accessor dgrida = pdgrid->getAccessor();
    for(tbb::enumerable_thread_specific<FloatGridT>::const_iterator q = tgrida.begin(); q != tgrida.end(); ++q){
        //
        for(openvdb::FloatGrid::ValueOnIter m = std::get<0>(*q)->beginValueOn(); m.test(); ++m){
			openvdb::Coord c = m.getCoord();
            float f = m.getValue();
            dgrida.setValue(c,f);
        }
    }
}

Node::IAdvection * IAdvection::Create(uint level, NodeTree *pnt){
	return new Advection(level,pnt);
}

}
