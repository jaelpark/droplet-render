#include "main.h"
#include "node.h"
#include "noise.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/ParticlesToLevelSet.h>
#include <openvdb/tools/LevelSetUtil.h> //sdfToFogVolume
#include <openvdb/tools/GridTransformer.h> //resampleToMatch
#include <openvdb/tools/Composite.h>

#include "scene.h"
#include "SceneDensity.h"

namespace SceneData{

PostFog::PostFog(Node::NodeTree *_pnt, openvdb::FloatGrid::Ptr _pdgrid) : BaseObject(_pnt), pdgrid(_pdgrid){
	//
}

PostFog::~PostFog(){
	//
}

}

namespace Node{

BaseFogNode1::BaseFogNode1(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseNode(_level,pnt){
	pdgrid = openvdb::FloatGrid::create();
	pdgrid->setGridClass(openvdb::GRID_FOG_VOLUME);
}

BaseFogNode1::~BaseFogNode1(){
	//
}

BaseFogNode * BaseFogNode::Create(uint level, NodeTree *pnt){
	return new BaseFogNode1(level,pnt);
}

BaseVectorFieldNode1::BaseVectorFieldNode1(uint _level, NodeTree *pnt) : BaseVectorFieldNode(_level,pnt), BaseNode(_level,pnt){
	pvgrid = openvdb::Vec3SGrid::create();
	pvgrid->setGridClass(openvdb::GRID_FOG_VOLUME);
}

BaseVectorFieldNode1::~BaseVectorFieldNode1(){
	//
}

BaseVectorFieldNode * BaseVectorFieldNode::Create(uint level, NodeTree *pnt){
	return new BaseVectorFieldNode1(level,pnt);
}

class ParticleInputList{
public:
	typedef openvdb::Vec3R PosType;
	ParticleInputList(SceneData::ParticleSystem *_pps, float _rscale) : pps(_pps), rscale(_rscale){}//, rscale(1.0f), vscale(1.0f){}
	//ParticleInputList(openvdb::Real rscale1 = 1.0f, openvdb::Real vscale1 = 1.0f) : rscale(rscale1), vscale(vscale1){}

	size_t size() const{
		return pps->pl.size();
	}

	void getPos(size_t x, PosType &p) const{
		dfloat3 &t = pps->pl[x];
		p = openvdb::Vec3R(t.x,t.y,t.z);
	}

	void getPosRad(size_t x, PosType &p, openvdb::Real &r) const{
		dfloat3 &t = pps->pl[x];
		p = openvdb::Vec3R(t.x,t.y,t.z);
		r = rscale;
	}

protected:
	SceneData::ParticleSystem *pps;
	float rscale;
};

ParticleInput::ParticleInput(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseFogNode1(_level,pnt), BaseNode(_level,pnt), IParticleInput(_level,pnt){
	//
	//DebugPrintf(">> ParticleInput()\n");
}

ParticleInput::~ParticleInput(){
	//
}

void ParticleInput::Evaluate(const void *pp){
	InputNodeParams *pd = (InputNodeParams*)pp;
	SceneData::BaseObject *pbob = std::get<INP_OBJECT>(*pd);
	SceneData::ParticleSystem *pps = dynamic_cast<SceneData::ParticleSystem*>(pbob);
	if(!pps){
		DebugPrintf("Warning: Invalid use of ParticleInput where %s is expected. Input forced to empty.",typeid(pbob).name());
		return;
	}

	BaseValueNode<float> *psizen = dynamic_cast<BaseValueNode<float>*>(pnodes[INPUT_SIZE]); //should lifetime determine the size?
	BaseValueNode<float> *pcoffn = dynamic_cast<BaseValueNode<float>*>(pnodes[INPUT_CUTOFF]);

	dfloat3 zr(0.0f);
	ValueNodeParams np(&zr,&zr,0.0f,0.0f,&zr,0.0f,0.0f,pd);
	pntree->EvaluateNodes0(&np,level+1,emask);

	float size = psizen->locr(indices[INPUT_SIZE]); //TODO: these should be static. Remove EvaluateNodes0 above.
	float coff = pcoffn->locr(indices[INPUT_CUTOFF]);

	openvdb::math::Transform::Ptr pgridtr = std::get<INP_TRANSFORM>(*pd);
	//pdgrid->setTransform(pgridtr);

	pdgrid = openvdb::FloatGrid::create(coff);
	pdgrid->setTransform(pgridtr);
	pdgrid->setGridClass(openvdb::GRID_LEVEL_SET);

	ParticleInputList pl(pps,size);

	DebugPrintf("> Rasterizing particle spheres...\n");
	openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> lsf(*pdgrid);
	lsf.setRmax(1e5f); //max radius in voxel units
	lsf.setGrainSize(1);
	lsf.rasterizeSpheres(pl);
	lsf.finalize();

	float r = size/pgridtr->voxelSize().x();
	if(r < lsf.getRmin() || r > lsf.getRmax())
		DebugPrintf("Warning: Particle size either to big or small for sphere rasterizer."); //TODO: check this earlier somewhere

	DebugPrintf("> Converting fog volume...\n");
	openvdb::tools::sdfToFogVolume(*pdgrid);

	pdgrid->tree().voxelizeActiveTiles(true); //voxelize the interior so that futher processing works
}

IParticleInput * IParticleInput::Create(uint level, NodeTree *pnt){
	return new ParticleInput(level,pnt);
}

FieldInput::FieldInput(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseFogNode1(_level,pnt), BaseVectorFieldNode(_level,pnt), BaseVectorFieldNode1(_level,pnt), BaseNode(_level,pnt), IFieldInput(_level,pnt){
	//
	//printf("FieldInput()\n");
}

FieldInput::~FieldInput(){
	//
}

void FieldInput::Evaluate(const void *pp){
	InputNodeParams *pd = (InputNodeParams*)pp;
	SceneData::BaseObject *pbob = std::get<INP_OBJECT>(*pd);
	SceneData::ParticleSystem *pps = dynamic_cast<SceneData::ParticleSystem*>(pbob);
	if(!pps){
		DebugPrintf("Warning: Invalid use of FieldInput where %s is expected. Input forced to empty.",typeid(pbob).name());
		return;
	}

	BaseValueNode<float> *prasres = dynamic_cast<BaseValueNode<float>*>(pnodes[INPUT_RASTERIZATIONRES]);
	BaseValueNode<float> *pweight = dynamic_cast<BaseValueNode<float>*>(pnodes[INPUT_WEIGHT]);

	dfloat3 zr(0.0f);
	ValueNodeParams np(&zr,&zr,0.0f,0.0f,&zr,0.0,0.0f,pd);
	pntree->EvaluateNodes0(&np,level+1,emask);

	openvdb::math::Transform::Ptr pgridtr = std::get<INP_TRANSFORM>(*pd);
	pdgrid->setTransform(pgridtr);
	pvgrid->setTransform(pgridtr);

	openvdb::FloatGrid::Ptr ptgridd;
	openvdb::Vec3SGrid::Ptr ptgridv;
	if(pgridtr->voxelSize().x() < prasres->locr(indices[INPUT_RASTERIZATIONRES])){ //TODO: const expression
		openvdb::math::Transform::Ptr pgridtr1 = openvdb::math::Transform::createLinearTransform(prasres->locr(indices[INPUT_RASTERIZATIONRES]));

		ptgridd = openvdb::FloatGrid::create();
		ptgridd->setGridClass(openvdb::GRID_FOG_VOLUME);
		ptgridd->setTransform(pgridtr1);

		ptgridv = openvdb::Vec3SGrid::create();
		ptgridv->setGridClass(openvdb::GRID_FOG_VOLUME);
		ptgridv->setTransform(pgridtr1);

	}else{
		ptgridd = pdgrid;
		ptgridv = pvgrid;
	}

	DebugPrintf("> Rasterizing particle fields...\n");

	//TODO: Multithreading maybe? Currently it's not too slow though
	openvdb::FloatGrid::Accessor dgrida = ptgridd->getAccessor();
	openvdb::Vec3SGrid::Accessor vgrida = ptgridv->getAccessor();
	for(uint i = 0; i < pps->pl.size(); ++i){
		openvdb::Vec3s posw(pps->pl[i].x,pps->pl[i].y,pps->pl[i].z);
		openvdb::Vec3s c = ptgridd->transform().worldToIndex(posw); //assume cell-centered indices
		openvdb::Vec3f f = openvdb::Vec3f(floorf(c.x()-0.5f),floorf(c.y()-0.5f),floorf(c.z()-0.5f));
		openvdb::Vec3f b = c-f;
		openvdb::Vec3f B = openvdb::Vec3f(1.0f)-b;

		ValueNodeParams np1((dfloat3*)posw.asPointer(),&zr,0.0f,0.0f,(dfloat3*)posw.asPointer(),0.0f,0.0f,pd);
		pntree->EvaluateNodes0(&np1,level+1,emask);

		openvdb::Coord q((int)f.x(),(int)f.y(),(int)f.z());

		{
			dgrida.modifyValue(q.offsetBy(0,0,0),[&](float &v){v += pweight->locr(indices[INPUT_WEIGHT])*B.x()*B.y()*B.z();});
			dgrida.modifyValue(q.offsetBy(1,0,0),[&](float &v){v += pweight->locr(indices[INPUT_WEIGHT])*b.x()*B.y()*B.z();});
			dgrida.modifyValue(q.offsetBy(0,1,0),[&](float &v){v += pweight->locr(indices[INPUT_WEIGHT])*B.x()*b.y()*B.z();});
			dgrida.modifyValue(q.offsetBy(1,1,0),[&](float &v){v += pweight->locr(indices[INPUT_WEIGHT])*b.x()*b.y()*B.z();});
			dgrida.modifyValue(q.offsetBy(0,0,1),[&](float &v){v += pweight->locr(indices[INPUT_WEIGHT])*B.x()*B.y()*b.z();});
			dgrida.modifyValue(q.offsetBy(1,0,1),[&](float &v){v += pweight->locr(indices[INPUT_WEIGHT])*b.x()*B.y()*b.z();});
			dgrida.modifyValue(q.offsetBy(0,1,1),[&](float &v){v += pweight->locr(indices[INPUT_WEIGHT])*B.x()*b.y()*b.z();});
			dgrida.modifyValue(q.offsetBy(1,1,1),[&](float &v){v += pweight->locr(indices[INPUT_WEIGHT])*b.x()*b.y()*b.z();});
		}

		if(omask & 0x2){
			openvdb::Vec3s velw(pps->vl[i].x,pps->vl[i].y,pps->vl[i].z);
			vgrida.modifyValue(q.offsetBy(0,0,0),[&](openvdb::Vec3s &v){v += velw*B.x()*B.y()*B.z();});
			vgrida.modifyValue(q.offsetBy(1,0,0),[&](openvdb::Vec3s &v){v += velw*b.x()*B.y()*B.z();});
			vgrida.modifyValue(q.offsetBy(0,1,0),[&](openvdb::Vec3s &v){v += velw*B.x()*b.y()*B.z();});
			vgrida.modifyValue(q.offsetBy(1,1,0),[&](openvdb::Vec3s &v){v += velw*b.x()*b.y()*B.z();});
			vgrida.modifyValue(q.offsetBy(0,0,1),[&](openvdb::Vec3s &v){v += velw*B.x()*B.y()*b.z();});
			vgrida.modifyValue(q.offsetBy(1,0,1),[&](openvdb::Vec3s &v){v += velw*b.x()*B.y()*b.z();});
			vgrida.modifyValue(q.offsetBy(0,1,1),[&](openvdb::Vec3s &v){v += velw*B.x()*b.y()*b.z();});
			vgrida.modifyValue(q.offsetBy(1,1,1),[&](openvdb::Vec3s &v){v += velw*b.x()*b.y()*b.z();});
		}
	}

	//normalize the vgrid by dgrid (weighted average)
	if(omask & 0x2){
		for(openvdb::Vec3SGrid::ValueOnIter m = ptgridv->beginValueOn(); m.test(); ++m){
			float w = dgrida.getValue(m.getCoord())/pweight->locr(indices[INPUT_WEIGHT]);
			m.setValue(*m/w);
		}
	}

	if(pgridtr->voxelSize().x() < prasres->locr(indices[INPUT_RASTERIZATIONRES])){
		DebugPrintf("> Upsampling particle fog...\n");
		if(omask & 0x1)
			openvdb::tools::resampleToMatch<openvdb::tools::BoxSampler>(*ptgridd,*pdgrid);
		if(omask & 0x2)
			openvdb::tools::resampleToMatch<openvdb::tools::BoxSampler>(*ptgridv,*pvgrid);
	}else DebugPrintf("Used native grid resolution for particle rasterization.\n");

	if(omask & 0x1)
		pdgrid->tree().prune();
	if(omask & 0x2)
		pvgrid->tree().prune();
}

IFieldInput * IFieldInput::Create(uint level, NodeTree *pnt){
	return new FieldInput(level,pnt);
}

SmokeCache::SmokeCache(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseVectorFieldNode(_level,pnt), BaseFogNode1(_level,pnt), BaseVectorFieldNode1(_level,pnt), BaseNode(_level,pnt), ISmokeCache(_level,pnt){
	//
	//DebugPrintf(">> ParticleInput()\n");
}

SmokeCache::~SmokeCache(){
	//
}

void SmokeCache::Evaluate(const void *pp){
	InputNodeParams *pd = (InputNodeParams*)pp;
	SceneData::BaseObject *pbob = std::get<INP_OBJECT>(*pd);
	SceneData::SmokeCache *psmc = dynamic_cast<SceneData::SmokeCache*>(pbob);
	if(!psmc){
		DebugPrintf("Warning: Invalid use of ParticleInput where %s is expected. Input forced to empty.",typeid(pbob).name());
		return;
	}

	//Experimental Blender VDB smoke cache support. As long as there's no decent way to get cache name references with Python, this won't work.
	//TODO: might want to get the velocity grid as well
	//openvdb::io::File vdbc("/tmp/blendcache_droplet_random1/fog_000070_00.vdb");
	openvdb::io::File vdbc(psmc->pvdb);
	try{
		vdbc.open(false);
		if(strlen(psmc->prho) > 0){
			openvdb::FloatGrid::Ptr ptgrid = openvdb::gridPtrCast<openvdb::FloatGrid>(vdbc.readGrid(psmc->prho));
			DebugPrintf("Read OpenVDB smoke cache: %s\n",ptgrid->getName().c_str());

			DebugPrintf("> Upsampling smoke cache...\n");
			openvdb::math::Transform::Ptr pgridtr = std::get<INP_TRANSFORM>(*pd);
			pdgrid->setTransform(pgridtr);

			openvdb::tools::resampleToMatch<openvdb::tools::BoxSampler>(*ptgrid,*pdgrid);
			pdgrid->tree().prune();
		}

		if(strlen(psmc->pvel) > 0){
			openvdb::Vec3SGrid::Ptr ptgrid = openvdb::gridPtrCast<openvdb::Vec3SGrid>(vdbc.readGrid(psmc->pvel));
			DebugPrintf("Read OpenVDB smoke cache: %s\n",ptgrid->getName().c_str());

			DebugPrintf("> Upsampling smoke cache...\n");
			openvdb::math::Transform::Ptr pgridtr = std::get<INP_TRANSFORM>(*pd);
			pvgrid->setTransform(pgridtr);

			openvdb::tools::resampleToMatch<openvdb::tools::BoxSampler>(*ptgrid,*pvgrid);
			pvgrid->tree().prune();
		}

		vdbc.close();

	}catch(const openvdb::IoError &e){
		DebugPrintf("SmokeCache: %s\n",e.what());
	}
}

ISmokeCache * ISmokeCache::Create(uint level, NodeTree *pnt){
	return new SmokeCache(level,pnt);
}

FogPostInput::FogPostInput(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseFogNode1(_level,pnt), BaseNode(_level,pnt), IFogPostInput(_level,pnt){
	//
}

FogPostInput::~FogPostInput(){
	//
}

void FogPostInput::Evaluate(const void *pp){
	InputNodeParams *pd = (InputNodeParams*)pp;
	SceneData::BaseObject *pbob = std::get<INP_OBJECT>(*pd);
	SceneData::PostFog *ppf = dynamic_cast<SceneData::PostFog*>(pbob);
	if(!ppf){
		DebugPrintf("Warning: Invalid use of PostFog where %s is expected. Input forced to empty.",typeid(pbob).name());
		return;
	}
	pdgrid = ppf->pdgrid;
}

IFogPostInput * IFogPostInput::Create(uint level, NodeTree *pnt){
	return new FogPostInput(level,pnt);
}

Composite::Composite(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseFogNode1(_level,pnt), BaseNode(_level,pnt), IComposite(_level,pnt){
	//
}

Composite::~Composite(){
	//
}

void Composite::Evaluate(const void *pp){
	InputNodeParams *pd = (InputNodeParams*)pp;

	BaseValueNode<float> *pvalue = dynamic_cast<BaseValueNode<float>*>(pnodes[IComposite::INPUT_VALUE]);
	BaseFogNode1 *pnode = dynamic_cast<BaseFogNode1*>(pnodes[IComposite::INPUT_FOG]);

	dfloat3 zr(0.0f);
	const FloatGridBoxSampler *psamplers[] = {std::get<INP_SDFSAMPLER>(*pd),std::get<INP_FOGSAMPLER>(*pd)};

	openvdb::math::Transform::Ptr pgridtr = std::get<INP_TRANSFORM>(*pd);
	pdgrid->setTransform(pgridtr);

	DebugPrintf("> Compositing fog volume...\n");

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
			const openvdb::FloatGrid::ValueOnIter &m = r.iterator();

			openvdb::Coord c = m.getCoord();
			openvdb::math::Vec3s posw = pnode->pdgrid->transform().indexToWorld(c); //should use the same pgridtr as here

			ValueNodeParams np1((dfloat3*)posw.asPointer(),&zr,0.0f,m.getValue(),(dfloat3*)posw.asPointer(),m.getValue(),0.0f,pd);
			pntree->EvaluateNodes0(&np1,level+1,emask);

			float f = pvalue->locr(indices[IComposite::INPUT_VALUE]);
			std::get<1>(fgt).setValue(c,f);
		}
	});

	for(tbb::enumerable_thread_specific<FloatGridT>::const_iterator q = tgrida.begin(); q != tgrida.end(); ++q)
		openvdb::tools::compSum(*pdgrid,*std::get<0>(*q));
}

IComposite * IComposite::Create(uint level, NodeTree *pnt){
	return new Composite(level,pnt);
}

Advection::Advection(uint _level, NodeTree *pnt, uint _flags) : BaseFogNode(_level,pnt), BaseFogNode1(_level,pnt), BaseNode(_level,pnt), IAdvection(_level,pnt), flags(_flags){
	//
}

Advection::~Advection(){
	//
}

void Advection::Evaluate(const void *pp){
	InputNodeParams *pd = (InputNodeParams*)pp;

	BaseValueNode<float> *pthrs = dynamic_cast<BaseValueNode<float>*>(pnodes[INPUT_THRESHOLD]);
	BaseValueNode<float> *pdist = dynamic_cast<BaseValueNode<float>*>(pnodes[INPUT_DISTANCE]);
	BaseValueNode<int> *piters = dynamic_cast<BaseValueNode<int>*>(pnodes[INPUT_ITERATIONS]);
	BaseValueNode<float> *pdn = dynamic_cast<BaseValueNode<float>*>(pnodes[INPUT_DENSITY]);
	BaseValueNode<dfloat3> *pvn = dynamic_cast<BaseValueNode<dfloat3>*>(pnodes[INPUT_VELOCITY]);
	BaseFogNode1 *pnode = dynamic_cast<BaseFogNode1*>(pnodes[INPUT_FOG]);

	dfloat3 zr(0.0f);

	openvdb::math::Transform::Ptr pgridtr = std::get<INP_TRANSFORM>(*pd);
	pdgrid->setTransform(pgridtr);

	DebugPrintf("> Advecting fog volume...\n");

	FloatGridBoxSampler sampler(*pnode->pdgrid); //advection sampler

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
			const openvdb::FloatGrid::ValueOnIter &m = r.iterator();

			openvdb::Coord c = m.getCoord();
			openvdb::math::Vec3s posw = pnode->pdgrid->transform().indexToWorld(c);

			float f = m.getValue();
			float4 rc = float4::load(posw.asPointer());

			ValueNodeParams np1((dfloat3*)posw.asPointer(),(dfloat3*)posw.asPointer(),0.0f,f,(dfloat3*)posw.asPointer(),f,0.0f,pd);
			pntree->EvaluateNodes0(&np1,level+1,emask);

			float th = pthrs->locr(indices[INPUT_THRESHOLD]);
			if(f > th){
				std::get<1>(fgt).setValue(c,0.0f); //std::get<1>(fgt).setValue(c,f);
				continue;
			}

			float s = pdist->locr(indices[INPUT_DISTANCE])/(float)piters->locr(indices[INPUT_ITERATIONS]); //step size
			float p; //density

			uint ic = piters->locr(indices[INPUT_ITERATIONS]);
			for(uint i = 0; i < ic; ++i){
				p = pdn->locr(indices[INPUT_DENSITY]);
				//TODO: do actual integration instead of sampling the last value?
				if(p > th && flags & 1<<BOOL_BREAK_ITERATION)
					break;
				dfloat3 vs = pvn->locr(indices[INPUT_VELOCITY]);
				float4 v = float4::load(&vs);
				if(float4::dot3(v,v).get<0>() < 1e-5f)
					break;
				rc += s*v;
				//

				openvdb::math::Vec3s poswa;
				float4::store((dfloat3*)poswa.asPointer(),rc);

				//if advection info / density used
				float f1 = sampler.wsSample(poswa);

				new(&np1) ValueNodeParams((dfloat3*)posw.asPointer(),&zr,0.0f,f,(dfloat3*)poswa.asPointer(),f1,(float)(i+1)/(float)ic,pd);
				pntree->EvaluateNodes0(&np1,level+1,emask);
			}

			//std::get<1>(fgt).setValue(c,f+p);
			std::get<1>(fgt).setValue(c,p);
		}
	});

	for(tbb::enumerable_thread_specific<FloatGridT>::const_iterator q = tgrida.begin(); q != tgrida.end(); ++q)
		openvdb::tools::compSum(*pdgrid,*std::get<0>(*q));
}

IAdvection * IAdvection::Create(uint level, NodeTree *pnt, uint flags){
	return new Advection(level,pnt,flags);
}

}
