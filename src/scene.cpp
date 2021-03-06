#include "main.h"
#include "node.h"
#include "noise.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h> //samplers
#include <openvdb/tools/Composite.h> //csg/comp
#include <openvdb/tools/GridOperators.h> //gradient

#include "scene.h"
#include "SceneSurface.h"
#include "SceneDensity.h"

#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/scalable_allocator.h>

#include <cfloat>

namespace Node{

ValueNodeParams::ValueNodeParams(const dfloat3 *_pvoxw, const dfloat3 *_pcptw, float _s, float _p, const dfloat3 *_pvoxwa, float _advdensity, float _advdist, const InputNodeParams *pnp) : pvoxw(_pvoxw), pcptw(_pcptw), distance(_s), density(_p), pvoxwa(_pvoxwa), advdensity(_advdensity), advdist(_advdist), pnodeparams(pnp){
	//
}

ValueNodeParams::~ValueNodeParams(){
	//
}

const dfloat3 * ValueNodeParams::GetObjectPosW() const{
	return &std::get<INP_OBJECT>(*pnodeparams)->location;
}

const dfloat3 * ValueNodeParams::GetVoxPosW() const{
	return pvoxw;
}

const dfloat3 * ValueNodeParams::GetCptPosW() const{
	return pcptw;
}

float ValueNodeParams::GetLocalDistance() const{
	return distance;
}

float ValueNodeParams::GetLocalDensity() const{
	return density;
}

const dfloat3 * ValueNodeParams::GetVoxPosWAdv() const{
	return pvoxwa;
}

float ValueNodeParams::GetAdvectionDistance() const{
	return advdist;
}

float ValueNodeParams::GetAdvectionDensity() const{
	return advdensity;
}

float ValueNodeParams::SampleGlobalDistance(const dfloat3 &p, bool q) const{
	const FloatGridBoxSampler *ptsampler = q?std::get<INP_QGRSAMPLER>(*pnodeparams):std::get<INP_SDFSAMPLER>(*pnodeparams);
	return ptsampler?ptsampler->wsSample(*(openvdb::Vec3f*)&p):FLT_MAX;
}

float ValueNodeParams::SampleGlobalDensity(const dfloat3 &p) const{
	const FloatGridBoxSampler *pfsampler = std::get<INP_FOGSAMPLER>(*pnodeparams);
	return pfsampler?pfsampler->wsSample(*(openvdb::Vec3f*)&p):0.0f;
}

dfloat3 ValueNodeParams::SampleGlobalVector(const dfloat3 &p) const{
	const VectorGridBoxSampler *pvsampler = std::get<INP_VELSAMPLER>(*pnodeparams);
	return pvsampler?*((dfloat3*)pvsampler->wsSample(*(openvdb::Vec3f*)&p).asPointer()):dfloat3(0.0f);
}

dfloat3 ValueNodeParams::SampleGlobalGradient(const dfloat3 &p) const{
	const VectorGridBoxSampler *pgsampler = std::get<INP_GRADSAMPLER>(*pnodeparams);
	return pgsampler?*((dfloat3*)pgsampler->wsSample(*(openvdb::Vec3f*)&p).asPointer()):dfloat3(0.0f);
}

}

BoundingBox::BoundingBox(){
	//
}

BoundingBox::BoundingBox(const float4 &c, const float4 &e){
	float4::store(&sc,c);
	float4::store(&se,e);
}

BoundingBox::~BoundingBox(){
	//
}

//Some intersection code ripped from DirectXCollision library
const uint32_t XM_PERMUTE_0X = 0;
const uint32_t XM_PERMUTE_0Y = 1;
const uint32_t XM_PERMUTE_0Z = 2;
const uint32_t XM_PERMUTE_0W = 3;
const uint32_t XM_PERMUTE_1X = 4;
const uint32_t XM_PERMUTE_1Y = 5;
const uint32_t XM_PERMUTE_1Z = 6;
const uint32_t XM_PERMUTE_1W = 7;

bool BoundingBox::Intersects(const float4 &v0, const float4 &v1, const float4 &v2) const{
	float4 z = float4::zero();
	float4 c = float4::load(&sc);
	float4 e = float4::load(&se);
	float4 bmin = c-e;
	float4 bmax = c+e;
	float4 tmin = float4::min(float4::min(v0,v1),v2);
	float4 tmax = float4::max(float4::max(v0,v1),v2);

	float4 dj = float4::Or(float4::Greater(tmin,bmax),float4::Greater(bmin,tmax));
	if(dj.swizzle<0,1,2,0>().AnyTrue())
		return false;

	float4 n = float4::cross(v1-v0,v2-v0);
	float4 d = float4::dot3(n,v0);

	//if(float4::AllTrue(float4::EqualR(n,z)))
	if(float4::Equal(n,z).AllTrue())
		return false;

	float4 ns = float4::Greater(n,z);
	float4 vmin = float4::select(bmax,bmin,ns);
	float4 vmax = float4::select(bmin,bmax,ns);

	float4 dmin = float4::dot3(vmin,n);
	float4 dmax = float4::dot3(vmax,n);

	float4 ni = float4::Greater(dmin,d);
	ni = float4::Or(ni,float4::Less(dmax,d));

	float4 tv0 = v0-c;
	float4 tv1 = v1-c;
	float4 tv2 = v2-c;
	float4 e0 = tv1-tv0;
	float4 e1 = tv2-tv1;
	float4 e2 = tv0-tv2;

	e0.set<3>(0.0f);
	e1.set<3>(0.0f);
	e2.set<3>(0.0f);

	float4 a, p0, p1, p2, min, max, r;
	//1
	a = float4::permute(e0,-e0,XM_PERMUTE_0W,XM_PERMUTE_1Z,XM_PERMUTE_0Y,XM_PERMUTE_0X);
	p0 = float4::dot3(tv0,a);
	p2 = float4::dot3(tv2,a);
	min = float4::min(p0,p2);
	max = float4::max(p0,p2);
	r = float4::dot3(e,float4::abs(a));
	ni = float4::Or(ni,float4::Greater(min,r));
	ni = float4::Or(ni,float4::Less(max,-r));

	a = float4::permute(e1,-e1,XM_PERMUTE_0W,XM_PERMUTE_1Z,XM_PERMUTE_0Y,XM_PERMUTE_0X);
	p0 = float4::dot3(tv0,a);
	p1 = float4::dot3(tv1,a);
	min = float4::min(p0,p1);
	max = float4::max(p0,p1);
	r = float4::dot3(e,float4::abs(a));
	ni = float4::Or(ni,float4::Greater(min,r));
	ni = float4::Or(ni,float4::Less(max,-r));

	a = float4::permute(e2,-e2,XM_PERMUTE_0W,XM_PERMUTE_1Z,XM_PERMUTE_0Y,XM_PERMUTE_0X);
	p0 = float4::dot3(tv0,a);
	p1 = float4::dot3(tv1,a);
	min = float4::min(p0,p1);
	max = float4::max(p0,p1);
	r = float4::dot3(e,float4::abs(a));
	ni = float4::Or(ni,float4::Greater(min,r));
	ni = float4::Or(ni,float4::Less(max,-r));
	//2
	a = float4::permute(e0,-e0,XM_PERMUTE_0Z,XM_PERMUTE_0W,XM_PERMUTE_1X,XM_PERMUTE_0Y);
	p0 = float4::dot3(tv0,a);
	p2 = float4::dot3(tv2,a);
	min = float4::min(p0,p2);
	max = float4::max(p0,p2);
	r = float4::dot3(e,float4::abs(a));
	ni = float4::Or(ni,float4::Greater(min,r));
	ni = float4::Or(ni,float4::Less(max,-r));

	a = float4::permute(e1,-e1,XM_PERMUTE_0Z,XM_PERMUTE_0W,XM_PERMUTE_1X,XM_PERMUTE_0Y);
	p0 = float4::dot3(tv0,a);
	p1 = float4::dot3(tv1,a);
	min = float4::min(p0,p1);
	max = float4::max(p0,p1);
	r = float4::dot3(e,float4::abs(a));
	ni = float4::Or(ni,float4::Greater(min,r));
	ni = float4::Or(ni,float4::Less(max,-r));

	a = float4::permute(e2,-e2,XM_PERMUTE_0Z,XM_PERMUTE_0W,XM_PERMUTE_1X,XM_PERMUTE_0Y);
	p0 = float4::dot3(tv0,a);
	p1 = float4::dot3(tv1,a);
	min = float4::min(p0,p1);
	max = float4::max(p0,p1);
	r = float4::dot3(e,float4::abs(a));
	ni = float4::Or(ni,float4::Greater(min,r));
	ni = float4::Or(ni,float4::Less(max,-r));
	//3
	a = float4::permute(e0,-e0,XM_PERMUTE_1Y,XM_PERMUTE_0X,XM_PERMUTE_0W,XM_PERMUTE_0Z);
	p0 = float4::dot3(tv0,a);
	p2 = float4::dot3(tv2,a);
	min = float4::min(p0,p2);
	max = float4::max(p0,p2);
	r = float4::dot3(e,float4::abs(a));
	ni = float4::Or(ni,float4::Greater(min,r));
	ni = float4::Or(ni,float4::Less(max,-r));

	a = float4::permute(e1,-e1,XM_PERMUTE_1Y,XM_PERMUTE_0X,XM_PERMUTE_0W,XM_PERMUTE_0Z);
	p0 = float4::dot3(tv0,a);
	p1 = float4::dot3(tv1,a);
	min = float4::min(p0,p1);
	max = float4::max(p0,p1);
	r = float4::dot3(e,float4::abs(a));
	ni = float4::Or(ni,float4::Greater(min,r));
	ni = float4::Or(ni,float4::Less(max,-r));

	a = float4::permute(e2,-e2,XM_PERMUTE_1Y,XM_PERMUTE_0X,XM_PERMUTE_0W,XM_PERMUTE_0Z);
	p0 = float4::dot3(tv0,a);
	p1 = float4::dot3(tv1,a);
	min = float4::min(p0,p1);
	max = float4::max(p0,p1);
	r = float4::dot3(e,float4::abs(a));
	ni = float4::Or(ni,float4::Greater(min,r));
	ni = float4::Or(ni,float4::Less(max,-r));

	return ni.AnyFalse();
}

bool BoundingBox::Intersects(const BoundingBox &bb) const{
	float4 ca = float4::load(&sc);
	float4 ea = float4::load(&se);
	float4 bmina = ca-ea;
	float4 bmaxa = ca+ea;

	float4 cb = float4::load(&bb.sc);
	float4 eb = float4::load(&bb.se);
	float4 bminb = cb-eb;
	float4 bmaxb = cb+eb;

	float4 dj = float4::Or(float4::Greater(bmina,bmaxb),float4::Greater(bminb,bmaxa));
	return dj.swizzle<0,1,2,0>().AllFalse();
}

OctreeStructure::OctreeStructure(){
	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i)
		volx[i] = ~0u;
	memset(chn,0,sizeof(chn));
}

OctreeStructure::~OctreeStructure(){
	//
}

Octree::Octree(uint _x) : x(_x){//, lock(ATOMIC_FLAG_INIT){
	memset(pch,0,sizeof(pch));
}

Octree::~Octree(){
	//
}

void Octree::BuildPath(const float4 &c, const float4 &e, const float4 &c1, const float4 &e1, uint level, uint mlevel, std::atomic<uint> *pindex, std::atomic<uint> *pleafx, tbb::concurrent_vector<Octree> *proot, tbb::concurrent_vector<OctreeStructure> *pob, VOLUME_BUFFER bx){
	//m.lock();
	float4::store(&(*pob)[x].ce,float4::select(c,e,float4::selectctrl(0,0,0,1)));
	memset((*pob)[x].qval,0,sizeof((*pob)[x].qval)); //these are set during the resampling phase
	//m.unlock();

	if(level >= mlevel-1){
		//m.lock();//for(; lock.test_and_set(std::memory_order_acquire););
		if((*pob)[x].volx[bx] == ~0u)
			(*pob)[x].volx[bx] = pleafx->fetch_add(1);
		//m.unlock();//lock.clear(std::memory_order_release);
		return;
	}//else (*pob)[x].volx = ~0;

	BoundingBox aabb1(c1,e1);

	float4 ee = 0.5f*e;
	for(uint i = 0; i < 8; ++i){
		float4 sv = 2.0f*float4((float)(i%2),(float)((i/2)%2),(float)(i/4),0.0f)-float4::one();
		float4 cc = c+sv*ee;
		BoundingBox aabb(cc,1.1f*ee); //scale by 1.1 to expand for cases where surface is goes near the leaf boundary

		//aabb.Intersects
		if(aabb.Intersects(aabb1)){
			//m.lock();//for(; lock.test_and_set(std::memory_order_acquire););
			if(!pch[i]){
				uint index = pindex->fetch_add(1)+1;
				pob->grow_to_at_least(index+1);

				//pch[i] = new(proot+index) Octree(index);
				pch[i] = new(scalable_malloc(sizeof(Octree))) Octree(index);//new(&(*proot)[index]) Octree(index);
				(*pob)[x].chn[i] = index;
			}
			//m.unlock();//lock.clear(std::memory_order_release);

			pch[i]->BuildPath(cc,ee,c1,e1,level+1,mlevel,pindex,pleafx,proot,pob,bx);
		}
	}
}

/*void Octree::BuildPath(const float4 &c, const float4 &e, const float4 &v0, const float4 &v1, const float4 &v2, uint level, uint mlevel, std::atomic<uint> *pindex, std::atomic<uint> *pleafx, tbb::concurrent_vector<Octree> *proot, tbb::concurrent_vector<OctreeStructure> *pob, VOLUME_BUFFER bx){
	//m.lock();
	float4::store(&(*pob)[x].ce,float4::select(c,e,float4::selectctrl(0,0,0,1)));
	memset((*pob)[x].qval,0,sizeof((*pob)[x].qval)); //these are set during the resampling phase
	//m.unlock();

	if(level >= mlevel-1){
		//m.lock();//for(; lock.test_and_set(std::memory_order_acquire););
		if((*pob)[x].volx[bx] == ~0u)
			(*pob)[x].volx[bx] = pleafx->fetch_add(1);
		//m.unlock();//lock.clear(std::memory_order_release);
		return;
	}//else (*pob)[x].volx = ~0;

	float4 ee = 0.5f*e;
	for(uint i = 0; i < 8; ++i){
		float4 sv = 2.0f*float4((float)(i%2),(float)((i/2)%2),(float)(i/4),0.0f)-float4::one();
		float4 cc = c+sv*ee;
		BoundingBox aabb(cc,ee);

		if(aabb.Intersects(v0,v1,v2)){
			//m.lock();//for(; lock.test_and_set(std::memory_order_acquire););
			if(!pch[i]){
				uint index = pindex->fetch_add(1)+1;
				pob->grow_to_at_least(index+1);

				pch[i] = new(scalable_malloc(sizeof(Octree))) Octree(index);//pch[i] = new(&(*proot)[index]) Octree(index);
				(*pob)[x].chn[i] = index;
			}
			//m.unlock();//lock.clear(std::memory_order_release);

			pch[i]->BuildPath(cc,ee,v0,v1,v2,level+1,mlevel,pindex,pleafx,proot,pob,bx);
		}
	}
}*/

void Octree::FreeRecursive(){
	for(uint i = 0; i < 8; ++i)
		if(pch[i]){
			pch[i]->FreeRecursive();
			scalable_free(pch[i]);
		}
}

static bool S_FindSceneInfo(SceneData::BaseObject *pto){
	for(uint j = 0; j < pto->pnt->nodes0.size(); ++j){
		Node::SceneInfo *psci = dynamic_cast<Node::SceneInfo *>(pto->pnt->nodes0[j]);
		if(!psci)
			continue;
		if(psci->omask & ((1<<Node::SceneInfo::OUTPUT_FLOAT_DISTANCE)
			|(1<<(Node::SceneInfo::OUTPUT_FLOAT_COUNT+Node::SceneInfo::OUTPUT_VECTOR_GRADIENT))))
			return true;
	}
	return false;
}

static openvdb::GridBase::Ptr S_ReadGridExcept(openvdb::io::File &vdbc, const char *pname){
	openvdb::GridBase::Ptr pgridb = vdbc.readGrid(pname);
	if(!pgridb){
		vdbc.close();
		throw(0);
	}
	return pgridb;
}

using PostFogParams = std::tuple<SceneData::BaseObject *, openvdb::FloatGrid::Ptr, uint>;
enum PFP{
	PFP_OBJECT,
	PFP_INPUTGRID,
	PFP_FLAGS
};

static void S_Create(float s, float qb, float lvc, float bvc, uint maxd, bool cache, const char *pcachedir, openvdb::FloatGrid::Ptr pgrid[VOLUME_BUFFER_COUNT], Scene *pscene){
	openvdb::math::Transform::Ptr pgridtr = openvdb::math::Transform::createLinearTransform(s);

	//Find if SceneInfo's distance or gradient output was used anywhere in the node trees and automatically determine if a query field should be constructed.
	//TODO: skip holdouts
	bool qfield =
		std::find_if(SceneData::Surface::objs.begin(),SceneData::Surface::objs.end(),S_FindSceneInfo) != SceneData::Surface::objs.end() ||
		std::find_if(SceneData::ParticleSystem::prss.begin(),SceneData::ParticleSystem::prss.end(),S_FindSceneInfo) != SceneData::ParticleSystem::prss.end() ||
		std::find_if(SceneData::SmokeCache::objs.begin(),SceneData::SmokeCache::objs.end(),S_FindSceneInfo) != SceneData::SmokeCache::objs.end();
	if(qfield)
		DebugPrintf("SceneInfo.distance or gradient in use, will construct a query field.\n");

	pgrid[VOLUME_BUFFER_SDF] = openvdb::FloatGrid::create(s*bvc);
	pgrid[VOLUME_BUFFER_SDF]->setTransform(pgridtr);
	pgrid[VOLUME_BUFFER_SDF]->setGridClass(openvdb::GRID_LEVEL_SET);

	pgrid[VOLUME_BUFFER_FOG] = openvdb::FloatGrid::create();
	pgrid[VOLUME_BUFFER_FOG]->setTransform(pgridtr);
	pgrid[VOLUME_BUFFER_FOG]->setGridClass(openvdb::GRID_FOG_VOLUME);

	//Low-res query field
	openvdb::math::Transform::Ptr pqsdftr = openvdb::math::Transform::createLinearTransform(qb/bvc);
	openvdb::FloatGrid::Ptr pqsdf = openvdb::FloatGrid::create(qb);
	pqsdf->setTransform(pqsdftr);
	pqsdf->setGridClass(openvdb::GRID_LEVEL_SET);

	float4 scaabbmin = float4(FLT_MAX);
	float4 scaabbmax = -scaabbmin;

	std::vector<BoundingBox,tbb::cache_aligned_allocator<BoundingBox>> gridbvs[VOLUME_BUFFER_COUNT];
	std::vector<PostFogParams,tbb::cache_aligned_allocator<PostFogParams>> fogppl; //input grids to be post-processed

	for(uint i = 0, n = SceneData::Surface::objs.size(); i < n; ++i){
		if(SceneData::Surface::objs[i]->flags & SCENEOBJ_HOLDOUT)
			continue;

		DebugPrintf("Processing surface %s (%u/%u)\n",SceneData::Surface::objs[i]->pname,i+1,n);

		openvdb::FloatGrid::Ptr ptgrid = 0, phgrid = 0, pdgrid;
		bool qonly = dynamic_cast<Node::OutputNode*>(SceneData::Surface::objs[i]->pnt->GetRoot())->qonly;

		char fn[256];
		snprintf(fn,sizeof(fn),"%s/droplet-surface-cache-%s.vdb",pcachedir,SceneData::Surface::objs[i]->pname);
		openvdb::io::File vdbc(fn);
		try{
			if(!cache || !(SceneData::Surface::objs[i]->flags & SCENEOBJ_CACHED))
				throw(0);
			vdbc.open(false);

			if(!qonly)
				ptgrid = openvdb::gridPtrCast<openvdb::FloatGrid>(S_ReadGridExcept(vdbc,"surface"));
			if(qfield)
				phgrid = openvdb::gridPtrCast<openvdb::FloatGrid>(S_ReadGridExcept(vdbc,"surface.query"));
			//
			pdgrid = openvdb::gridPtrCast<openvdb::FloatGrid>(S_ReadGridExcept(vdbc,"surface.density"));
			if(SceneData::Surface::objs[i]->pnt->GetRoot()->imask & 1<<Node::OutputNode::INPUT_FOGPOST)
				fogppl.push_back(PostFogParams(SceneData::Surface::objs[i],pdgrid->deepCopy(),SceneData::Surface::objs[i]->flags));

			//DebugPrintf("Read cached surface (VDB %f MB)\n",(float)ptgrid->memUsage()/1e6f);
			DebugPrintf("Read cached surface\n");

			vdbc.close();

		}catch(...){
			Node::InputNodeParams snp(SceneData::Surface::objs[i],pgridtr,0,0,0,0,0);
			SceneData::Surface::objs[i]->pnt->EvaluateNodes1(&snp,0,1<<Node::OutputNode::INPUT_SURFACE|1<<Node::OutputNode::INPUT_FOG);

			Node::BaseSurfaceNode1 *pdsn = dynamic_cast<Node::BaseSurfaceNode1*>(SceneData::Surface::objs[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_SURFACE]);

			if(!qonly)
				ptgrid = pdsn->ComputeLevelSet(pgridtr,bvc,bvc);
			if(qfield)
				phgrid = pdsn->ComputeLevelSet(pqsdftr,bvc,bvc);

			pdgrid = dynamic_cast<Node::BaseFogNode1*>(SceneData::Surface::objs[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_FOG])->pdgrid;
			if(SceneData::Surface::objs[i]->pnt->GetRoot()->imask & 1<<Node::OutputNode::INPUT_FOGPOST)
				fogppl.push_back(PostFogParams(SceneData::Surface::objs[i],pdgrid->deepCopy(),SceneData::Surface::objs[i]->flags));

			if(cache){
				if(vdbc.isOpen())
					vdbc.close();
				openvdb::GridCPtrVec gvec;
				if(ptgrid){
					ptgrid->setName("surface");
					gvec.push_back(ptgrid);
				}
				if(phgrid){
					phgrid->setName("surface.query");
					gvec.push_back(phgrid);
				}
				pdgrid->setName("surface.density");
				gvec.push_back(pdgrid);
				vdbc.write(gvec);
				vdbc.close();
			}

			//DebugPrintf("Completed surface calculations (VDB %f MB)\n",(float)ptgrid->memUsage()/1e6f);
			DebugPrintf("Completed surface calculations\n");
		}

		if(ptgrid)
			openvdb::tools::csgUnion(*pgrid[VOLUME_BUFFER_SDF],*ptgrid);
		if(phgrid)
			openvdb::tools::csgUnion(*pqsdf,*phgrid);

		openvdb::tools::compMax(*pgrid[VOLUME_BUFFER_FOG],*pdgrid);
	}

	openvdb::Vec3SGrid::Ptr ptvel = openvdb::Vec3SGrid::create();
	ptvel->setTransform(pgridtr);
	ptvel->setGridClass(openvdb::GRID_FOG_VOLUME);

	for(uint i = 0, n = SceneData::SmokeCache::objs.size(); i < n; ++i){
		openvdb::FloatGrid::Ptr pdgrid;

		DebugPrintf("Processing smoke cache %s (%u/%u)\n",SceneData::SmokeCache::objs[i]->pname,i+1,n);

		char fn[256];
		snprintf(fn,sizeof(fn),"%s/droplet-smoke-cache-%s.vdb",pcachedir,SceneData::SmokeCache::objs[i]->pname);
		openvdb::io::File vdbc(fn);
		try{
			if(!cache || !(SceneData::SmokeCache::objs[i]->flags & SCENEOBJ_CACHED))
				throw(0);
			vdbc.open(false);

			pdgrid = openvdb::gridPtrCast<openvdb::FloatGrid>(S_ReadGridExcept(vdbc,"fog"));
			if(SceneData::SmokeCache::objs[i]->pnt->GetRoot()->imask & 1<<Node::OutputNode::INPUT_FOGPOST)
				fogppl.push_back(PostFogParams(SceneData::SmokeCache::objs[i],pdgrid->deepCopy(),SceneData::SmokeCache::objs[i]->flags));

			DebugPrintf("Read cached smoke cache (VDB %f MB)\n",(float)pdgrid->memUsage()/1e6f);

			vdbc.close();

		}catch(...){
			Node::InputNodeParams snp(SceneData::SmokeCache::objs[i],pgridtr,0,0,0,0,0);
			SceneData::SmokeCache::objs[i]->pnt->EvaluateNodes1(&snp,0,1<<Node::OutputNode::INPUT_FOG);

			pdgrid = dynamic_cast<Node::BaseFogNode1*>(SceneData::SmokeCache::objs[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_FOG])->pdgrid;
			if(SceneData::SmokeCache::objs[i]->pnt->GetRoot()->imask & 1<<Node::OutputNode::INPUT_FOGPOST)
				fogppl.push_back(PostFogParams(SceneData::SmokeCache::objs[i],pdgrid->deepCopy(),0));

			if(cache){
				if(vdbc.isOpen())
					vdbc.close();
				pdgrid->setName("fog");
				openvdb::GridCPtrVec gvec{pdgrid};
				vdbc.write(gvec);
				vdbc.close();
			}

			DebugPrintf("Completed smoke cache calculations (VDB %f MB)\n",(float)pdgrid->memUsage()/1e6f);
		}

		openvdb::tools::compMax(*pgrid[VOLUME_BUFFER_FOG],*pdgrid);
	}

	for(uint i = 0, n = SceneData::ParticleSystem::prss.size(); i < n; ++i){
		openvdb::FloatGrid::Ptr pdgrid;
		openvdb::Vec3SGrid::Ptr pvgrid;

		DebugPrintf("Processing particle fog %s (%u/%u)\n",SceneData::ParticleSystem::prss[i]->pname,i+1,n);

		char fn[256];
		snprintf(fn,sizeof(fn),"%s/droplet-fog-cache-%s.vdb",pcachedir,SceneData::ParticleSystem::prss[i]->pname);
		openvdb::io::File vdbc(fn);
		try{
			if(!cache || !(SceneData::ParticleSystem::prss[i]->flags & SCENEOBJ_CACHED))
				throw(0);
			vdbc.open(false);

			pdgrid = openvdb::gridPtrCast<openvdb::FloatGrid>(S_ReadGridExcept(vdbc,"fog"));
			pvgrid = openvdb::gridPtrCast<openvdb::Vec3SGrid>(S_ReadGridExcept(vdbc,"vel"));
			if(SceneData::ParticleSystem::prss[i]->pnt->GetRoot()->imask & 1<<Node::OutputNode::INPUT_FOGPOST)
				fogppl.push_back(PostFogParams(SceneData::ParticleSystem::prss[i],pdgrid->deepCopy(),SceneData::ParticleSystem::prss[i]->flags));

			DebugPrintf("Read cached particle fog (VDB %f MB)\n",(float)pdgrid->memUsage()/1e6f);

			vdbc.close();

		}catch(...){
			Node::InputNodeParams snp(SceneData::ParticleSystem::prss[i],pgridtr,0,0,0,0,0);
			SceneData::ParticleSystem::prss[i]->pnt->EvaluateNodes1(&snp,0,1<<Node::OutputNode::INPUT_FOG|1<<Node::OutputNode::INPUT_VECTOR);

			pdgrid = dynamic_cast<Node::BaseFogNode1*>(SceneData::ParticleSystem::prss[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_FOG])->pdgrid;
			pvgrid = dynamic_cast<Node::BaseVectorFieldNode1*>(SceneData::ParticleSystem::prss[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_VECTOR])->pvgrid;
			if(SceneData::ParticleSystem::prss[i]->pnt->GetRoot()->imask & 1<<Node::OutputNode::INPUT_FOGPOST)
				fogppl.push_back(PostFogParams(SceneData::ParticleSystem::prss[i],pdgrid->deepCopy(),0));

			if(cache){
				if(vdbc.isOpen())
					vdbc.close();
				pdgrid->setName("fog");
				pvgrid->setName("vel");
				openvdb::GridCPtrVec gvec{pdgrid,pvgrid};
				vdbc.write(gvec);
				vdbc.close();
			}

			DebugPrintf("Completed particle fog calculations (VDB %f MB)\n",(float)pdgrid->memUsage()/1e6f);
		}

		openvdb::tools::compMax(*pgrid[VOLUME_BUFFER_FOG],*pdgrid);
		openvdb::tools::compSum(*ptvel,*pvgrid);
	}

	//fog post-processor
	if(fogppl.size() > 0){
		openvdb::Vec3SGrid::Ptr pgrad = openvdb::tools::gradient(*pqsdf);
		openvdb::FloatGrid::Ptr pdgrid;

		for(uint i = 0, n = fogppl.size(); i < n; ++i){
			DebugPrintf("Post processing fog %s (%u/%u)\n",std::get<PFP_OBJECT>(fogppl[i])->pname,i+1,n);

			char fn[256];
			snprintf(fn,sizeof(fn),"%s/droplet-post-cache-%s.vdb",pcachedir,std::get<PFP_OBJECT>(fogppl[i])->pname);
			openvdb::io::File vdbc(fn);
			try{
				if(!cache || !(std::get<PFP_FLAGS>(fogppl[i]) & SCENEOBJ_CACHED))
					throw(0);
				vdbc.open(false);

				pdgrid = openvdb::gridPtrCast<openvdb::FloatGrid>(S_ReadGridExcept(vdbc,"fog"));
				DebugPrintf("Read post processing results (VDB %f MB)\n",(float)pdgrid->memUsage()/1e6f);

				vdbc.close();

			}catch(...){
				FloatGridBoxSampler *pqsampler = new FloatGridBoxSampler(*pqsdf);
				FloatGridBoxSampler *ppsampler = new FloatGridBoxSampler(*pgrid[VOLUME_BUFFER_FOG]);
				VectorGridBoxSampler *pvsampler = new VectorGridBoxSampler(*ptvel);
				VectorGridBoxSampler *pgsampler = new VectorGridBoxSampler(*pgrad);
				FloatGridBoxSampler *pdsampler = new FloatGridBoxSampler(*pgrid[VOLUME_BUFFER_SDF]);

				SceneData::PostFog fobj(std::get<PFP_OBJECT>(fogppl[i])->pnt,std::get<PFP_INPUTGRID>(fogppl[i]),
					&std::get<PFP_OBJECT>(fogppl[i])->location,std::get<PFP_FLAGS>(fogppl[i]));
				Node::InputNodeParams snp(&fobj,pgridtr,pdsampler,pqsampler,ppsampler,pvsampler,pgsampler);
				fobj.pnt->EvaluateNodes1(&snp,0,1<<Node::OutputNode::INPUT_FOGPOST);

				delete pqsampler;
				delete ppsampler;
				delete pvsampler;
				delete pgsampler;
				delete pdsampler;

				pdgrid = dynamic_cast<Node::BaseFogNode1*>(fobj.pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_FOGPOST])->pdgrid;

				if(cache){
					pdgrid->setName("fog");
					openvdb::GridCPtrVec gvec{pdgrid};
					vdbc.write(gvec);
					vdbc.close();
				}

				DebugPrintf("Completed post processing step (VDB %f MB)\n",(float)pdgrid->memUsage()/1e6f);
			}

			switch(dynamic_cast<Node::OutputNode*>(std::get<PFP_OBJECT>(fogppl[i])->pnt->GetRoot())->opch){
			default:
			case 'M':
				openvdb::tools::compMax(*pgrid[VOLUME_BUFFER_FOG],*pdgrid);
				break;
			case 'm':
				openvdb::tools::compMin(*pgrid[VOLUME_BUFFER_FOG],*pdgrid);
				break;
			case '+':
				openvdb::tools::compSum(*pgrid[VOLUME_BUFFER_FOG],*pdgrid);
				break;
			case '*':
				openvdb::tools::compMul(*pgrid[VOLUME_BUFFER_FOG],*pdgrid);
				break;
			case '=':
				openvdb::tools::compReplace(*pgrid[VOLUME_BUFFER_FOG],*pdgrid);
				break;
			}
		}
	}

	/*
	for(prss){
		t = EvaluateNodes1(FOG); //ParticleInput
		if(fogpost){
			g = DeepCopy(t);
			pplist.push(g);
			combine(a,t);
		}else combine(final,t);
	}

	b = final?combineCopy(final,a):a; //global grid node parameter (there's Tree::combine2(), see cookbook)
	//workaround: b = DeepCopy(final); combine(b,a); //it's okay to lose a-grid
	for(pplist){
		t = EvaluateNodes1(FOGPOST); //FogPostInput
		combine(final,t);
	}

	aabb = EvaluateBounds(final);
	*/

	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i){
		for(openvdb::FloatGrid::TreeType::LeafCIter m = pgrid[i]->tree().cbeginLeaf(); m; ++m){
			const openvdb::FloatGrid::TreeType::LeafNodeType *pl = m.getLeaf();

			openvdb::math::CoordBBox bbox;// = pl->getNodeBoundingBox();
			pl->evalActiveBoundingBox(bbox);

			openvdb::math::Coord bdim = bbox.dim();
			openvdb::Vec3s dimi = bdim.asVec3s(); //index-space dimension
			openvdb::Vec3s extw = 0.5f*dimi*s; //pgrid1->transform().indexToWorld(dimi);

			openvdb::Vec3d posi = bbox.getCenter();
			openvdb::Vec3d posw = pgridtr->indexToWorld(posi);

			BoundingBox aabb;
			aabb.sc = dfloat3(posw.x(),posw.y(),posw.z());
			aabb.se = dfloat3(extw.x(),extw.y(),extw.z());

			gridbvs[i].push_back(aabb);

			float4 c = float4::load(&aabb.sc);
			float4 e = float4::load(&aabb.se);
			scaabbmin = float4::min(c-e,scaabbmin);
			scaabbmax = float4::max(c+e,scaabbmax);
		}
	}

	float4 c = 0.5f*(scaabbmax+scaabbmin);
	float4 e = 0.5f*(scaabbmax-scaabbmin);

	BoundingBox scaabb(c,e);
	float4 a = float4::max(float4::max(e.splat<0>(),e.splat<1>()),e.splat<2>());

	float d = 2.0f*a.get<0>();
	uint mlevel = (uint)ceilf(logf(d/(lvc*s))/SM_LN2);

	//check the depth limit
	if(mlevel > maxd){
		lvc *= powf(2.0f,(float)(mlevel-maxd));
		mlevel = maxd;
	}

	float msdf = (float)pgrid[VOLUME_BUFFER_SDF]->memUsage()/1e6f, mfog = (float)pgrid[VOLUME_BUFFER_FOG]->memUsage()/1e6f;
	DebugPrintf("Center = (%.3f, %.3f, %.3f)\nExtents = (%.3f, %.3f, %.3f) (max w = %.3f)\nVDB total = %f MB\n  SDF = %f MB\n  Fog = %f MB\n",
		scaabb.sc.x,scaabb.sc.y,scaabb.sc.z,scaabb.se.x,scaabb.se.y,scaabb.se.z,d,msdf+mfog,msdf,mfog);

	float y = powf(2.0f,(float)mlevel);
	float r = y*lvc; //sparse voxel resolution
	//The actual voxel size is now v = d/(2^k*N), where k is the octree depth and N=lvc.
	DebugPrintf("> Constructing octree (depth = %u, leaf = %u*%f = %f, sparse res = %u^3)...\n",mlevel,(uint)lvc,d/r,lvc*d/r,(uint)r);

	Octree *proot = new(scalable_malloc(sizeof(Octree))) Octree(0);

	pscene->ob.reserve(500000);
	pscene->ob.emplace_back();

	std::atomic<uint> indexa(0);
	std::array<std::atomic<uint>,VOLUME_BUFFER_COUNT> leafxn = {};

#if 0
	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i){
		tbb::parallel_for(tbb::blocked_range<size_t>(0,gridbvs[i].size()),[&](const tbb::blocked_range<size_t> &nr){
			for(uint j = nr.begin(); j < nr.end(); ++j){
				float4 c1 = float4::load(&gridbvs[i][j].sc);
				float4 e1 = float4::load(&gridbvs[i][j].se);
				proot->BuildPath(c,a,c1,e1,0,mlevel,&indexa,&leafxb,&pscene->root,&pscene->ob,VOLUME_BUFFER_FOG);
			}
		});
	}
#else
	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i){
		for(uint j = 0, n = gridbvs[i].size(); j < n; ++j){
			float4 c1 = float4::load(&gridbvs[i][j].sc);
			float4 e1 = float4::load(&gridbvs[i][j].se);
			proot->BuildPath(c,a,c1,e1,0,mlevel,&indexa,&leafxn[i],&pscene->root,&pscene->ob,(VOLUME_BUFFER)i);
		}
	}
#endif

	proot->FreeRecursive();
	scalable_free(proot);

	pscene->lvoxc = (uint)lvc;
	pscene->index = indexa;
	pscene->leafx[VOLUME_BUFFER_SDF] = leafxn[VOLUME_BUFFER_SDF];
	pscene->leafx[VOLUME_BUFFER_FOG] = leafxn[VOLUME_BUFFER_FOG];
}

namespace SceneData{
#define STRDUP(s) strcpy(new char[strlen(s)+1],s)

BaseObject::BaseObject(Node::NodeTree *_pnt, const char *_pname, const dfloat3 *ploc, uint _flags) : pnt(_pnt), location(*ploc), flags(_flags){
	pname = STRDUP(_pname);
}

BaseObject::~BaseObject(){
	delete []pname;
}

ParticleSystem::ParticleSystem(Node::NodeTree *_pnt, const char *pname, const dfloat3 *ploc, uint flags) : BaseObject(_pnt,pname,ploc,flags){
	ParticleSystem::prss.push_back(this);
}

ParticleSystem::~ParticleSystem(){
	//
}

void ParticleSystem::DeleteAll(){
	for(uint i = 0; i < prss.size(); ++i)
		delete prss[i];
	prss.clear();
}

std::vector<ParticleSystem *> ParticleSystem::prss;

SmokeCache::SmokeCache(Node::NodeTree *_pnt, const char *pname, const dfloat3 *ploc, uint flags, const char *_pvdb, const char *_prho, const char *_pvel) : BaseObject(_pnt,pname,ploc,flags){
	SmokeCache::objs.push_back(this);
	pvdb = STRDUP(_pvdb);
	prho = STRDUP(_prho);
	pvel = STRDUP(_pvel);
}

SmokeCache::~SmokeCache(){
	delete []pvdb;
	delete []prho;
	delete []pvel;
}

void SmokeCache::DeleteAll(){
	for(uint i = 0; i < objs.size(); ++i)
		delete objs[i];
	objs.clear();
}

std::vector<SmokeCache *> SmokeCache::objs;

Surface::Surface(Node::NodeTree *_pnt, const char *pname, const dfloat3 *ploc, uint flags) : BaseObject(_pnt,pname,ploc,flags){
	Surface::objs.push_back(this);
}

Surface::~Surface(){
	//
}

void Surface::DeleteAll(){
	for(uint i = 0; i < objs.size(); ++i)
		delete objs[i];
	objs.clear();
}

std::vector<Surface *> Surface::objs;

}

#undef STRDUP

Scene::Scene(){
	//
}

Scene::~Scene(){
	//
}

void Scene::Initialize(float s, uint maxd, float qb, uint smask, bool cache, const char *pcachedir){
	openvdb::initialize();

	const float lvc = 8.0f; //minimum number of voxels in an octree leaf
	const float bvc = 4.0f; //number of narrow band voxels counting from the surface

	openvdb::FloatGrid::Ptr pgrid[VOLUME_BUFFER_COUNT];// = {0};
	FloatGridBoxSampler *psampler[VOLUME_BUFFER_COUNT];

	S_Create(s,qb,lvc,bvc,maxd,cache,pcachedir,pgrid,this);

	DebugPrintf("> Resampling volume data...\n");

	try{
		lvoxc3 = lvoxc*lvoxc*lvoxc;
		for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i){
			psampler[i] = new FloatGridBoxSampler(*pgrid[i]); //non-cached, thread safe version
			pvol[i] = new float[lvoxc3*leafx[i]];
		}
	}catch(std::bad_alloc &ba){
		DebugPrintf("FATAL: bad allocation: %s\n",ba.what());
	}

	const uint uN = lvoxc;//BLCLOUD_uN;
	float4 nv = float4((float)uN);

	/*typedef openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor, openvdb::tools::BoxSampler> FastGridSampler;
	tbb::enumerable_thread_specific<FastGridSampler> fsampler([&]()->FastGridSampler{
		return FastGridSampler(pgrid->getConstAccessor(),pgrid->transform()); //Seems to cause severe problems sometimes.
	});*/
	tbb::parallel_for(tbb::blocked_range<size_t>(0,index),[&](const tbb::blocked_range<size_t> &nr){
		//FastGridSampler &ffs = fsampler.local();
		openvdb::Vec3f posw;
		for(uint i = nr.begin(); i < nr.end(); ++i){
			if(ob[i].volx[VOLUME_BUFFER_SDF] == ~0u){
				if(ob[i].volx[VOLUME_BUFFER_FOG] == ~0u)
					continue; //not a leaf; exit early
				float d = psampler[VOLUME_BUFFER_SDF]->wsSample(*(openvdb::Vec3f*)&ob[i].ce);
				if(d < 0.0f){
					//If the fog leaf is completely inside the sdf surface (no overlapping sdf leaf -> volx == ~0u),
					//remove it. It's useless there, and removing it simplifies the space skipping algorithm.
					ob[i].volx[VOLUME_BUFFER_FOG] = ~0u;
					continue;
				}
			}
			//
			float4 nc = float4::load(&ob[i].ce);
			float4 ne = nc.splat<3>();
			//
			ob[i].qval[VOLUME_BUFFER_SDF] = FLT_MAX;
			ob[i].qval[VOLUME_BUFFER_FOG] = 0.01f;
			//
			for(uint j = 0; j < uN*uN*uN; ++j){
				float4 nn = (nv-2.0f*float4((float)(j%uN),(float)((j/uN)%uN),(float)(j/(uN*uN)),1.0f)-float4::one())/(nv-float4::one());
				float4 nw = nc-ne*nn;
				float4::store((dfloat3*)posw.asPointer(),nw);

				if(ob[i].volx[VOLUME_BUFFER_SDF] != ~0u){
					pvol[VOLUME_BUFFER_SDF][ob[i].volx[VOLUME_BUFFER_SDF]*lvoxc3+j] = psampler[VOLUME_BUFFER_SDF]->wsSample(posw);
					ob[i].qval[VOLUME_BUFFER_SDF] = openvdb::math::Min(ob[i].qval[VOLUME_BUFFER_SDF],pvol[VOLUME_BUFFER_SDF][ob[i].volx[VOLUME_BUFFER_SDF]*lvoxc3+j]);
				}

				if(ob[i].volx[VOLUME_BUFFER_FOG] != ~0u){
					pvol[VOLUME_BUFFER_FOG][ob[i].volx[VOLUME_BUFFER_FOG]*lvoxc3+j] = openvdb::math::Max(openvdb::math::Min(psampler[VOLUME_BUFFER_FOG]->wsSample(posw),1.0f),0.0f);
					ob[i].qval[VOLUME_BUFFER_FOG] = openvdb::math::Max(ob[i].qval[VOLUME_BUFFER_FOG],pvol[VOLUME_BUFFER_FOG][ob[i].volx[VOLUME_BUFFER_FOG]*lvoxc3+j]);
				}
			}
		}
	});

	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i)
		delete psampler[i];

	float msdf = (float)(leafx[VOLUME_BUFFER_SDF]*lvoxc3*sizeof(float))/1e6f;
	float mfog = (float)(leafx[VOLUME_BUFFER_FOG]*lvoxc3*sizeof(float))/1e6f;
	DebugPrintf("Volume size = %f MB\n  SDF = %f MB\n  Fog = %f MB\n",msdf+mfog,msdf,mfog);
	//
}

void Scene::Destroy(){
	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i)
		delete []pvol[i];
	ob.clear();
}
