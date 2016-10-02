#include "main.h"
#include "node.h"
#include "noise.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h> //samplers
#include <openvdb/tools/Composite.h> //csg/comp

#include "scene.h"
#include "SceneSurface.h"
#include "SceneDensity.h"

#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/scalable_allocator.h>

namespace Node{

ValueNodeParams::ValueNodeParams(const dfloat3 *_pvoxw, const dfloat3 *_pcptw, float _s, float _p, const FloatGridBoxSampler *_psampler[VOLUME_BUFFER_COUNT], const FloatGridBoxSampler *_pqsampler, const VectorGridBoxSampler *_pvsampler)
	: pvoxw(_pvoxw), pcptw(_pcptw), distance(_s), density(_p){
	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i)
		psampler[i] = _psampler[i];
	pqsampler = _pqsampler;
	pvsampler = _pvsampler;
}

ValueNodeParams::~ValueNodeParams(){
	//
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

float ValueNodeParams::SampleGlobalDistance(const dfloat3 &p, bool q) const{
	const FloatGridBoxSampler *ptsampler = q?pqsampler:psampler[VOLUME_BUFFER_SDF];
	return ptsampler?ptsampler->wsSample(*(openvdb::Vec3f*)&p):FLT_MAX;
}

float ValueNodeParams::SampleGlobalDensity(const dfloat3 &p) const{
	return psampler[VOLUME_BUFFER_FOG]?psampler[VOLUME_BUFFER_FOG]->wsSample(*(openvdb::Vec3f*)&p):0.0f;
}

dfloat3 ValueNodeParams::SampleGlobalVector(const dfloat3 &p) const{
	return pvsampler?*((dfloat3*)pvsampler->wsSample(*(openvdb::Vec3f*)&p).asPointer()):dfloat3(0.0f);
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
	float4 z = sfloat1::zero();
	float4 c = float4::load(&sc);
	float4 e = float4::load(&se);
	float4 bmin = c-e;
	float4 bmax = c+e;
	float4 tmin = float4::min(float4::min(v0,v1),v2);
	float4 tmax = float4::max(float4::max(v0,v1),v2);

	float4 dj = float4::Or(float4::Greater(tmin,bmax),float4::Greater(bmin,tmax));
	//if(float4::AnyTrue(float4::EqualR(float4::swizzle(dj,0,1,2,0),float4::trueI())))
	//if(float4::AnyTrue(float4::EqualR(dj.swizzle<0,1,2,0>(),sint1::trueI())))
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

	//AnyFalse (EqualR(ni,true))
	//return float4::AnyFalse(float4::EqualR(ni,sint1::trueI()));
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
}

OctreeStructure::~OctreeStructure(){
	//
}

Octree::Octree(uint _x) : x(_x){//, lock(ATOMIC_FLAG_INIT){
	memset(pch,0,sizeof(pch));
}

Octree::Octree(){
	//memset(pch,0,sizeof(pch));
}

Octree::~Octree(){
	//
}

//Some template setup to combine BoundingBox/Triangle cases?
void Octree::BuildPath(const float4 &c, const float4 &e, const float4 &c1, const float4 &e1, uint level, uint mlevel, std::atomic<uint> *pindex, std::atomic<uint> *pleafx, tbb::concurrent_vector<Octree> *proot, tbb::concurrent_vector<OctreeStructure> *pob, VOLUME_BUFFER bx){
	m.lock();
	float4::store(&(*pob)[x].ce,float4::select(c,e,float4::selectctrl(0,0,0,1)));
	memset((*pob)[x].qval,0,sizeof((*pob)[x].qval)); //these are set during the resampling phase
	m.unlock();

	if(level >= mlevel-1){
		m.lock();//for(; lock.test_and_set(std::memory_order_acquire););
		if((*pob)[x].volx[bx] == ~0u)
			(*pob)[x].volx[bx] = pleafx->fetch_add(1);
		m.unlock();//lock.clear(std::memory_order_release);
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
			m.lock();//for(; lock.test_and_set(std::memory_order_acquire););
			if(!pch[i]){
				uint index = pindex->fetch_add(1)+1;
				pob->grow_to_at_least(index+1);

				//pch[i] = new(proot+index) Octree(index);
				pch[i] = new(scalable_malloc(sizeof(Octree))) Octree(index);//new(&(*proot)[index]) Octree(index);
				(*pob)[x].chn[i] = index;
			}
			m.unlock();//lock.clear(std::memory_order_release);

			pch[i]->BuildPath(cc,ee,c1,e1,level+1,mlevel,pindex,pleafx,proot,pob,bx);
		}
	}
}

void Octree::BuildPath(const float4 &c, const float4 &e, const float4 &v0, const float4 &v1, const float4 &v2, uint level, uint mlevel, std::atomic<uint> *pindex, std::atomic<uint> *pleafx, tbb::concurrent_vector<Octree> *proot, tbb::concurrent_vector<OctreeStructure> *pob, VOLUME_BUFFER bx){
	m.lock();
	float4::store(&(*pob)[x].ce,float4::select(c,e,float4::selectctrl(0,0,0,1)));
	memset((*pob)[x].qval,0,sizeof((*pob)[x].qval)); //these are set during the resampling phase
	m.unlock();

	if(level >= mlevel-1){
		m.lock();//for(; lock.test_and_set(std::memory_order_acquire););
		if((*pob)[x].volx[bx] == ~0u)
			(*pob)[x].volx[bx] = pleafx->fetch_add(1);
		m.unlock();//lock.clear(std::memory_order_release);
		return;
	}//else (*pob)[x].volx = ~0;

	float4 ee = 0.5f*e;
	for(uint i = 0; i < 8; ++i){
		float4 sv = 2.0f*float4((float)(i%2),(float)((i/2)%2),(float)(i/4),0.0f)-float4::one();
		float4 cc = c+sv*ee;
		BoundingBox aabb(cc,ee);

		if(aabb.Intersects(v0,v1,v2)){
			m.lock();//for(; lock.test_and_set(std::memory_order_acquire););
			if(!pch[i]){
				uint index = pindex->fetch_add(1)+1;
				pob->grow_to_at_least(index+1);

				pch[i] = new(scalable_malloc(sizeof(Octree))) Octree(index);//pch[i] = new(&(*proot)[index]) Octree(index);
				(*pob)[x].chn[i] = index;
			}
			m.unlock();//lock.clear(std::memory_order_release);

			pch[i]->BuildPath(cc,ee,v0,v1,v2,level+1,mlevel,pindex,pleafx,proot,pob,bx);
		}
	}
}

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
		if(psci->omask & 1<<Node::SceneInfo::OUTPUT_FLOAT_DISTANCE)
			return true;
	}
	return false;
}

using PostFogParams = std::tuple<SceneData::BaseObject *, openvdb::FloatGrid::Ptr>;
enum PFP{
	PFP_OBJECT,
	PFP_INPUTGRID
};

static void S_Create(float s, float qb, float lvc, float bvc, uint maxd, openvdb::FloatGrid::Ptr pgrid[VOLUME_BUFFER_COUNT], Scene *pscene){
	openvdb::math::Transform::Ptr pgridtr = openvdb::math::Transform::createLinearTransform(s);
	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i)
		pgrid[i] = 0;

	//Find if SceneInfo distance output was used anywhere in the node trees and automatically determine if a query field should be constructed.
	bool qfield =
		std::find_if(SceneData::Surface::objs.begin(),SceneData::Surface::objs.end(),S_FindSceneInfo) != SceneData::Surface::objs.end() ||
		std::find_if(SceneData::ParticleSystem::prss.begin(),SceneData::ParticleSystem::prss.end(),S_FindSceneInfo) != SceneData::ParticleSystem::prss.end() ||
		std::find_if(SceneData::SmokeCache::objs.begin(),SceneData::SmokeCache::objs.end(),S_FindSceneInfo) != SceneData::SmokeCache::objs.end();
	if(qfield)
		DebugPrintf("SceneInfo.distance in use, will construct a query field.\n");

	//Low-res query field
	openvdb::math::Transform::Ptr pqsdftr = openvdb::math::Transform::createLinearTransform(qb/bvc);
	openvdb::FloatGrid::Ptr pqsdf = openvdb::FloatGrid::create(qb);
	pqsdf->setTransform(pqsdftr);
	pqsdf->setGridClass(openvdb::GRID_LEVEL_SET);

	float4 scaabbmin = float4(FLT_MAX);
	float4 scaabbmax = -scaabbmin;

	//store the resulting surfaces - node trees are shared among objects, so the data is lost after each evaluation
	std::vector<dfloat3> vl;
	std::vector<duint3> tl;
	std::vector<duint4> ql;
	uint vca = 0;

	std::vector<BoundingBox> fogbvs;
	std::vector<PostFogParams> fogppl; //input grids to be post-processed

	for(uint i = 0; i < SceneData::Surface::objs.size(); ++i){
		Node::InputNodeParams snp(SceneData::Surface::objs[i],pgridtr,0,0,0,0);
		SceneData::Surface::objs[i]->pnt->EvaluateNodes1(&snp,0,1<<Node::OutputNode::INPUT_SURFACE);

		Node::BaseSurfaceNode1 *pdsn = dynamic_cast<Node::BaseSurfaceNode1*>(SceneData::Surface::objs[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_SURFACE]);
		if(pdsn->vl.size() > 0){
			openvdb::FloatGrid::Ptr ptgrid = pdsn->ComputeLevelSet(pgridtr,bvc,bvc);

			if(qfield){
				openvdb::FloatGrid::Ptr phgrid = pdsn->ComputeLevelSet(pqsdftr,bvc,2.0f);
				openvdb::tools::csgUnion(*pqsdf,*phgrid);
			}

			for(uint j = 0; j < pdsn->vl.size(); ++j){
				float4 p = float4::load((dfloat3*)&pdsn->vl[j]);
				scaabbmin = float4::min(p,scaabbmin);
				scaabbmax = float4::max(p,scaabbmax);

				vl.push_back(*(dfloat3*)pdsn->vl[j].asPointer());
			}

			for(uint j = 0; j < pdsn->tl.size(); ++j){
				openvdb::Vec3I t = pdsn->tl[j]+openvdb::Vec3I(vca);
				tl.push_back(*(duint3*)t.asPointer());
			}
			for(uint j = 0; j < pdsn->ql.size(); ++j){
				openvdb::Vec4I q = pdsn->ql[j]+openvdb::Vec4I(vca);
				ql.push_back(*(duint4*)q.asPointer());
			}

			vca += vl.size();

			if(pgrid[VOLUME_BUFFER_SDF])
			   openvdb::tools::csgUnion(*pgrid[VOLUME_BUFFER_SDF],*ptgrid);
			else pgrid[VOLUME_BUFFER_SDF] = ptgrid;
		}
	}

	openvdb::FloatGrid::Ptr ptfog = openvdb::FloatGrid::create(); //temporary grid to store the fog to be post-processed
	ptfog->setTransform(pgridtr);
	ptfog->setGridClass(openvdb::GRID_FOG_VOLUME);

	openvdb::Vec3SGrid::Ptr ptvel = openvdb::Vec3SGrid::create();
	ptvel->setTransform(pgridtr);
	ptvel->setGridClass(openvdb::GRID_FOG_VOLUME);

	for(uint i = 0; i < SceneData::SmokeCache::objs.size(); ++i){
		Node::InputNodeParams snp(SceneData::SmokeCache::objs[i],pgridtr,0,0,0,0);
		SceneData::SmokeCache::objs[i]->pnt->EvaluateNodes1(&snp,0,1<<Node::OutputNode::INPUT_FOG);

		Node::BaseFogNode1 *pdfn = dynamic_cast<Node::BaseFogNode1*>(SceneData::SmokeCache::objs[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_FOG]);
		if(pdfn->pdgrid->activeVoxelCount() > 0){
			if(SceneData::SmokeCache::objs[i]->pnt->GetRoot()->imask & 1<<Node::OutputNode::INPUT_FOGPOST){
				fogppl.push_back(PostFogParams(SceneData::SmokeCache::objs[i],pdfn->pdgrid->deepCopy()));
				openvdb::tools::compMax(*ptfog,*pdfn->pdgrid);
			}else
			if(pgrid[VOLUME_BUFFER_FOG])
				openvdb::tools::compMax(*pgrid[VOLUME_BUFFER_FOG],*pdfn->pdgrid); //compSum
			else pgrid[VOLUME_BUFFER_FOG] = pdfn->pdgrid;
		}
	}

	for(uint i = 0; i < SceneData::ParticleSystem::prss.size(); ++i){
		Node::InputNodeParams snp(SceneData::ParticleSystem::prss[i],pgridtr,0,0,0,0);
		SceneData::ParticleSystem::prss[i]->pnt->EvaluateNodes1(&snp,0,1<<Node::OutputNode::INPUT_FOG|1<<Node::OutputNode::INPUT_VECTOR);

		Node::BaseFogNode1 *pdfn = dynamic_cast<Node::BaseFogNode1*>(SceneData::ParticleSystem::prss[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_FOG]);
		if(pdfn->pdgrid->activeVoxelCount() > 0){
			if(SceneData::ParticleSystem::prss[i]->pnt->GetRoot()->imask & 1<<Node::OutputNode::INPUT_FOGPOST){
				fogppl.push_back(PostFogParams(SceneData::ParticleSystem::prss[i],pdfn->pdgrid->deepCopy()));
				openvdb::tools::compMax(*ptfog,*pdfn->pdgrid);
			}else
			if(pgrid[VOLUME_BUFFER_FOG])
				openvdb::tools::compMax(*pgrid[VOLUME_BUFFER_FOG],*pdfn->pdgrid);
			else pgrid[VOLUME_BUFFER_FOG] = pdfn->pdgrid;
		}

		Node::BaseVectorFieldNode1 *pvfn = dynamic_cast<Node::BaseVectorFieldNode1*>(SceneData::ParticleSystem::prss[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_VECTOR]);
		if(pvfn->pvgrid->activeVoxelCount() > 0)
			openvdb::tools::compSum(*ptvel,*pvfn->pvgrid);
	}

	//fog post-processor
	if(fogppl.size() > 0){
		openvdb::FloatGrid::Ptr pqfog;
		if(pgrid[VOLUME_BUFFER_FOG]){
			pqfog = pgrid[VOLUME_BUFFER_FOG]->deepCopy();
			openvdb::tools::compMax(*pqfog,*ptfog);
		}else pqfog = ptfog;

		FloatGridBoxSampler *pqsampler = new FloatGridBoxSampler(*pqsdf);
		FloatGridBoxSampler *ppsampler = new FloatGridBoxSampler(*pqfog);
		VectorGridBoxSampler *pvsampler = new VectorGridBoxSampler(*ptvel);
		FloatGridBoxSampler *pdsampler = pgrid[VOLUME_BUFFER_SDF]?new FloatGridBoxSampler(*pgrid[VOLUME_BUFFER_SDF]):pqsampler;

		for(uint i = 0; i < fogppl.size(); ++i){
			SceneData::PostFog fobj(std::get<PFP_OBJECT>(fogppl[i])->pnt,std::get<PFP_INPUTGRID>(fogppl[i]));
			Node::InputNodeParams snp(&fobj,pgridtr,pdsampler,pqsampler,ppsampler,pvsampler);
			fobj.pnt->EvaluateNodes1(&snp,0,1<<Node::OutputNode::INPUT_FOGPOST);

			Node::BaseFogNode1 *pdfn = dynamic_cast<Node::BaseFogNode1*>(fobj.pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_FOGPOST]);
			if(pgrid[VOLUME_BUFFER_FOG])
				openvdb::tools::compMax(*pgrid[VOLUME_BUFFER_FOG],*pdfn->pdgrid);
			else pgrid[VOLUME_BUFFER_FOG] = pdfn->pdgrid;
		}

		delete pqsampler;
		delete ppsampler;
		delete pvsampler;
		if(pgrid[VOLUME_BUFFER_SDF])
			delete pdsampler;
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

	if(pgrid[VOLUME_BUFFER_FOG]){
		for(openvdb::FloatGrid::TreeType::LeafCIter m = pgrid[VOLUME_BUFFER_FOG]->tree().cbeginLeaf(); m; ++m){
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

			fogbvs.push_back(aabb);

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

	DebugPrintf("Center = (%.3f, %.3f, %.3f)\nExtents = (%.3f, %.3f, %.3f) (max w = %.3f)\n",
		scaabb.sc.x,scaabb.sc.y,scaabb.sc.z,scaabb.se.x,scaabb.se.y,scaabb.se.z,d);

	float y = powf(2.0f,(float)mlevel);
	float r = y*lvc; //sparse voxel resolution
	//The actual voxel size is now v = d/(2^k*N), where k is the octree depth and N=lvc.
	DebugPrintf("> Constructing octree (depth = %u, leaf = %u*%f = %f, sparse res = %u^3)...\n",mlevel,(uint)lvc,d/r,lvc*d/r,(uint)r);

	Octree *proot = new(scalable_malloc(sizeof(Octree))) Octree(0);

	pscene->ob.reserve(500000);
	pscene->ob.emplace_back();

	std::atomic<uint> indexa(0);
	std::atomic<uint> leafxa(0); //sdf
	std::atomic<uint> leafxb(0); //fog

#if 0
	tbb::parallel_for(tbb::blocked_range<size_t>(0,tl.size()),[&](const tbb::blocked_range<size_t> &nr){
		for(uint i = nr.begin(); i < nr.end(); ++i){
			float4 v0 = float4::load(&vl[tl[i].x]);
			float4 v1 = float4::load(&vl[tl[i].y]);
			float4 v2 = float4::load(&vl[tl[i].z]);
			proot->BuildPath(c,a,v0,v1,v2,0,mlevel,&indexa,&leafxa,&pscene->root,&pscene->ob,VOLUME_BUFFER_SDF);
		}
	});

	tbb::parallel_for(tbb::blocked_range<size_t>(0,ql.size()),[&](const tbb::blocked_range<size_t> &nr){
		for(uint i = nr.begin(); i < nr.end(); ++i){
			float4 v0, v1, v2;

			v0 = float4::load(&vl[ql[i].x]);
			v1 = float4::load(&vl[ql[i].y]);
			v2 = float4::load(&vl[ql[i].z]);
			proot->BuildPath(c,a,v0,v1,v2,0,mlevel,&indexa,&leafxa,&pscene->root,&pscene->ob,VOLUME_BUFFER_SDF);

			v0 = float4::load(&vl[ql[i].z]);
			v1 = float4::load(&vl[ql[i].w]);
			v2 = float4::load(&vl[ql[i].x]);
			proot->BuildPath(c,a,v0,v1,v2,0,mlevel,&indexa,&leafxa,&pscene->root,&pscene->ob,VOLUME_BUFFER_SDF);
		}
	});

	tbb::parallel_for(tbb::blocked_range<size_t>(0,fogbvs.size()),[&](const tbb::blocked_range<size_t> &nr){
		for(uint i = nr.begin(); i < nr.end(); ++i){
			float4 c1 = float4::load(&fogbvs[i].sc);
			float4 e1 = float4::load(&fogbvs[i].se);
			proot->BuildPath(c,a,c1,e1,0,mlevel,&indexa,&leafxb,&pscene->root,&pscene->ob,VOLUME_BUFFER_FOG);
		}
	});
#else
	for(uint i = 0; i < tl.size(); ++i){
		float4 v0 = float4::load(&vl[tl[i].x]);
		float4 v1 = float4::load(&vl[tl[i].y]);
		float4 v2 = float4::load(&vl[tl[i].z]);
		proot->BuildPath(c,a,v0,v1,v2,0,mlevel,&indexa,&leafxa,&pscene->root,&pscene->ob,VOLUME_BUFFER_SDF);
	}

	for(uint i = 0; i < ql.size(); ++i){
		float4 v0, v1, v2;

		v0 = float4::load(&vl[ql[i].x]);
		v1 = float4::load(&vl[ql[i].y]);
		v2 = float4::load(&vl[ql[i].z]);
		proot->BuildPath(c,a,v0,v1,v2,0,mlevel,&indexa,&leafxa,&pscene->root,&pscene->ob,VOLUME_BUFFER_SDF);

		v0 = float4::load(&vl[ql[i].z]);
		v1 = float4::load(&vl[ql[i].w]);
		v2 = float4::load(&vl[ql[i].x]);
		proot->BuildPath(c,a,v0,v1,v2,0,mlevel,&indexa,&leafxa,&pscene->root,&pscene->ob,VOLUME_BUFFER_SDF);
	}

	for(uint i = 0; i < fogbvs.size(); ++i){
		float4 c1 = float4::load(&fogbvs[i].sc);
		float4 e1 = float4::load(&fogbvs[i].se);
		proot->BuildPath(c,a,c1,e1,0,mlevel,&indexa,&leafxb,&pscene->root,&pscene->ob,VOLUME_BUFFER_FOG);
	}
#endif

	proot->FreeRecursive();
	scalable_free(proot);

	pscene->lvoxc = (uint)lvc;
	pscene->index = indexa;
	pscene->leafx[VOLUME_BUFFER_SDF] = leafxa;
	pscene->leafx[VOLUME_BUFFER_FOG] = leafxb;
}

namespace SceneData{

BaseObject::BaseObject(Node::NodeTree *_pnt) : pnt(_pnt){
	//
}

BaseObject::~BaseObject(){
	//
}

ParticleSystem::ParticleSystem(Node::NodeTree *_pnt) : BaseObject(_pnt){
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

SmokeCache::SmokeCache(Node::NodeTree *_pnt) : BaseObject(_pnt){
	SmokeCache::objs.push_back(this);
}

SmokeCache::~SmokeCache(){
	//
}

void SmokeCache::DeleteAll(){
	for(uint i = 0; i < objs.size(); ++i)
		delete objs[i];
	objs.clear();
}

std::vector<SmokeCache *> SmokeCache::objs;

Surface::Surface(Node::NodeTree *_pnt) : BaseObject(_pnt){
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

Scene::Scene(){
	//
}

Scene::~Scene(){
	//
}

void Scene::Initialize(float s, uint maxd, float qb, SCENE_CACHE_MODE cm){
	openvdb::initialize();

	const float lvc = 8.0f; //minimum number of voxels in an octree leaf
	const float bvc = 4.0f; //number of narrow band voxels counting from the surface

	openvdb::FloatGrid::Ptr pgrid[VOLUME_BUFFER_COUNT];// = {0};
	FloatGridBoxSampler *psampler[VOLUME_BUFFER_COUNT];

	openvdb::io::File vdbc("/tmp/droplet-fileid.vdb");
	try{
		if(cm != SCENE_CACHE_READ)
			throw(0);

		vdbc.open(false);
		pgrid[VOLUME_BUFFER_SDF] = openvdb::gridPtrCast<openvdb::FloatGrid>(vdbc.readGrid("surface-levelset"));
		pgrid[VOLUME_BUFFER_FOG] = 0;
		vdbc.close();

		{
			FILE *pf = fopen("/tmp/droplet-fileid.bin","rb");
			if(!pf)
				throw(0);
			fread(&lvoxc,1,4,pf);
			fread(&index,1,4,pf);
			fread(&leafx[VOLUME_BUFFER_SDF],1,4,pf);
			leafx[VOLUME_BUFFER_FOG] = 0;

			uint pobl = index+1;
			ob.grow_to_at_least(pobl);
			for(uint i = 0; i < pobl; ++i)
				fread(&ob[i],1,sizeof(OctreeStructure),pf);

			fclose(pf);
		}

	}catch(...){
		if(cm == SCENE_CACHE_READ){
			DebugPrintf("Attempt to read VDB-cache failed. Writing to a new one.\n");
			cm = SCENE_CACHE_WRITE;
		}

		S_Create(s,qb,lvc,bvc,maxd,pgrid,this);
		if(pgrid[VOLUME_BUFFER_SDF])
			pgrid[VOLUME_BUFFER_SDF]->setName("surface-levelset");

		if(cm == SCENE_CACHE_WRITE){
			openvdb::GridCPtrVec gvec{pgrid[VOLUME_BUFFER_SDF]}; //include the fog grid also
			vdbc.write(gvec);
			vdbc.close();

			FILE *pf = fopen("/tmp/droplet-fileid.bin","wb");
			fwrite(&lvoxc,1,4,pf);
			fwrite(&index,1,4,pf);
			fwrite(&leafx[VOLUME_BUFFER_SDF],1,4,pf);
			//fwrite(pob,1,(index+1)*sizeof(OctreeStructure),pf);
			for(uint i = 0; i < index+1; ++i)
				fwrite(&ob[i],1,sizeof(OctreeStructure),pf);

			fclose(pf);
		}
	}

	DebugPrintf("> Resampling volume data...\n");
	
	lvoxc3 = lvoxc*lvoxc*lvoxc;
	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i){
		if(pgrid[i]){
			pvol[i] = new float[lvoxc3*leafx[i]];
			psampler[i] = new FloatGridBoxSampler(*pgrid[i]); //non-cached, thread safe version
		}else{
			pvol[i] = 0;
			psampler[i] = 0;
		}
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
				if(psampler[VOLUME_BUFFER_SDF]){
					float d = psampler[VOLUME_BUFFER_SDF]->wsSample(*(openvdb::Vec3f*)&ob[i].ce);
					if(d < 0.0f){
						//If the fog leaf is completely inside the sdf surface (no overlapping sdf leaf -> volx == ~0u),
						//remove it. It's useless there, and removing it simplifies the space skipping algorithm.
						ob[i].volx[VOLUME_BUFFER_FOG] = ~0u;
						continue;
					}
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
		if(psampler[i])
			delete psampler[i];

	uint sdfs = leafx[VOLUME_BUFFER_SDF]*lvoxc3*sizeof(float);
	uint fogs = leafx[VOLUME_BUFFER_FOG]*lvoxc3*sizeof(float);
	DebugPrintf("Volume size = %f MB\n  SDF = %f MB\n  Fog = %f MB\n",(float)(sdfs+fogs)/1e6f,(float)sdfs/1e6f,(float)fogs/1e6f);
	//
}

void Scene::Destroy(){
	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i)
		if(pvol[i])
			delete []pvol[i];
	ob.clear();
}
