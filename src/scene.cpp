#include "main.h"
#include "node.h"
#include "scene.h"
#include "noise.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h> //samplers
#include <openvdb/tools/Composite.h> //csg/comp

#include "SceneSurface.h"
#include "SceneDensity.h"

#include <tbb/parallel_for.h>

namespace Node{

using InputNodeParams = std::tuple<SceneData::BaseObject *, openvdb::math::Transform::Ptr>;
enum INP{
	INP_OBJECT,
	INP_TRANSFORM
};

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

//Some intersection code ripped (borrowed) from DirectXCollision library
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
    if(float4::AnyTrue(float4::EqualR(dj.swizzle<0,1,2,0>(),sint1::trueI())))
        return false;

    float4 n = float4::cross(v1-v0,v2-v0);
    float4 d = float4::dot3(n,v0);

    if(float4::AllTrue(float4::EqualR(n,z)))
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
    return float4::AnyFalse(float4::EqualR(ni,sint1::trueI()));
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
    //return !float4::AnyTrue(float4::EqualR(float4::swizzle(dj,0,1,2,0),float4::trueI()));
    return !float4::AnyTrue(float4::EqualR(dj.swizzle<0,1,2,0>(),sint1::trueI()));
}

Octree::Octree(uint _x) : x(_x){//, lock(ATOMIC_FLAG_INIT){
	memset(pch,0,sizeof(pch));
}

Octree::~Octree(){
    //
}

//Some template setup to combine BoundingBox/Triangle cases?
void Octree::BuildPath(const float4 &c, const float4 &e, const float4 &c1, const float4 &e1, uint level, uint mlevel, std::atomic<uint> *pindex, std::atomic<uint> *pleafx, Octree *proot, OctreeStructure *pob, VOLUME_BUFFER bx){
    float4::store(&pob[x].ce,float4::select(c,e,float4::selectctrl(0,0,0,1)));
    memset(pob[x].qval,0,sizeof(pob[x].qval)); //these are set during the resampling phase

    if(level >= mlevel-1){
        m.lock();//for(; lock.test_and_set(std::memory_order_acquire););
        if(pob[x].volx[bx] == ~0u)
            pob[x].volx[bx] = pleafx->fetch_add(1);
        m.unlock();//lock.clear(std::memory_order_release);
        return;
    }//else pob[x].volx = ~0;

    BoundingBox aabb1(c1,e1);

    float4 ee = 0.5f*e;
    for(uint i = 0; i < 8; ++i){
        float4 sv = 2.0f*float4((float)(i%2),(float)((i/2)%2),(float)(i/4),0.0f)-float4::one();
        float4 cc = c+sv*ee;
        BoundingBox aabb(cc,ee);

        //aabb.Intersects
        if(aabb.Intersects(aabb1)){
            m.lock();//for(; lock.test_and_set(std::memory_order_acquire););
            if(!pch[i]){
                uint index = pindex->fetch_add(1)+1;
                pch[i] = new(proot+index) Octree(index);
                pob[x].chn[i] = index;
            }
            m.unlock();//lock.clear(std::memory_order_release);

            pch[i]->BuildPath(cc,ee,c1,e1,level+1,mlevel,pindex,pleafx,proot,pob,bx);
        }
    }
}

void Octree::BuildPath(const float4 &c, const float4 &e, const float4 &v0, const float4 &v1, const float4 &v2, uint level, uint mlevel, std::atomic<uint> *pindex, std::atomic<uint> *pleafx, Octree *proot, OctreeStructure *pob, VOLUME_BUFFER bx){
    float4::store(&pob[x].ce,float4::select(c,e,float4::selectctrl(0,0,0,1)));
    memset(pob[x].qval,0,sizeof(pob[x].qval)); //these are set during the resampling phase

    if(level >= mlevel-1){
        m.lock();//for(; lock.test_and_set(std::memory_order_acquire););
        if(pob[x].volx[bx] == ~0u)
            pob[x].volx[bx] = pleafx->fetch_add(1);
        m.unlock();//lock.clear(std::memory_order_release);
        return;
    }//else pob[x].volx = ~0;

    float4 ee = 0.5f*e;
    for(uint i = 0; i < 8; ++i){
        float4 sv = 2.0f*float4((float)(i%2),(float)((i/2)%2),(float)(i/4),0.0f)-float4::one();
        float4 cc = c+sv*ee;
        BoundingBox aabb(cc,ee);

        if(aabb.Intersects(v0,v1,v2)){
            m.lock();//for(; lock.test_and_set(std::memory_order_acquire););
            if(!pch[i]){
                uint index = pindex->fetch_add(1)+1;
                pch[i] = new(proot+index) Octree(index);
                pob[x].chn[i] = index;
            }
            m.unlock();//lock.clear(std::memory_order_release);

            pch[i]->BuildPath(cc,ee,v0,v1,v2,level+1,mlevel,pindex,pleafx,proot,pob,bx);
        }
    }
}

/*class ParticleList{
public:
	typedef openvdb::Vec3R value_type;
	ParticleList(openvdb::Real rscale1 = 1.0f, openvdb::Real vscale1 = 1.0f) : rscale(rscale1), vscale(vscale1){}
	ParticleList(std::vector<openvdb::Vec4s> *psph1) : psph(psph1), rscale(1.0f), vscale(1.0f){}

	size_t size() const{
		return psph->size();
	}

	void getPos(size_t n, openvdb::Vec3R &p) const{
		p = (*psph)[n].getVec3();
	}

	void getPosRad(size_t n, openvdb::Vec3R &p, openvdb::Real &r) const{
		p = (*psph)[n].getVec3();
		r = (*psph)[n].w()*rscale;
	}

	void getAtt(size_t n, openvdb::Int32 &a){
		a = openvdb::Int32(n);
	}

protected:
	std::vector<openvdb::Vec4s> *psph;
	openvdb::Real rscale;
	openvdb::Real vscale;
};*/

static void S_Create(float s, float lff, openvdb::FloatGrid::Ptr pgrid[VOLUME_BUFFER_COUNT], Scene *pscene){
    openvdb::math::Transform::Ptr pgridtr = openvdb::math::Transform::createLinearTransform(s);
	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i)
		pgrid[i] = 0;

    float4 scaabbmin = float4(FLT_MAX);
    float4 scaabbmax = -scaabbmin;

    //store the resulting surfaces - node trees are shared among objects, so the data is lost after each evaluation
    std::vector<dfloat3> vl;
    std::vector<duint3> tl;
    std::vector<duint4> ql;
    uint vca = 0;

	std::vector<BoundingBox> fogbvs;

	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i)
		pgrid[i] = 0;

    for(uint i = 0; i < SceneData::Surface::objs.size(); ++i){
		Node::InputNodeParams snp(SceneData::Surface::objs[i],pgridtr);
		SceneData::Surface::objs[i]->pnt->EvaluateNodes1(&snp,0,1<<Node::OutputNode::INPUT_SURFACE);

        //dynamic cast to BaseSurfaceNode1 - getting empty
        Node::BaseSurfaceNode1 *pdsn = dynamic_cast<Node::BaseSurfaceNode1*>(SceneData::Surface::objs[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_SURFACE]);
        if(pdsn->vl.size() > 0){
            //openvdb::FloatGrid::Ptr ptgrid = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(*pgridtr,pdsn->vl,pdsn->tl,pdsn->ql,lff,lff);
			openvdb::FloatGrid::Ptr ptgrid = pdsn->ComputeLevelSet(pgridtr,lff);

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

        //TODO: get also the field sdf here, if available.
    }

	for(uint i = 0; i < SceneData::SmokeCache::objs.size(); ++i){
		Node::InputNodeParams snp(SceneData::SmokeCache::objs[i],pgridtr);
		SceneData::SmokeCache::objs[i]->pnt->EvaluateNodes1(&snp,0,1<<Node::OutputNode::INPUT_FOG);

		Node::BaseFogNode1 *pdfn = dynamic_cast<Node::BaseFogNode1*>(SceneData::SmokeCache::objs[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_FOG]);
		if(pdfn->pdgrid->activeVoxelCount() > 0){
			for(openvdb::FloatGrid::TreeType::LeafCIter m = pdfn->pdgrid->tree().cbeginLeaf(); m; ++m){
				const openvdb::FloatGrid::TreeType::LeafNodeType *pl = m.getLeaf();

				openvdb::math::CoordBBox bbox;// = pl->getNodeBoundingBox();
				pl->evalActiveBoundingBox(bbox);

				openvdb::math::Coord bdim = bbox.dim();
				openvdb::Vec3s dimi = bdim.asVec3s(); //index-space dimension
				openvdb::Vec3s extw = 0.5f*dimi*s; //pgrid1->transform().indexToWorld(dimi);

				openvdb::Vec3d posi = bbox.getCenter();
				openvdb::Vec3d posw = pdfn->pdgrid->transformPtr()->indexToWorld(posi);

				BoundingBox aabb;
				aabb.sc = dfloat3(posw.x(),posw.y(),posw.z());
				aabb.se = dfloat3(extw.x(),extw.y(),extw.z());

				fogbvs.push_back(aabb);

				float4 c = float4::load(&aabb.sc);
				float4 e = float4::load(&aabb.se);
				scaabbmin = float4::min(c-e,scaabbmin);
                scaabbmax = float4::max(c+e,scaabbmax);
			}

			if(pgrid[VOLUME_BUFFER_FOG])
				openvdb::tools::compSum(*pgrid[VOLUME_BUFFER_FOG],*pdfn->pdgrid);
			else pgrid[VOLUME_BUFFER_FOG] = pdfn->pdgrid;
		}
	}

	for(uint i = 0; i < SceneData::ParticleSystem::prss.size(); ++i){
		Node::InputNodeParams snp(SceneData::ParticleSystem::prss[i],pgridtr);
        SceneData::ParticleSystem::prss[i]->pnt->EvaluateNodes1(&snp,0,1<<Node::OutputNode::INPUT_FOG);

		Node::BaseFogNode1 *pdfn = dynamic_cast<Node::BaseFogNode1*>(SceneData::ParticleSystem::prss[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_FOG]);
		if(pdfn->pdgrid->activeVoxelCount() > 0){
			for(openvdb::FloatGrid::TreeType::LeafCIter m = pdfn->pdgrid->tree().cbeginLeaf(); m; ++m){
				const openvdb::FloatGrid::TreeType::LeafNodeType *pl = m.getLeaf();

				openvdb::math::CoordBBox bbox;// = pl->getNodeBoundingBox();
				pl->evalActiveBoundingBox(bbox);

				openvdb::math::Coord bdim = bbox.dim();
				openvdb::Vec3s dimi = bdim.asVec3s(); //index-space dimension
				openvdb::Vec3s extw = 0.5f*dimi*s; //pgrid1->transform().indexToWorld(dimi);

				openvdb::Vec3d posi = bbox.getCenter();
				openvdb::Vec3d posw = pdfn->pdgrid->transformPtr()->indexToWorld(posi);

				BoundingBox aabb;
				aabb.sc = dfloat3(posw.x(),posw.y(),posw.z());
				aabb.se = dfloat3(extw.x(),extw.y(),extw.z());

				fogbvs.push_back(aabb);

				float4 c = float4::load(&aabb.sc);
				float4 e = float4::load(&aabb.se);
				scaabbmin = float4::min(c-e,scaabbmin);
                scaabbmax = float4::max(c+e,scaabbmax);
			}

			if(pgrid[VOLUME_BUFFER_FOG])
				openvdb::tools::compSum(*pgrid[VOLUME_BUFFER_FOG],*pdfn->pdgrid);
			else pgrid[VOLUME_BUFFER_FOG] = pdfn->pdgrid;
			//pgrid[VOLUME_BUFFER_FOG] = pdfn->pdgrid;
		}
	}

    float4 c = 0.5f*(scaabbmax+scaabbmin);
    float4 e = 0.5f*(scaabbmax-scaabbmin);

    BoundingBox scaabb(c,e);
    float4 a = float4::max(float4::max(e.splat<0>(),e.splat<1>()),e.splat<2>());

    float d = 2.0f*a.get<0>();
    float N = (float)BLCLOUD_uN;
    uint mlevel = (uint)ceilf(logf(d/(N*s))/0.69315f);
    float r = powf(2.0f,(float)mlevel)*N; //sparse voxel resolution
    //The actual voxel size is now v = d/(2^k*N), where k is the octree depth.

    DebugPrintf("Center = (%f, %f, %f)\nExtents = (%f, %f, %f) (max w = %f)\n",scaabb.sc.x,scaabb.sc.y,scaabb.sc.z,scaabb.se.x,scaabb.se.y,scaabb.se.z,d);
    DebugPrintf("> Constructing octree (depth = %u, voxel = %f, sparse res = %u^3)...\n",mlevel,d/r,(uint)r);

#define BLCLOUD_MAXNODES 500000 //there should be some user-defined limit here or some way to estimate this
    uint octl = BLCLOUD_MAXNODES*sizeof(Octree);
    Octree *proot = (Octree*)_mm_malloc(octl,16);
    new(proot) Octree(0);

    uint pobl = BLCLOUD_MAXNODES*sizeof(OctreeStructure);
    pscene->pob = (OctreeStructure*)_mm_malloc(pobl,16);
    memset(pscene->pob,0,pobl);

    for(uint i = 0; i < pobl/sizeof(OctreeStructure); ++i)
		for(uint j = 0; j < VOLUME_BUFFER_COUNT; ++j)
        	pscene->pob[i].volx[j] = ~0u;

    std::atomic<uint> indexa(0);
    std::atomic<uint> leafxa(0); //sdf
	std::atomic<uint> leafxb(0); //fog
    tbb::parallel_for(tbb::blocked_range<size_t>(0,tl.size()),[&](const tbb::blocked_range<size_t> &nr){
        for(uint i = nr.begin(); i < nr.end(); ++i){
            float4 v0 = float4::load(&vl[tl[i].x]);
            float4 v1 = float4::load(&vl[tl[i].y]);
            float4 v2 = float4::load(&vl[tl[i].z]);
            proot->BuildPath(c,a,v0,v1,v2,0,mlevel,&indexa,&leafxa,proot,pscene->pob,VOLUME_BUFFER_SDF);
        }
    });

    tbb::parallel_for(tbb::blocked_range<size_t>(0,ql.size()),[&](const tbb::blocked_range<size_t> &nr){
        for(uint i = nr.begin(); i < nr.end(); ++i){
            float4 v0, v1, v2;

            v0 = float4::load(&vl[ql[i].x]);
            v1 = float4::load(&vl[ql[i].y]);
            v2 = float4::load(&vl[ql[i].z]);
            proot->BuildPath(c,a,v0,v1,v2,0,mlevel,&indexa,&leafxa,proot,pscene->pob,VOLUME_BUFFER_SDF);

            v0 = float4::load(&vl[ql[i].z]);
            v1 = float4::load(&vl[ql[i].w]);
            v2 = float4::load(&vl[ql[i].x]);
            proot->BuildPath(c,a,v0,v1,v2,0,mlevel,&indexa,&leafxa,proot,pscene->pob,VOLUME_BUFFER_SDF);
        }
    });

	tbb::parallel_for(tbb::blocked_range<size_t>(0,fogbvs.size()),[&](const tbb::blocked_range<size_t> &nr){
        for(uint i = nr.begin(); i < nr.end(); ++i){
			float4 c1 = float4::load(&fogbvs[i].sc);
			float4 e1 = float4::load(&fogbvs[i].se);
            proot->BuildPath(c,a,c1,e1,0,mlevel,&indexa,&leafxb,proot,pscene->pob,VOLUME_BUFFER_FOG);
        }
    });

    _mm_free(proot);

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

void Scene::Initialize(float s, SCENE_CACHE_MODE cm){
	openvdb::initialize();

    const float lff = 4.0f; //levelset offset (narrow band voxels counting from the surface)

    openvdb::FloatGrid::Ptr pgrid[VOLUME_BUFFER_COUNT];// = {0};
    //openvdb::io::File vdbc = openvdb::io::File("/tmp/droplet-fileid.vdb");
	openvdb::io::File vdbc("/tmp/droplet-fileid.vdb");
    try{
        if(cm != SCENE_CACHE_READ)
            throw(0);

        vdbc.open(false);
        pgrid[VOLUME_BUFFER_SDF] = openvdb::gridPtrCast<openvdb::FloatGrid>(vdbc.readGrid("surface-levelset"));
        vdbc.close();

        {
            FILE *pf = fopen("/tmp/droplet-fileid.bin","rb"); //TODO: use grid metadata to store this info
            if(!pf)
                throw(0);
            fread(&index,1,4,pf);
            fread(&leafx[VOLUME_BUFFER_SDF],1,4,pf);

            uint pobl = (index+1)*sizeof(OctreeStructure);
            pob = (OctreeStructure*)_mm_malloc(pobl,16);
            fread(pob,1,pobl,pf);

            fclose(pf);
        }

    }catch(...){
        if(cm == SCENE_CACHE_READ){
            DebugPrintf("Attempt to read VDB-cache failed. Writing to a new one.\n");
            cm = SCENE_CACHE_WRITE;
        }

        //S_Create(pvl,ptl,s,lff,&pgrid,&pob,&index,&leafx);
        S_Create(s,lff,pgrid,this);
		if(pgrid[VOLUME_BUFFER_SDF])
        	pgrid[VOLUME_BUFFER_SDF]->setName("surface-levelset");

        if(cm == SCENE_CACHE_WRITE){
            openvdb::GridCPtrVec gvec{pgrid[VOLUME_BUFFER_SDF]}; //include the fog grid also
            vdbc.write(gvec);
            vdbc.close();

            FILE *pf = fopen("/tmp/droplet-fileid.bin","wb");
            fwrite(&index,1,4,pf);
            fwrite(&leafx[VOLUME_BUFFER_SDF],1,4,pf);
            fwrite(pob,1,(index+1)*sizeof(OctreeStructure),pf);

            fclose(pf);
        }
    }

	DebugPrintf("> Resampling volume data...\n");

	//openvdb::FloatGrid::ConstAccessor
	//openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor, openvdb::tools::BoxSampler> fsampler(pgrid->getConstAccessor(),pgrid->transform());

	openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> *psampler[VOLUME_BUFFER_COUNT];
	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i){
		if(pgrid[i]){
			pbuf[i] = new LeafVolume[leafx[i]];
			psampler[i] = new openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler>(*pgrid[i]);
		}else{
			pbuf[i] = 0;
			psampler[i] = 0;
		}
	}
    //openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> samplerd(*pgrid[VOLUME_BUFFER_SDF]); //non-cached, thread safe version
	//openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> samplerf(*pgrid[VOLUME_BUFFER_FOG]);

    //float4 nv = float4(N);
    const uint uN = BLCLOUD_uN;
    float4 nv = float4((float)uN);

    /*typedef openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor, openvdb::tools::BoxSampler> FastGridSampler;
    tbb::enumerable_thread_specific<FastGridSampler> fsampler([&]()->FastGridSampler{
        return FastGridSampler(pgrid->getConstAccessor(),pgrid->transform()); //Seems to cause severe problems sometimes.
    });*/
	tbb::parallel_for(tbb::blocked_range<size_t>(0,index),[&](const tbb::blocked_range<size_t> &nr){
        //FastGridSampler &ffs = fsampler.local();

		openvdb::Vec3f posw;
		for(uint i = nr.begin(); i < nr.end(); ++i){
			if(pob[i].volx[VOLUME_BUFFER_SDF] == ~0u){
				if(pob[i].volx[VOLUME_BUFFER_FOG] == ~0u)
					continue; //not a leaf; exit early
				if(psampler[VOLUME_BUFFER_SDF]){
					float d = psampler[VOLUME_BUFFER_SDF]->wsSample(*(openvdb::Vec3f*)&pob[i].ce);
					if(d < 0.0f){
						//If the fog leaf is completely inside the sdf surface (no overlapping sdf leaf -> volx == ~0u),
						//remove it. It's useless there, and it causes incorrect space skipping.
						pob[i].volx[VOLUME_BUFFER_FOG] = ~0u;
						continue;
					}
				}
			}
			//
            float4 nc = float4::load(&pob[i].ce);
            float4 ne = nc.splat<3>();
            //
			pob[i].qval[VOLUME_BUFFER_SDF] = FLT_MAX;
			pob[i].qval[VOLUME_BUFFER_FOG] = 0.0f;
			//
            for(uint j = 0; j < uN*uN*uN; ++j){
                float4 nn = (nv-2.0f*float4((float)(j%uN),(float)((j/uN)%uN),(float)(j/(uN*uN)),1.0f)-float4::one())/(nv-float4::one());
                float4 nw = nc-ne*nn;
                float4::store((dfloat3*)posw.asPointer(),nw);

				if(pob[i].volx[VOLUME_BUFFER_SDF] != ~0u){
	                pbuf[VOLUME_BUFFER_SDF][pob[i].volx[VOLUME_BUFFER_SDF]].pvol[j] = psampler[VOLUME_BUFFER_SDF]->wsSample(posw);
					pob[i].qval[VOLUME_BUFFER_SDF] = openvdb::math::Min(pob[i].qval[VOLUME_BUFFER_SDF],pbuf[VOLUME_BUFFER_SDF][pob[i].volx[VOLUME_BUFFER_SDF]].pvol[j]);
				}

				if(pob[i].volx[VOLUME_BUFFER_FOG] != ~0u){
					pbuf[VOLUME_BUFFER_FOG][pob[i].volx[VOLUME_BUFFER_FOG]].pvol[j] = openvdb::math::Min(psampler[VOLUME_BUFFER_FOG]->wsSample(posw),0.07f);
					pob[i].qval[VOLUME_BUFFER_FOG] = openvdb::math::Max(pob[i].qval[VOLUME_BUFFER_FOG],pbuf[VOLUME_BUFFER_FOG][pob[i].volx[VOLUME_BUFFER_FOG]].pvol[j]);
				}
			}
		}
	});

	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i)
		if(psampler[i])
			delete psampler[i];

	uint sdfs = leafx[VOLUME_BUFFER_SDF]*sizeof(LeafVolume);
	uint fogs = leafx[VOLUME_BUFFER_FOG]*sizeof(LeafVolume);
	DebugPrintf("Volume size = %f MB\n  SDF = %f MB\n  Fog = %f MB\n",(float)(sdfs+fogs)/1e6f,(float)sdfs/1e6f,(float)fogs/1e6f);
    //
}

void Scene::Destroy(){
#if 0
	_aligned_free(pvolume);
#endif
	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i)
		if(pbuf[i])
    		delete []pbuf[i];
    _mm_free(pob);
}
