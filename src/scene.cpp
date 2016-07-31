#include "main.h"
#include "node.h"
#include "scene.h"
#include "noise.h"

#include <openvdb/openvdb.h>
//#include <openvdb/tools/LevelSetSphere.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h> //sdf rebuilding
#include <openvdb/tools/Interpolation.h>
//#include <openvdb/tools/LevelSetRebuild.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/tools/GridTransformer.h>

//#include <OpenVDB/tools/VolumeToSpheres.h>
//#include <OpenVDB/tools/ParticlesToLevelSet.h>

//#include <OpenVDB/io/Stream.h>

#include <openvdb/math/Operators.h> //CPT_RANGE

#include <tbb/tbb.h>
#include <tbb/parallel_for.h>

namespace Node{

//Additional layer to keep openvdb-specific data unexposed (the compile times are ridiculous)
class BaseFogNode1 : public virtual BaseFogNode{
public:
    BaseFogNode1(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt){
        pdgrid = openvdb::FloatGrid::create();
        pdgrid->setGridClass(openvdb::GRID_FOG_VOLUME);
    }

    ~BaseFogNode1(){
        //
    }

    openvdb::FloatGrid::Ptr pdgrid;
	//openvdb::Vec3SGrid::Ptr pvgrid; //sadly, there's no 4d vector grid
};

BaseFogNode * BaseFogNode::Create(uint level, NodeTree *pnt){
    return new BaseFogNode1(level,pnt);
}

template<class T>
using InputNodeParams = std::tuple<T *, openvdb::math::Transform::Ptr>;
enum INP{
	INP_OBJECT,
	INP_TRANSFORM
};

class BaseSurfaceNode1 : public virtual BaseSurfaceNode{
public:
    BaseSurfaceNode1(uint _level, NodeTree *pnt) : BaseSurfaceNode(_level,pnt){
        //
        pdgrid = openvdb::FloatGrid::create();
        pdgrid->setGridClass(openvdb::GRID_FOG_VOLUME);
        //transform shouldn't matter here, only index coordinates are used
    }

    ~BaseSurfaceNode1(){
        //
    }

    openvdb::FloatGrid::Ptr pdgrid; //billowing grid
    std::vector<openvdb::Vec3s> vl;
    std::vector<openvdb::Vec3I> tl;
    std::vector<openvdb::Vec4I> ql;
};

BaseSurfaceNode * BaseSurfaceNode::Create(uint level, NodeTree *pnt){
    return new BaseSurfaceNode1(level,pnt);
}

class SurfaceInput : public BaseSurfaceNode1, public ISurfaceInput{
public:
    SurfaceInput(uint _level, NodeTree *pnt) : BaseSurfaceNode(_level,pnt), BaseSurfaceNode1(_level,pnt), ISurfaceInput(_level,pnt){
        //
        //DebugPrintf(">> SurfaceInput()\n");
    }

    ~SurfaceInput(){
        //
    }

    void Evaluate(const void *pp){
        //DebugPrintf("---Surface::Evaluate()\n");
        InputNodeParams<SceneObject> *pd = (InputNodeParams<SceneObject>*)pp;
        SceneObject *pobj = std::get<INP_OBJECT>(*pd);
        //const std::vector<dfloat3> *pvl = std::get<SNP_VERTICES>(*pd);
        //const std::vector<uint> *ptl = std::get<SNP_TRIANGLES>(*pd);

        //vl.assign((openvdb::Vec3s*)pobj->vl.data(),(openvdb::Vec3s*)pobj->vl.data()+pobj->vl.size());
        //tl.assign((openvdb::Vec3I*)pobj->tl.data(),(openvdb::Vec3I*)pobj->tl.data()+pobj->tl.size()/3);

        //clear possible data from previous evaluation
        vl.clear();
        tl.clear();
        ql.clear();

        vl.reserve(pobj->vl.size());
        for(uint i = 0; i < pobj->vl.size(); ++i){
            openvdb::Vec3s t(pobj->vl[i].x,pobj->vl[i].y,pobj->vl[i].z);
            vl.push_back(t);
        }
        tl.reserve(pobj->tl.size());
        for(uint i = 0; i < pobj->tl.size(); i += 3){
            openvdb::Vec3I t(pobj->tl[i+0],pobj->tl[i+1],pobj->tl[i+2]);
            tl.push_back(t);
        }
    }
};

Node::ISurfaceInput * ISurfaceInput::Create(uint level, NodeTree *pnt){
    return new SurfaceInput(level,pnt);
}

class ParticleInput : public BaseFogNode1, public IParticleInput{
public:
	ParticleInput(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseFogNode1(_level,pnt), IParticleInput(_level,pnt){
        //
        DebugPrintf(">> ParticleInput()\n");
    }

    ~ParticleInput(){
        //
    }

    void Evaluate(const void *pp){
		InputNodeParams<ParticleSystem> *pd = (InputNodeParams<ParticleSystem>*)pp;
		ParticleSystem *pps = std::get<INP_OBJECT>(*pd);

		openvdb::math::Transform::Ptr pgridtr = std::get<INP_TRANSFORM>(*pd);
		pdgrid->setTransform(pgridtr);

		openvdb::math::Transform::Ptr pgridtr1 = openvdb::math::Transform::createLinearTransform(0.1f); //TODO: global setting (similar to detail size)
		openvdb::FloatGrid::Ptr ptgrid = openvdb::FloatGrid::create();
        ptgrid->setGridClass(openvdb::GRID_FOG_VOLUME);
		ptgrid->setTransform(pgridtr1);

		DebugPrintf("> Weighting particles...\n");

		openvdb::FloatGrid::Accessor grida = ptgrid->getAccessor();
		for(uint i = 0; i < pps->vl.size(); ++i){
			openvdb::Vec3d t(pps->vl[i].x,pps->vl[i].y,pps->vl[i].z);
			openvdb::Vec3d c = pgridtr1->worldToIndex(t); //assume cell-centered indices
			openvdb::Vec3d f = openvdb::Vec3d(floorf(c.x()-0.5f),floorf(c.y()-0.5f),floorf(c.z()-0.5f));
			openvdb::Vec3f b = c-f;

			//TODO: particle weight factor
			openvdb::Coord q((int)f.x(),(int)f.y(),(int)f.z());
			grida.modifyValue(q.offsetBy(0,0,0),[&](float &v){v += (1.0f-b.x())*(1.0f-b.y())*(1.0f-b.z());});
			grida.modifyValue(q.offsetBy(1,0,0),[&](float &v){v += b.x()*(1.0f-b.y())*(1.0f-b.z());});
			grida.modifyValue(q.offsetBy(0,1,0),[&](float &v){v += (1.0f-b.x())*b.y()*(1.0f-b.z());});
			grida.modifyValue(q.offsetBy(1,1,0),[&](float &v){v += b.x()*b.y()*(1.0f-b.z());});
			grida.modifyValue(q.offsetBy(0,0,1),[&](float &v){v += (1.0f-b.x())*(1.0f-b.y())*b.z();});
			grida.modifyValue(q.offsetBy(1,0,1),[&](float &v){v += b.x()*(1.0f-b.y())*b.z();});
			grida.modifyValue(q.offsetBy(0,1,1),[&](float &v){v += (1.0f-b.x())*b.y()*b.z();});
			grida.modifyValue(q.offsetBy(1,1,1),[&](float &v){v += b.x()*b.y()*b.z();});
		}

		DebugPrintf("> Upsampling particle fog...\n");

		//upsample the result
		openvdb::tools::resampleToMatch<openvdb::tools::BoxSampler>(*ptgrid,*pdgrid);
    }
};

Node::IParticleInput * IParticleInput::Create(uint level, NodeTree *pnt){
	return new ParticleInput(level,pnt);
}

class Displacement : public BaseSurfaceNode1, public IDisplacement{
public:
    Displacement(uint _level, NodeTree *pnt) : BaseSurfaceNode(_level,pnt), BaseSurfaceNode1(_level,pnt), IDisplacement(_level,pnt){
        //
        //DebugPrintf(">> Displacement()\n");
    }

    ~Displacement(){
        //
    }

    void Evaluate(const void *pp){
        //
        //DebugPrintf("---Displacement::Evaluate()\n");
    }
};

Node::IDisplacement * IDisplacement::Create(uint level, NodeTree *pnt){
    return new Displacement(level,pnt);
}

class fBmPerlinNoise : public BaseSurfaceNode1, public IfBmPerlinNoise{
public:
    fBmPerlinNoise(uint _level, NodeTree *pnt) : BaseSurfaceNode(_level,pnt), BaseSurfaceNode1(_level,pnt), IfBmPerlinNoise(_level,pnt){
        //
        //DebugPrintf(">> fBmPerlinNoise()\n");
    }

    ~fBmPerlinNoise(){
        //
    }

    void Evaluate(const void *pp){
        InputNodeParams<SceneObject> *pd = (InputNodeParams<SceneObject>*)pp;
        //DebugPrintf("---PerlinNoise::Evaluate()\n");

        const float lff = 4.0f;

        BaseValueNode<int> *poctn = dynamic_cast<BaseValueNode<int>*>(pnodes[IfBmPerlinNoise::INPUT_OCTAVES]);

        BaseValueNode<float> *pfreqn = dynamic_cast<BaseValueNode<float>*>(pnodes[IfBmPerlinNoise::INPUT_FREQ]);
        BaseValueNode<float> *pampn = dynamic_cast<BaseValueNode<float>*>(pnodes[IfBmPerlinNoise::INPUT_AMP]);
        BaseValueNode<float> *pfjumpn = dynamic_cast<BaseValueNode<float>*>(pnodes[IfBmPerlinNoise::INPUT_FJUMP]);
        BaseValueNode<float> *pgainn = dynamic_cast<BaseValueNode<float>*>(pnodes[IfBmPerlinNoise::INPUT_GAIN]);
        BaseValueNode<float> *pbilln = dynamic_cast<BaseValueNode<float>*>(pnodes[IfBmPerlinNoise::INPUT_BILLOW]);

        BaseSurfaceNode1 *pnode = dynamic_cast<BaseSurfaceNode1*>(pnodes[IfBmPerlinNoise::INPUT_SURFACE]);
        //BaseSurfaceNode1 *pgridn = dynamic_cast<BaseSurfaceNode1*>(pnodes[IfBmPerlinNoise::INPUT_GRID]); //could be just using input_surface

        //Temp: get max displacement distance. This won't work with non-constant node input values.
        //BaseValueNode<float>::EvaluateAll(0,level+1);
        //EvaluateValueGroup(level+1);
        //NodeTree::nts[0]->EvaluateLNodes(0,level+1);
        pntree->EvaluateNodes0(0,level+1,emask);
        float amp = fBm::GetAmplitudeMax(poctn->result,pampn->result,pgainn->result);

        openvdb::math::Transform::Ptr pgridtr = std::get<INP_TRANSFORM>(*pd);//openvdb::math::Transform::createLinearTransform(s);
        openvdb::FloatGrid::Ptr psgrid = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(*pgridtr,pnode->vl,pnode->tl,pnode->ql,ceilf(amp/pgridtr->voxelSize().x()+lff),lff);
        //openvdb::FloatGrid::Ptr pgrid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*pgridtr,points,tris,4.0f);
        //openvdb::FloatGrid::Ptr pgrid = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*pgridtr,points,tris,quads,4.0f);
        //DebugPrintf("Volume conversion face count = %u+0 (initial)\n",(uint)tris.size());

        DebugPrintf("Allocated disp. exterior narrow band for amp = %f+%f (%u voxels)\n",amp,pgridtr->voxelSize().x()*lff,(uint)ceilf(amp/pgridtr->voxelSize().x()+lff));

        DebugPrintf("> Displacing SDF...\n");

        /*openvdb::FloatGrid::Ptr pdgrid = openvdb::FloatGrid::create();
        pdgrid->setTransform(pgridtr);
        pdgrid->setGridClass(openvdb::GRID_FOG_VOLUME);*/

        /*openvdb::FloatGrid::Accessor grida = psgrid->getAccessor();
        openvdb::FloatGrid::Accessor gridd = pdgrid->getAccessor();
        openvdb::FloatGrid::ConstAccessor gridac = psgrid->getConstAccessor();
        openvdb::FloatGrid::ConstAccessor griddc = pdgrid->getConstAccessor();

        openvdb::math::CPT_RANGE<openvdb::math::UniformScaleMap,openvdb::math::CD_2ND> cptr;
        for(openvdb::FloatGrid::ValueOnIter m = psgrid->beginValueOn(); m.test(); ++m){
            //BaseValueNode<float>::EvaluateAll(0,level+1);
            EvaluateValueGroup(level+1);

            openvdb::Coord c = m.getCoord();
            openvdb::math::Vec3s posw = cptr.result(*pgridtr->map<openvdb::math::UniformScaleMap>(),gridac,c);
            sfloat4 sposw = sfloat1(posw.x(),posw.y(),posw.z(),0.0f); //TODO: utilize vectorization

            float f = fabsf(fBm::noise(sposw,poctn->result,pfreqn->result,pampn->result,pfjumpn->result,pgainn->result).get<0>());
            //if(i > 0)
                //f *= powf(std::min(gridd.getValue(c)/(*pdls)[i].qscale,1.0f),(*pdls)[i].billow);
            gridd.setValue(c,f);
        }

        for(openvdb::FloatGrid::ValueOnIter m = psgrid->beginValueOn(); m.test(); ++m){
            openvdb::Coord c = m.getCoord();
            float d = m.getValue();
            float f = griddc.getValue(c);

            //grida.setValue(c,d-f);
            m.setValue(d-f);
        }*/

        typedef std::tuple<openvdb::FloatGrid::Ptr, openvdb::FloatGrid::Accessor, openvdb::FloatGrid::ConstAccessor> FloatGridT;
        tbb::enumerable_thread_specific<FloatGridT> tgrida([&]()->FloatGridT{
            openvdb::FloatGrid::Ptr ptgrid = openvdb::FloatGrid::create();
            ptgrid->setTransform(pgridtr);
            ptgrid->setGridClass(openvdb::GRID_FOG_VOLUME);
            return FloatGridT(ptgrid,ptgrid->getAccessor(),psgrid->getConstAccessor());
        });
        openvdb::math::CPT_RANGE<openvdb::math::UniformScaleMap,openvdb::math::CD_2ND> cptr;
        tbb::parallel_for(openvdb::tree::IteratorRange<openvdb::FloatGrid::ValueOnIter>(psgrid->beginValueOn()),
            [&](openvdb::tree::IteratorRange<openvdb::FloatGrid::ValueOnIter> &r){
            FloatGridT &fgt = tgrida.local();
            for(; r; ++r){
                //TODO: need a thread-safe node evaluation
                //EvaluateValueGroup(level+1);
                const openvdb::FloatGrid::ValueOnIter &m = r.iterator();

                openvdb::Coord c = m.getCoord(); //, m.getValue()
                openvdb::math::Vec3s posw = cptr.result(*pgridtr->map<openvdb::math::UniformScaleMap>(),std::get<2>(fgt),c);
                sfloat4 sposw = sfloat1(posw.x(),posw.y(),posw.z(),0.0f); //TODO: utilize vectorization

                float f = fabsf(fBm::noise(sposw,poctn->result,pfreqn->result,pampn->result,pfjumpn->result,pgainn->result).get<0>());
                //f *= powf(std::min(dgridap.getValue(c),1.0f),pbilln->result);
                //TODO: normalize the previous layer displacement (f/amp when writing) so that qscale isn't needed at all.
                //if(i > 0)
                    //f *= powf(std::min(aa.second.getValue(c)/(*pdls)[i].qscale,1.0f),(*pdls)[i].billow);
                std::get<1>(fgt).setValue(c,f);
            }
        });

        //This could probably be faster with parallel reduction
        //for(tbb::enumerable_thread_specific<FloatGridT>::const_iterator m = tgrida.begin(); m != tgrida.end(); ++m)
        //openvdb::tools::compSum(*psgrid,*std::get<0>(*m));
        openvdb::FloatGrid::Accessor sgrida = psgrid->getAccessor();
        openvdb::FloatGrid::Accessor dgrida = pdgrid->getAccessor();
        openvdb::FloatGrid::ConstAccessor dgrida0 = pnode->pdgrid->getConstAccessor();
        for(tbb::enumerable_thread_specific<FloatGridT>::const_iterator q = tgrida.begin(); q != tgrida.end(); ++q){
            //
            for(openvdb::FloatGrid::ValueOnIter m = std::get<0>(*q)->beginValueOn(); m.test(); ++m){
                openvdb::Coord c = m.getCoord();

                float d = sgrida.getValue(c);
                float f = m.getValue();

                dgrida.setValue(c,f/amp);
                f *= powf(std::min(dgrida0.getValue(c),1.0f),pbilln->result); //this is here because of the normalization

                sgrida.setValue(c,d-f);

                //if(omask & (1<<OUTPUT_GRID))
                    //dgrida.setValue(c,d/amp);
            }
        }

        DebugPrintf("> Rebuilding...\n");

        openvdb::tools::volumeToMesh<openvdb::FloatGrid>(*psgrid,vl,tl,ql,0.0);
    }
};

IfBmPerlinNoise * IfBmPerlinNoise::Create(uint level, NodeTree *pnt){
    return new fBmPerlinNoise(level,pnt);
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
    if(float4::AnyTrue(float4::EqualR(dj.swizzle<0,1,2,0>(),float4::trueI())))
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
    return float4::AnyFalse(float4::EqualR(ni,float4::trueI()));
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
    return !float4::AnyTrue(float4::EqualR(dj.swizzle<0,1,2,0>(),float4::trueI()));
}

Octree::Octree(uint _x) : x(_x), lock(ATOMIC_FLAG_INIT){
	memset(pch,0,sizeof(pch));
}

Octree::~Octree(){
    //
}

//Some template setup to combine BoundingBox/Triangle cases?
void Octree::BuildPath(const float4 &c, const float4 &e, const float4 &c1, const float4 &e1, uint level, uint mlevel, std::atomic<uint> *pindex, std::atomic<uint> *pleafx, Octree *proot, OctreeStructure *pob, VOLUME_BUFFER bx){
    float4::store(&pob[x].ce,float4::select(c,e,float4::selectctrl(0,0,0,1)));
    pob[x].pmax = 0.0f;

    if(level >= mlevel-1){
        for(; lock.test_and_set(std::memory_order_acquire););
        if(pob[x].volx[bx] == ~0u)
            pob[x].volx[bx] = pleafx->fetch_add(1);//(*pleafx)++;
        lock.clear(std::memory_order_release);
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
            for(; lock.test_and_set(std::memory_order_acquire););
            if(!pch[i]){
                //pch[i] = NEWM(Octree)(++(*pindex));
                //pch[i] = new Octree(++(*pindex));
                //++(*pindex);
                uint index = pindex->fetch_add(1)+1;
                pch[i] = new(proot+index) Octree(index);
                pob[x].chn[i] = index;
            }
            lock.clear(std::memory_order_release);

            pch[i]->BuildPath(cc,ee,c1,e1,level+1,mlevel,pindex,pleafx,proot,pob,bx);
        }
    }
}

void Octree::BuildPath(const float4 &c, const float4 &e, const float4 &v0, const float4 &v1, const float4 &v2, uint level, uint mlevel, std::atomic<uint> *pindex, std::atomic<uint> *pleafx, Octree *proot, OctreeStructure *pob, VOLUME_BUFFER bx){
    float4::store(&pob[x].ce,float4::select(c,e,float4::selectctrl(0,0,0,1)));
    pob[x].pmax = 0.0f;

    if(level >= mlevel-1){
        for(; lock.test_and_set(std::memory_order_acquire););
        if(pob[x].volx[bx] == ~0u)
            pob[x].volx[bx] = pleafx->fetch_add(1);//(*pleafx)++;
        lock.clear(std::memory_order_release);
        return;
    }//else pob[x].volx = ~0;

    float4 ee = 0.5f*e;
    for(uint i = 0; i < 8; ++i){
        float4 sv = 2.0f*float4((float)(i%2),(float)((i/2)%2),(float)(i/4),0.0f)-float4::one();
        float4 cc = c+sv*ee;
        BoundingBox aabb(cc,ee);

        if(aabb.Intersects(v0,v1,v2)){
            for(; lock.test_and_set(std::memory_order_acquire););
            if(!pch[i]){
                //pch[i] = NEWM(Octree)(++(*pindex));
                //pch[i] = new Octree(++(*pindex));
                //++(*pindex);
                uint index = pindex->fetch_add(1)+1;
                pch[i] = new(proot+index) Octree(index);
                pob[x].chn[i] = index;
            }
            lock.clear(std::memory_order_release);

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

static void S_Create(float s, float lff, openvdb::FloatGrid::Ptr *pgrid, OctreeStructure **pob, uint *pindex, uint *pleafx){
    openvdb::math::Transform::Ptr pgridtr = openvdb::math::Transform::createLinearTransform(s);
    *pgrid = 0;

    float4 scaabbmin = float4(FLT_MAX);
    float4 scaabbmax = -scaabbmin;

    //store the resulting surfaces - node trees are shared among objects, so the data is lost after each evaluation
    std::vector<dfloat3> vl;
    std::vector<duint3> tl;
    std::vector<duint4> ql;
    uint vca = 0;

	std::vector<BoundingBox> fogbvs;

    for(uint i = 0; i < SceneObject::objs.size(); ++i){
		Node::InputNodeParams<SceneObject> snp(SceneObject::objs[i],pgridtr);
		SceneObject::objs[i]->pnt->EvaluateNodes1(&snp,0,1<<Node::OutputNode::INPUT_SURFACE);

        //dynamic cast to BaseSurfaceNode1 - getting empty
        Node::BaseSurfaceNode1 *pdsn = dynamic_cast<Node::BaseSurfaceNode1*>(SceneObject::objs[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_SURFACE]);
        if(pdsn->vl.size() > 0){
            openvdb::FloatGrid::Ptr ptgrid = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(*pgridtr,pdsn->vl,pdsn->tl,pdsn->ql,lff,lff);

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

            if(*pgrid)
               openvdb::tools::csgUnion(**pgrid,*ptgrid);
            else *pgrid = ptgrid;
        }

        /*Node::BaseSurfaceNode1 *pdfn = dynamic_cast<Node::BaseSurfaceNode1*>(SceneObject::objs[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_FIELD]);
        if(pdfn->vl.size() > 0){
            openvdb::FloatGrid::Ptr ptgrid = openvdb::tools::meshToSignedDistanceField<openvdb::FloatGrid>(*pgridtr,pdsn->vl,pdsn->tl,pdsn->ql,lff,lff);

            //openvdb::tools::sdfToFogVolume
            //combine the levelsets, then do ^^ at the end?
        }*/

        //TODO: get also the field sdf here, if available.
    }

	//Particle systems probably need their own material reference and node tree.
	for(uint i = 0; i < ParticleSystem::prss.size(); ++i){
		Node::InputNodeParams<ParticleSystem> snp(ParticleSystem::prss[i],pgridtr);
        ParticleSystem::prss[i]->pnt->EvaluateNodes1(&snp,0,1<<Node::OutputNode::INPUT_FOG);

		Node::BaseFogNode1 *pdfn = dynamic_cast<Node::BaseFogNode1*>(ParticleSystem::prss[i]->pnt->GetRoot()->pnodes[Node::OutputNode::INPUT_FOG]);
		{ //TODO: check if not empty
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
    *pob = (OctreeStructure*)_mm_malloc(pobl,16);
    memset(*pob,0,pobl);

    for(uint i = 0; i < pobl/sizeof(OctreeStructure); ++i)
		for(uint j = 0; j < VOLUME_BUFFER_COUNT; ++j)
        	(*pob)[i].volx[j] = ~0u;

    std::atomic<uint> indexa(0);
    std::atomic<uint> leafxa(0); //sdf
	std::atomic<uint> leafxb(0); //fog
    tbb::parallel_for(tbb::blocked_range<size_t>(0,tl.size()),[&](const tbb::blocked_range<size_t> &nr){
        for(uint i = nr.begin(); i < nr.end(); ++i){
            float4 v0 = float4::load(&vl[tl[i].x]);
            float4 v1 = float4::load(&vl[tl[i].y]);
            float4 v2 = float4::load(&vl[tl[i].z]);
            proot->BuildPath(c,a,v0,v1,v2,0,mlevel,&indexa,&leafxa,proot,*pob,VOLUME_BUFFER_SDF);
        }
    });

    tbb::parallel_for(tbb::blocked_range<size_t>(0,ql.size()),[&](const tbb::blocked_range<size_t> &nr){
        for(uint i = nr.begin(); i < nr.end(); ++i){
            float4 v0, v1, v2;

            v0 = float4::load(&vl[ql[i].x]);
            v1 = float4::load(&vl[ql[i].y]);
            v2 = float4::load(&vl[ql[i].z]);
            proot->BuildPath(c,a,v0,v1,v2,0,mlevel,&indexa,&leafxa,proot,*pob,VOLUME_BUFFER_SDF);

            v0 = float4::load(&vl[ql[i].z]);
            v1 = float4::load(&vl[ql[i].w]);
            v2 = float4::load(&vl[ql[i].x]);
            proot->BuildPath(c,a,v0,v1,v2,0,mlevel,&indexa,&leafxa,proot,*pob,VOLUME_BUFFER_SDF);
        }
    });

	tbb::parallel_for(tbb::blocked_range<size_t>(0,fogbvs.size()),[&](const tbb::blocked_range<size_t> &nr){
        for(uint i = nr.begin(); i < nr.end(); ++i){
			float4 c1 = float4::load(&fogbvs[i].sc);
			float4 e1 = float4::load(&fogbvs[i].se);
            proot->BuildPath(c,a,c1,e1,0,mlevel,&indexa,&leafxb,proot,*pob,VOLUME_BUFFER_FOG);
        }
    });

    _mm_free(proot);

    *pindex = indexa;
    *pleafx = leafxa;
}

ParticleSystem::ParticleSystem(Node::NodeTree *_pnt) : pnt(_pnt){
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

SceneObject::SceneObject(Node::NodeTree *_pnt) : pnt(_pnt){
    SceneObject::objs.push_back(this);
}

SceneObject::~SceneObject(){
    //
}

void SceneObject::DeleteAll(){
    for(uint i = 0; i < objs.size(); ++i)
        delete objs[i];
    objs.clear();
}

std::vector<SceneObject *> SceneObject::objs;

Scene::Scene(){
	//
}

Scene::~Scene(){
	//
}

void Scene::Initialize(float s, SCENE_CACHE_MODE cm){
	openvdb::initialize();

    const float lff = 4.0f; //levelset offset (narrow band voxels counting from the surface)

    openvdb::FloatGrid::Ptr pgrid = 0;
    openvdb::io::File vdbc = openvdb::io::File("/tmp/droplet-fileid.vdb");
    try{
        if(cm != SCENE_CACHE_READ)
            throw(0);

        vdbc.open(false);
        pgrid = openvdb::gridPtrCast<openvdb::FloatGrid>(vdbc.readGrid("surface-levelset"));
        vdbc.close();

        {
            FILE *pf = fopen("/tmp/droplet-fileid.bin","rb");
            if(!pf)
                throw(0);
            fread(&index,1,4,pf);
            fread(&leafx,1,4,pf);

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
        S_Create(s,lff,&pgrid,&pob,&index,&leafx);
        pgrid->setName("surface-levelset");

        if(cm == SCENE_CACHE_WRITE){
            openvdb::GridCPtrVec gvec{pgrid}; //include the fog grid also
            vdbc.write(gvec);
            vdbc.close();

            FILE *pf = fopen("/tmp/droplet-fileid.bin","wb");
            fwrite(&index,1,4,pf);
            fwrite(&leafx,1,4,pf);
            fwrite(pob,1,(index+1)*sizeof(OctreeStructure),pf);

            fclose(pf);
        }
    }
#if 0

//#define BLCLOUD_FOGVOLUME
#ifdef BLCLOUD_FOGVOLUME
    //Blender OpenVDB smoke cache usage testing. It works, but currently has no place here.
    openvdb::io::File vdbf = openvdb::io::File("/tmp/smoke_000130_00.vdb");
	vdbf.open();

	openvdb::FloatGrid::Ptr pgrid1 = openvdb::gridPtrCast<openvdb::FloatGrid>(vdbf.readGrid("density"));
	DebugPrintf("Grid: %s, FloatGrid: %u, Class: %u\n",pgrid1->getName().c_str(),pgrid1->isType<openvdb::FloatGrid>(),pgrid1->getGridClass());
	openvdb::Vec3d vs = pgrid1->voxelSize();
	DebugPrintf("Voxel size = %f, %f, %f\n",vs.x(),vs.y(),vs.z());

	openvdb::math::CoordBBox bbox = pgrid1->evalActiveVoxelBoundingBox();
	openvdb::math::Coord crd = bbox.dim();
	openvdb::Vec3d dimi = crd.asVec3d();
	openvdb::Vec3d extw = 0.5f*dimi*vs;//pgridtr->indexToWorld(dimi);

	openvdb::Vec3d posi = bbox.getCenter();
	openvdb::Vec3d posw = pgrid1->transformPtr()->indexToWorld(posi);//posi*vs;//pgridtr->indexToWorld(posi);
	//hopefully index(0,0,0) is defined to be the origin^^

	pgrid1->setGridClass(openvdb::GridClass::GRID_FOG_VOLUME);

	DebugPrintf("center = [%f, %f, %f], extents = [%f, %f, %f]\n",posw.x(),posw.y(),posw.z(),extw.x(),extw.y(),extw.z());

	vdbf.close();

	std::vector<BoundingBox> fogbvs; //TODO: preallocate by leaf count, if possible
	//openvdb::tools::fillWithSpheres(*pgrid1,spheres,10000,true,4.0f,16.0f,0.5f);
	for(openvdb::FloatGrid::TreeType::LeafCIter m = pgrid1->tree().cbeginLeaf(); m; ++m){
		const openvdb::FloatGrid::TreeType::LeafNodeType *pl = m.getLeaf();
		openvdb::math::CoordBBox bbox;// = pl->getNodeBoundingBox();
		pl->evalActiveBoundingBox(bbox);
		openvdb::math::Coord crd = bbox.dim();
		openvdb::Vec3d dimi = crd.asVec3d();
		openvdb::Vec3d extw = 0.5f*dimi*vs;//pgrid1->transform().indexToWorld(dimi);

		openvdb::Vec3d posi = bbox.getCenter();
		openvdb::Vec3d posw = pgrid1->transformPtr()->indexToWorld(posi);
		//openvdb::Vec4s s(posw.x(),posw.y(),posw.z(),sqrt(3.0f)*dimw.x()); //inclusive
		//openvdb::Vec4s s(posw.x(),posw.y(),posw.z(),dimw.x()); //exclusive

		BoundingBox aabb;
		aabb.Center = XMFLOAT3(posw.x(),posw.y(),posw.z());
		aabb.Extents = XMFLOAT3(extw.x(),extw.y(),extw.z());

		fogbvs.push_back(aabb);
	}

	//openvdb::FloatGrid::Ptr pgrida = openvdb::createLevelSet<openvdb::FloatGrid>(s,4.0f);
	/*openvdb::tools::ParticlesToLevelSet<openvdb::FloatGrid> lsf(*pgrid);
	ParticleList pl(&spheres);
	lsf.setGrainSize(1);
	lsf.rasterizeSpheres(pl);*/

	//DebugPrintf("Fog sphere count = %u\n",spheres.size());
	DebugPrintf("Fog subvolume count = %u\n",fogbvs.size());
#endif

    /*float4 scaabbmin = float4::load((dfloat3*)&pdn->vl[0]);//XMVECTOR scaabbmin = XMLoadFloat3((XMFLOAT3*)&points[0]);
    float4 scaabbmax = scaabbmin;//XMVECTOR scaabbmax = scaabbmin;
    uint pointc = pdn->vl.size();
    for(uint i = 0; i < pointc; ++i){
        float4 p = float4::load((dfloat3*)&pdn->vl[i]);//XMVECTOR p = XMLoadFloat3((XMFLOAT3*)&points[i]);
        scaabbmin = float4::min(p,scaabbmin);
        scaabbmax = float4::max(p,scaabbmax);
    }*/

#ifdef BLCLOUD_FOGVOLUME
	uint fogbvc = fogbvs.size();
	for(uint i = 0; i < fogbvc; ++i){
		XMVECTOR c = XMLoadFloat3(&fogbvs[i].Center);
		XMVECTOR e = XMLoadFloat3(&fogbvs[i].Extents);
		scaabbmin = XMVectorMin(c-e,scaabbmin);
		scaabbmax = XMVectorMax(c+e,scaabbmax);
	}
#endif

	/*uint spherec = spheres.size();
	XMFLOAT4 *psdata = (XMFLOAT4*)spheres.data();
	for(uint i = 0; i < spherec; ++i){
		XMVECTOR p = XMLoadFloat4(&psdata[i]);
		XMVECTOR r = XMVectorSplatW(p);
		scaabbmin = XMVectorMin(p-r,scaabbmin);
		scaabbmax = XMVectorMax(p+r,scaabbmax);
	}*/

#ifdef BLCLOUD_FOGVOLUME
	for(uint i = 0; i < fogbvc; ++i){
		XMVECTOR c1 = XMLoadFloat3(&fogbvs[i].Center);
		XMVECTOR e1 = XMLoadFloat3(&fogbvs[i].Extents);

		proot->BuildPath(c,a,c1,e1,0,mlevel,&index,&leafx,pob);
	}
#endif

#endif

	DebugPrintf("> Resampling volume data...\n");

    pbuf[VOLUME_BUFFER_SDF] = new LeafVolume[leafx];
    pbuf[VOLUME_BUFFER_FOG] = new LeafVolume[1];
	//openvdb::FloatGrid::ConstAccessor
	//openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor, openvdb::tools::BoxSampler> fsampler(pgrid->getConstAccessor(),pgrid->transform());

    openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> samplerd(*pgrid); //non-cached, thread safe version
#ifdef BLCLOUD_FOGVOLUME
	openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> samplerf(*pgrid1);
#endif

    //float4 nv = float4(N);
    const uint uN = BLCLOUD_uN;
    float4 nv = float4((float)uN);

    /*typedef openvdb::tools::GridSampler<openvdb::FloatGrid::ConstAccessor, openvdb::tools::BoxSampler> FastGridSampler;
    tbb::enumerable_thread_specific<FastGridSampler> fsampler([&]()->FastGridSampler{
        return FastGridSampler(pgrid->getConstAccessor(),pgrid->transform()); //Seems to cause severe problems sometimes.
    });*/
	tbb::parallel_for(tbb::blocked_range<size_t>(0,index),[&](const tbb::blocked_range<size_t> &nr){
        //FastGridSampler &ffs = fsampler.local();

		for(uint i = nr.begin(); i < nr.end(); ++i){
            if(pob[i].volx[VOLUME_BUFFER_SDF] == ~0u)
				continue;
            float4 nc = float4::load(&pob[i].ce);
            float4 ne = nc.splat<3>();
            //
            for(uint j = 0; j < uN*uN*uN; ++j){
                float4 nn = (nv-float4(2.0f)*float4((float)(j%uN),(float)((j/uN)%uN),(float)(j/(uN*uN)),1.0f)-float4::one())/(nv-float4::one());
                float4 nw = nc-ne*nn;

				openvdb::Vec3f posw;
                float4::store((dfloat3*)posw.asPointer(),nw);

/*#ifdef BLCLOUD_FOGVOLUME
				float d0 = samplerd.wsSample(posw);
				float p0 = samplerf.wsSample(posw);

				pvol[pob[i].volx].pvol[0][j] = d0 < s?d0:(p0 > 0.0f?s:d0);
				pvol[pob[i].volx].pvol[1][j] = (pvol[pob[i].volx].pvol[0][j] <= 0.0f)?1.0f:p0;
				pob[i].pmax = openvdb::math::Max(pob[i].pmax,pvol[pob[i].volx].pvol[1][j]);
#else
                pvol[pob[i].volx].pvol[0][j] = samplerd.wsSample(posw);//samplerd.wsSample(posw);
				pvol[pob[i].volx].pvol[1][j] = (pvol[pob[i].volx].pvol[0][j] <= 0.0f); //-1.0f
				pob[i].pmax = openvdb::math::Max(pob[i].pmax,pvol[pob[i].volx].pvol[1][j]);
#endif*/
                pbuf[VOLUME_BUFFER_SDF][pob[i].volx[VOLUME_BUFFER_SDF]].pvol[j] = samplerd.wsSample(posw);

				/*pvol[pob[i].volx].pvol[0][j] = samplerd.wsSample(posw);
				pvol[pob[i].volx].pvol[1][j] = (pvol[pob[i].volx].pvol[0][j] <= 0.0f); //-1.0f
				pob[i].pmax = openvdb::math::Max(pob[i].pmax,pvol[pob[i].volx].pvol[1][j]);*/

                //TODO: find max value and store this with other leaf params. Evaluate density at the given posw.
			}
		}
	});

	DebugPrintf("Volume size = %f MB\n",(float)(leafx*sizeof(LeafVolume))/1e6f);
    //
}

void Scene::Destroy(){
#if 0
	_aligned_free(pvolume);
#endif
	for(uint i = 0; i < VOLUME_BUFFER_COUNT; ++i)
    	delete []pbuf[i];
    _mm_free(pob);
}
