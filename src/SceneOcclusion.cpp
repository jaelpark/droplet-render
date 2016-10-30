#include "main.h"

//#ifdef USE_EMBREE
#include <embree2/rtcore.h>
#include <embree2/rtcore_ray.h>

//#endif

#include "scene.h"
#include "SceneOcclusion.h"

SceneOcclusion::SceneOcclusion(){
	//
}

SceneOcclusion::~SceneOcclusion(){
	//
}

void SceneOcclusion::Initialize(){
#ifdef USE_EMBREE
	DebugPrintf("> Preparing occlusion geometry...\n");

	pdev = rtcNewDevice(0);
	pscene = rtcDeviceNewScene(pdev,RTC_SCENE_STATIC|RTC_SCENE_INCOHERENT|RTC_SCENE_HIGH_QUALITY,RTC_INTERSECT4);
	//rtcIntersect4?

	for(uint i = 0; i < SceneData::Surface::objs.size(); ++i){
		if(!SceneData::Surface::objs[i]->holdout)
			continue;
		uint m = rtcNewTriangleMesh(pscene,RTC_GEOMETRY_STATIC,SceneData::Surface::objs[i]->tl.size()/3,SceneData::Surface::objs[i]->vl.size());

		dfloat4 *pvs = (dfloat4*)rtcMapBuffer(pscene,m,RTC_VERTEX_BUFFER);
		for(uint i = 0, n = SceneData::Surface::objs[i]->vl.size(); i < n; ++i){
			dfloat3 &v = SceneData::Surface::objs[i]->vl[i];
			pvs[i] = dfloat4(v.x,v.y,v.z,1.0f);
		}
		rtcUnmapBuffer(pscene,m,RTC_VERTEX_BUFFER);

		dint3 *pts = (dint3*)rtcMapBuffer(pscene,m,RTC_INDEX_BUFFER);
		for(uint i = 0, n = SceneData::Surface::objs[i]->tl.size()/3; i < n; ++i){
			for(uint j = 0; j < 3; ++j)
				((int*)&pts[i])[j] = SceneData::Surface::objs[i]->tl[3*i+j];
		}
		rtcUnmapBuffer(pscene,m,RTC_INDEX_BUFFER);
	}

	rtcCommit(pscene);
#endif
}

void SceneOcclusion::Intersect(const sfloat4 &ro, const sfloat4 &rd, const sfloat1 &gm, dintN *pmask, dfloatN *pdist) const{
	//TODO: check if number of occluders > 0 to skip the initialization when not needed
	//RTCRayNt<BLCLOUD_VSIZE> ray;
#ifdef USE_EMBREE
	RTCRay4 ray;
	for(uint i = 0; i < BLCLOUD_VSIZE; ++i){
		dfloat4 RO = dfloat4(ro.get(i));
		ray.orgx[i] = RO.x;
		ray.orgy[i] = RO.y;
		ray.orgz[i] = RO.z;
		dfloat4 RD = dfloat4(rd.get(i));
		ray.dirx[i] = RD.x;
		ray.diry[i] = RD.y;
		ray.dirz[i] = RD.z;
		ray.tnear[i] = 0.0f;
		ray.tfar[i] = MAX_OCCLUSION_DIST;
		ray.mask[i] = -1;
		ray.time[i] = 0;
		ray.instID[i] = RTC_INVALID_GEOMETRY_ID;
		ray.geomID[i] = RTC_INVALID_GEOMETRY_ID;
		ray.primID[i] = RTC_INVALID_GEOMETRY_ID;
	}

	rtcIntersect4(&gm.v,pscene,ray);

	for(uint i = 0; i < BLCLOUD_VSIZE; ++i){
		pmask->v[i] = ray.geomID[i] != RTC_INVALID_GEOMETRY_ID?-1:0;
		pdist->v[i] = ray.tfar[i];
	}
#else
	*pmask = dintN(0);
#endif
}

void SceneOcclusion::Destroy(){
#ifdef USE_EMBREE
	//TODO: delete the geometry?
	rtcDeleteScene(pscene);
	rtcDeleteDevice(pdev);
#endif
}
