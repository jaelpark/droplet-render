#include "main.h"

#ifdef USE_EMBREE
#include <embree2/rtcore.h>
#include <embree2/rtcore_ray.h>
#endif

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

	for(uint i = 0; i < SceneData::Surface::objs.size(); ++i){
		if(!SceneData::Surface::objs[i]->flags & SCENEOBJ_HOLDOUT)
			continue;
		uint m = rtcNewTriangleMesh(pscene,RTC_GEOMETRY_STATIC,SceneData::Surface::objs[i]->tl.size()/3,SceneData::Surface::objs[i]->vl.size());

		struct Vertex{
			float x, y, z, w;
		} *pvs = (Vertex*)rtcMapBuffer(pscene,m,RTC_VERTEX_BUFFER);
		for(uint j = 0, n = SceneData::Surface::objs[i]->vl.size(); j < n; ++j){
			dfloat3 &v = SceneData::Surface::objs[i]->vl[j];
			pvs[j].x = v.x;
			pvs[j].y = v.y;
			pvs[j].z = v.z;
		}
		rtcUnmapBuffer(pscene,m,RTC_VERTEX_BUFFER);

		struct Triangle{
			int v0, v1, v2;
		} *pts = (Triangle*)rtcMapBuffer(pscene,m,RTC_INDEX_BUFFER);
		for(uint j = 0, n = SceneData::Surface::objs[i]->tl.size()/3; j < n; ++j){
			pts[j].v0 = SceneData::Surface::objs[i]->tl[3*j+0];
			pts[j].v1 = SceneData::Surface::objs[i]->tl[3*j+1];
			pts[j].v2 = SceneData::Surface::objs[i]->tl[3*j+2];
		}
		rtcUnmapBuffer(pscene,m,RTC_INDEX_BUFFER);
	}

	rtcCommit(pscene);
#endif
}

sint1 SceneOcclusion::Intersect(const sfloat4 &ro, const sfloat4 &rd, const sfloat1 &gm, sfloat1 &dist) const{
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

	dintN mask;
	for(uint i = 0; i < BLCLOUD_VSIZE; ++i){
		mask.v[i] = ray.geomID[i] != RTC_INVALID_GEOMETRY_ID?-1:0;
		((float*)&dist.v)[i] = ray.tfar[i];
	}
	return sint1::load(&mask);
#else
	return sint1(0);
#endif
}

void SceneOcclusion::Destroy(){
#ifdef USE_EMBREE
	//TODO: delete the geometry?
	rtcDeleteScene(pscene);
	rtcDeleteDevice(pdev);
#endif
}
