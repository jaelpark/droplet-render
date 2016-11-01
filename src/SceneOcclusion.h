#ifndef SCENE_OCCLUSION_H
#define SCENE_OCCLUSION_H

#define MAX_OCCLUSION_DIST 1e6f

class SceneOcclusion{
public:
	SceneOcclusion();
	~SceneOcclusion();
	void Initialize();
	sint1 Intersect(const sfloat4 &, const sfloat4 &, const sfloat1 &, sfloat1 &) const;
	void Destroy();
#ifdef USE_EMBREE
	struct __RTCDevice *pdev;
	struct __RTCScene *pscene;
#endif
};

#endif
