#ifndef SCENE_H
#define SCENE_H

enum SCENE_CACHE_MODE{
	SCENE_CACHE_DISABLED,
	SCENE_CACHE_READ,
	SCENE_CACHE_WRITE
};

enum VOLUME_BUFFER{
	VOLUME_BUFFER_SDF,
	VOLUME_BUFFER_FOG,
	VOLUME_BUFFER_COUNT
};

class BoundingBox{
public:
	BoundingBox();
	BoundingBox(const float4 &, const float4 &);
	~BoundingBox();
	bool Intersects(const float4 &, const float4 &, const float4 &) const;
	bool Intersects(const BoundingBox &) const;
	dfloat3 sc;
	dfloat3 se;
};

class OctreeStructure{
public:
	OctreeStructure();
	~OctreeStructure();
	dfloat4 ce; //(center.xyz,extent)
	uint chn[8];
	uint volx[VOLUME_BUFFER_COUNT]; //leaf volume index
	//
	//min/max query value to speed up rendering; sdf: min distance, fog: max density
	float qval[VOLUME_BUFFER_COUNT];
};

class Octree{
public:
	Octree(uint);
	Octree();
	~Octree();
	void BuildPath(const float4 &, const float4 &, const float4 &, const float4 &, uint, uint, std::atomic<uint> *, std::atomic<uint> *, tbb::concurrent_vector<Octree> *, tbb::concurrent_vector<OctreeStructure> *, VOLUME_BUFFER);
	void BuildPath(const float4 &, const float4 &, const float4 &, const float4 &, const float4 &, uint, uint, std::atomic<uint> *, std::atomic<uint> *, tbb::concurrent_vector<Octree> *, tbb::concurrent_vector<OctreeStructure> *, VOLUME_BUFFER);
	void FreeRecursive();
	Octree *pch[8];
	uint x; //node index
	//std::atomic_flag lock; //MT node write-access
	tbb::spin_mutex m;
};

namespace Node{
class NodeTree;
}

namespace SceneData{

class BaseObject{
public:
	BaseObject(Node::NodeTree *);
	virtual ~BaseObject();
	class Node::NodeTree *pnt;
};

class ParticleSystem : public BaseObject{
public:
	ParticleSystem(Node::NodeTree *);
	~ParticleSystem();
	static void DeleteAll();
	std::vector<dfloat3> pl; //position
	std::vector<dfloat3> vl; //velocity
	static std::vector<ParticleSystem *> prss;
};

class SmokeCache : public BaseObject{
public:
	SmokeCache(Node::NodeTree *);
	~SmokeCache();
	static void DeleteAll();
	static std::vector<SmokeCache *> objs;
};

class Surface : public BaseObject{
public:
	Surface(Node::NodeTree *);
	~Surface();
	static void DeleteAll();
	std::vector<dfloat3> vl;
	std::vector<uint> tl;
	static std::vector<Surface *> objs;
};

}

#ifdef OPENVDB_OPENVDB_HAS_BEEN_INCLUDED
typedef openvdb::tools::GridSampler<openvdb::FloatGrid, openvdb::tools::BoxSampler> FloatGridBoxSampler;
typedef openvdb::tools::GridSampler<openvdb::Vec3SGrid, openvdb::tools::BoxSampler> VectorGridBoxSampler;

namespace Node{

using InputNodeParams = std::tuple<SceneData::BaseObject *, openvdb::math::Transform::Ptr, const FloatGridBoxSampler *, const FloatGridBoxSampler *, const FloatGridBoxSampler *, const VectorGridBoxSampler *>;
enum INP{
	INP_OBJECT,
	INP_TRANSFORM,
	INP_SDFSAMPLER,
	INP_QGRSAMPLER,
	INP_FOGSAMPLER,
	INP_VELSAMPLER,
};

class ValueNodeParams : public IValueNodeParams{
public:
	ValueNodeParams(const dfloat3 *, const dfloat3 *, float, float, const dfloat3 *, const FloatGridBoxSampler *[VOLUME_BUFFER_COUNT], const FloatGridBoxSampler *, const VectorGridBoxSampler *);
	~ValueNodeParams();
	const dfloat3 * GetVoxPosW() const;
	const dfloat3 * GetCptPosW() const;
	float GetLocalDistance() const;
	float GetLocalDensity() const;
	const dfloat3 * GetVoxPosWAdv() const;
	float SampleGlobalDistance(const dfloat3 &, bool) const;
	float SampleGlobalDensity(const dfloat3 &) const;
	dfloat3 SampleGlobalVector(const dfloat3 &) const;
	//global
	const FloatGridBoxSampler *psampler[VOLUME_BUFFER_COUNT];
	const FloatGridBoxSampler *pqsampler;
	const VectorGridBoxSampler *pvsampler;
	//local
	const dfloat3 *pvoxw;
	const dfloat3 *pcptw;
	float distance;
	float density;
	//advection
	const dfloat3 *pvoxwa;
};

}
#endif

class Scene{
public:
	Scene();
	~Scene();
	void Initialize(float, uint, float, SCENE_CACHE_MODE);
	void Destroy();
	float *pvol[VOLUME_BUFFER_COUNT];
	uint lvoxc;
	uint index;
	uint leafx[VOLUME_BUFFER_COUNT];
	uint lvoxc3;
	tbb::concurrent_vector<Octree> root;
	tbb::concurrent_vector<OctreeStructure> ob;
};

#endif
