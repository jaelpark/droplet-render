#ifndef SCENE_H
#define SCENE_H

#define BLCLOUD_uN 8

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

struct OctreeStructure{ //GPU octree
	//XMFLOAT3 c;
    dfloat4 ce; //(center.xyz,extent)
	uint chn[8];
	uint volx[VOLUME_BUFFER_COUNT]; //leaf volume index
	float pmax; //max density
};

class Octree{
public:
	Octree(uint);
	~Octree();
	void BuildPath(const float4 &, const float4 &, const float4 &, const float4 &, uint, uint, std::atomic<uint> *, std::atomic<uint> *, Octree *, OctreeStructure *, VOLUME_BUFFER);
    void BuildPath(const float4 &, const float4 &, const float4 &, const float4 &, const float4 &, uint, uint, std::atomic<uint> *, std::atomic<uint> *, Octree *, OctreeStructure *, VOLUME_BUFFER);
	Octree *pch[8];
	uint x; //node index
    std::atomic_flag lock; //MT node write-access
};

//CPU version works with array of volumes
struct LeafVolume{
    float pvol[BLCLOUD_uN*BLCLOUD_uN*BLCLOUD_uN]; //pdst
	//float max; //should probably go to the OctreeStructure, to let this correspond more to the actual 3d volume texture
    /*
       TODO: Two separate LeafVolume arrays, one for distance and another for density. This allows resolutions to be independent
       of each other. Additionally one can allocate one type of leaf while the other doesn't exist. This requires two volx pointer
       indices. Whenever there's a need to allocate a leaf for just one type of voxels, the other volx points to an empty section (zeros)
       and the end of the voxel buffer array.
    */
};

namespace Node{
class NodeTree;
}

class ParticleSystem{
public:
    ParticleSystem(Node::NodeTree *);
    ~ParticleSystem();
    static void DeleteAll();
	class Node::NodeTree *pnt;
    std::vector<dfloat3> vl;
    static std::vector<ParticleSystem *> prss;
};

class SceneObject{
public:
    SceneObject(Node::NodeTree *);
    ~SceneObject();
    static void DeleteAll();
    class Node::NodeTree *pnt;
    std::vector<dfloat3> vl;
    std::vector<uint> tl;
    static std::vector<SceneObject *> objs;
};

//scene builder
class Scene{
public:
	Scene();
	~Scene();
    void Initialize(float, SCENE_CACHE_MODE);
	//void AddObject();
	//void BuildScene();
	void Destroy();
	OctreeStructure *pob;
    LeafVolume *pbuf[VOLUME_BUFFER_COUNT]; //-> psdfb, pfogb
    uint index;
    uint leafx;
    //uint octreesize;
    /*enum CACHE_MODE{
        //
    };*/
};

//octree

#endif
