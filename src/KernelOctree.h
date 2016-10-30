#ifndef KERNEL_OCTREE_H
#define KERNEL_OCTREE_H

class BaseOctreeTraverser{
public:
	BaseOctreeTraverser();
	~BaseOctreeTraverser();
	virtual void Initialize(const sfloat4 &, const sfloat4 &, const dintN &, const dfloatN &, const tbb::concurrent_vector<OctreeStructure> *) = 0;
	virtual void GetLeaf(uint, uint, uint *, float *, float *) = 0;
};

class OctreeFullTraverser : public BaseOctreeTraverser{
public:
	OctreeFullTraverser();
	~OctreeFullTraverser();
	void Initialize(const sfloat4 &, const sfloat4 &, const dintN &, const dfloatN &, const tbb::concurrent_vector<OctreeStructure> *); //traverse full path
	void GetLeaf(uint, uint, uint *, float *, float *);

	typedef std::tuple<uint, float, float> Node;
	std::vector<Node> ls[BLCLOUD_VSIZE];
};

class OctreeStepTraverser : public BaseOctreeTraverser{
public:
	OctreeStepTraverser();
	~OctreeStepTraverser();
	void Initialize(const sfloat4 &, const sfloat4 &, const dintN &, const dfloatN &, const tbb::concurrent_vector<OctreeStructure> *); //do nothing, just store the ray
	void GetLeaf(uint, uint, uint *, float *, float *);
};

#endif
