#ifndef KERNEL_OCTREE_H
#define KERNEL_OCTREE_H

namespace KernelOctree{

class BaseOctreeTraverser{
public:
	BaseOctreeTraverser();
	~BaseOctreeTraverser();
	virtual void Initialize(const sfloat4 &, const sfloat4 &, const dintN &, const dfloatN &, const tbb::concurrent_vector<OctreeStructure> *) = 0;
	virtual dintN GetLeaf(uint, duintN *, sfloat1 &, sfloat1 &) = 0;
};

class OctreeFullTraverser : public BaseOctreeTraverser{
public:
	OctreeFullTraverser();
	~OctreeFullTraverser();
	void Initialize(const sfloat4 &, const sfloat4 &, const dintN &, const dfloatN &, const tbb::concurrent_vector<OctreeStructure> *); //traverse full path
	dintN GetLeaf(uint, duintN *, sfloat1 &, sfloat1 &);

	typedef std::tuple<uint, float, float> Node;
	std::vector<Node> ls[BLCLOUD_VSIZE];
};

class OctreeStepTraverser : public BaseOctreeTraverser{
public:
	OctreeStepTraverser();
	~OctreeStepTraverser();
	void Initialize(const sfloat4 &, const sfloat4 &, const dintN &, const dfloatN &, const tbb::concurrent_vector<OctreeStructure> *); //do nothing, just store the ray
	dintN GetLeaf(uint, duintN *, sfloat1 &, sfloat1 &);
};

}

#endif
