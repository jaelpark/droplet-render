#ifndef KERNEL_OCTREE_H
#define KERNEL_OCTREE_H

namespace KernelOctree{

class BaseOctreeTraverser{
public:
	BaseOctreeTraverser();
	~BaseOctreeTraverser();
	virtual void Initialize(const sfloat4 &, const sfloat4 &, const dintN &, const dfloatN &, const tbb::concurrent_vector<OctreeStructure> *) = 0;
	virtual dintN GetLeaf(uint, duintN *, sfloat1 &, sfloat1 &) = 0;
protected:
	typedef std::tuple<uint, float, float> Node;
	const tbb::concurrent_vector<OctreeStructure> *pob;
};

class OctreeFullTraverser : public BaseOctreeTraverser{
public:
	OctreeFullTraverser();
	~OctreeFullTraverser();
	void Initialize(const sfloat4 &, const sfloat4 &, const dintN &, const dfloatN &, const tbb::concurrent_vector<OctreeStructure> *);
	dintN GetLeaf(uint, duintN *, sfloat1 &, sfloat1 &);
private:
	std::vector<Node> ls[BLCLOUD_VSIZE];
	bool OctreeProcessSubtree(const dfloat3 &, const dfloat3 &, uint, uint, uint, float, std::vector<Node> *);
};

class OctreeStepTraverser : public BaseOctreeTraverser{
public:
	OctreeStepTraverser();
	~OctreeStepTraverser();
	void Initialize(const sfloat4 &, const sfloat4 &, const dintN &, const dfloatN &, const tbb::concurrent_vector<OctreeStructure> *);
	dintN GetLeaf(uint, duintN *, sfloat1 &, sfloat1 &);
private:
#define MAX_DEPTH 16
	typedef struct{
		dfloat3 t0[MAX_DEPTH], t1[MAX_DEPTH], tm[MAX_DEPTH];
		uint n[MAX_DEPTH], tn[MAX_DEPTH], a, l;
		bool p[MAX_DEPTH]; //stack pointer
	} Stack;
#undef MAX_DEPTH
	Stack stack[BLCLOUD_VSIZE];
	dintN mask;
	dfloatN maxd;
	bool OctreeProcessSubtree(dfloat3 *, dfloat3 *, dfloat3 *, uint *, uint *, bool *, uint &, uint, float, Node *);
};

}

#endif
