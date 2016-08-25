#ifndef NODE_H
#define NODE_H

//class PyObject;

namespace Node{

using ValueNodeParams = std::tuple<dfloat3 *, dfloat3 *, float, float>;
enum VNP{
	VNP_VOXPOSW,
	VNP_CPTPOSW,
	VNP_DISTANCE,
	VNP_DENSITY
};

class NodeTree;

class BaseNode{
//protected: //should be protected
    //BaseNode();
protected:
    BaseNode(uint, NodeTree *);
	//BaseNode();
public:
    virtual ~BaseNode();
    //virtual void SetConstants(const void *);
    virtual void Evaluate(const void *) = 0;
    //virtual BaseNode * NewNode(void *, uint) const = 0;
    BaseNode *pnodes[12];
    NodeTree *pntree; //can be null
	uint indices[12]; //index to the input node (pnodes[x]) output socket, per-socket-type basis
	uint imask; //input mask - a bitmask to indicate if any input nodes are connected (not using default base node)
    uint omask; //output mask to help optimize storage in some cases
    uint emask; //ouput root node branch mask (e.g. 0x1 if required for the 1st input, 0x2 for the second, 0x1|0x2 for both, etc.)
    uint level;
};

template<class T>
class BaseValueResult{
public:
	BaseValueResult(const T &); //used by tbb to set the default when constructing BaseValueNode
	BaseValueResult();
	T value[4]; //values for the different output sockets
};

template<class T>
class BaseValueNode : public virtual BaseNode{
//NOTE: the virtual inheritance requires an explicit constructor call at higher levels of hierarchy also whenever BaseValueNode<T>(...) is called.
public:
    BaseValueNode(T, uint, NodeTree *);
    BaseValueNode(uint, NodeTree *);
    virtual ~BaseValueNode();
    virtual void Evaluate(const void *);
	inline T locr(uint outx){
		return this->BaseValueNode<T>::result.local().value[outx];
	}
	tbb::enumerable_thread_specific<BaseValueResult<T>> result;
};

//NOTE: this could all be replaced with some common arithmetic node implementation which calls some operand
template<class T>
class AddNode : public BaseValueNode<T>{
public:
    AddNode(uint, NodeTree *);
    ~AddNode();
    void Evaluate(const void *);
};

template<class T>
class SubNode : public BaseValueNode<T>{
public:
    SubNode(uint, NodeTree *);
    ~SubNode();
    void Evaluate(const void *);
};

template<class T>
class MulNode : public BaseValueNode<T>{
public:
    MulNode(uint, NodeTree *);
    ~MulNode();
    void Evaluate(const void *);
    //BaseNode * NewNode() const;
    //T result;
};

template<class T>
class DivNode : public BaseValueNode<T>{
public:
    DivNode(uint, NodeTree *);
    ~DivNode();
    void Evaluate(const void *);
};

template<class T>
class PowNode : public BaseValueNode<T>{
public:
    PowNode(uint, NodeTree *);
    ~PowNode();
    void Evaluate(const void *);
};

template<class T>
class MinNode : public BaseValueNode<T>{
public:
	MinNode(uint, NodeTree *);
	~MinNode();
	void Evaluate(const void *);
};

template<class T>
class MaxNode : public BaseValueNode<T>{
public:
	MaxNode(uint, NodeTree *);
	~MaxNode();
	void Evaluate(const void *);
};

class IFbmNoise : public virtual BaseValueNode<float>, public virtual BaseValueNode<dfloat3>{
protected:
	IFbmNoise(uint, NodeTree *);
	~IFbmNoise();
public:
	virtual void Evaluate(const void *);
	static IFbmNoise * Create(uint, NodeTree *);
	enum INPUT{
        INPUT_OCTAVES,
        INPUT_FREQ,
        INPUT_AMP,
        INPUT_FJUMP,
        INPUT_GAIN,
        INPUT_POSITION,
        INPUT_COUNT
    };
	enum OUTPUT_FLOAT{
		OUTPUT_FLOAT_NOISE,
		OUTPUT_FLOAT_MAXIMUM,
		OUTPUT_FLOAT_COUNT
	};
	enum OUTPUT_VECTOR{
		OUTPUT_VECTOR_NOISE,
		OUTPUT_VECTOR_COUNT
	};
};

class BaseFogNode : public BaseNode{
protected:
    BaseFogNode(uint, NodeTree *);
    virtual ~BaseFogNode();
public:
    virtual void Evaluate(const void *);
    static BaseFogNode * Create(uint, NodeTree *);
    static void SortNodes();
    static void EvaluateAll(const void *, uint);
};

class BaseSurfaceNode : public BaseNode{
protected:
    //BaseSurfaceNode();
    BaseSurfaceNode(uint, NodeTree *);
    virtual ~BaseSurfaceNode();
public:
    virtual void Evaluate(const void *);
    static BaseSurfaceNode * Create(uint, NodeTree *);
};

class BaseVectorFieldNode : public BaseNode{
protected:
	BaseVectorFieldNode(uint, NodeTree *);
	virtual ~BaseVectorFieldNode();
public:
	virtual void Evaluate(const void *);
	static BaseVectorFieldNode * Create(uint, NodeTree *);
};

class VoxelInfo : public BaseValueNode<float>, public BaseValueNode<dfloat3>{
public:
	VoxelInfo(uint, NodeTree *);
	~VoxelInfo();
	void Evaluate(const void *);
	enum OUTPUT_VECTOR{
		OUTPUT_VECTOR_VOXPOSW,
		OUTPUT_VECTOR_CPTPOSW,
		OUTPUT_VECTOR_COUNT
	};
	enum OUTPUT_FLOAT{
		OUTPUT_FLOAT_DISTANCE,
		OUTPUT_FLOAT_DENSITY,
		OUTPUT_FLOAT_COUNT
	};
};

class ISurfaceInput : public virtual BaseSurfaceNode{
protected:
    ISurfaceInput(uint, NodeTree *);
    ~ISurfaceInput();
public:
    virtual void Evaluate(const void *) = 0;
    static ISurfaceInput * Create(uint, NodeTree *);
};

class IParticleInput : public virtual BaseFogNode{
protected:
	IParticleInput(uint, NodeTree *);
	~IParticleInput();
public:
	virtual void Evaluate(const void *) = 0;
	static IParticleInput * Create(uint, NodeTree *);
	enum INPUT{
		/*INPUT_RASTERIZATIONRES,
		INPUT_WEIGHT,*/
		INPUT_SIZE,
		INPUT_CUTOFF,
		INPUT_COUNT
	};
};

class ISmokeCache : public virtual BaseFogNode{
protected:
	ISmokeCache(uint, NodeTree *);
	~ISmokeCache();
public:
	virtual void Evaluate(const void *) = 0;
	static ISmokeCache * Create(uint, NodeTree *);
};

class IFogPostInput : public virtual BaseFogNode{
protected:
	IFogPostInput(uint, NodeTree *);
	~IFogPostInput();
public:
	virtual void Evaluate(const void *) = 0;
	static IFogPostInput * Create(uint, NodeTree *);
};

class IComposite : public virtual BaseFogNode{
protected:
	IComposite(uint, NodeTree *);
	~IComposite();
public:
	virtual void Evaluate(const void *) = 0;
	static IComposite * Create(uint, NodeTree *);
	enum INPUT{
		INPUT_VALUE,
		INPUT_FOG,
		INPUT_COUNT
	};
};

class IAdvection : public virtual BaseFogNode{
protected:
	IAdvection(uint, NodeTree *);
	~IAdvection();
public:
	virtual void Evaluate(const void *) = 0;
	static IAdvection * Create(uint, NodeTree *);
	enum INPUT{
		INPUT_THRESHOLD,
		INPUT_DISTANCE,
		INPUT_ITERATIONS,
		INPUT_VELOCITY,
		INPUT_FOG,
		INPUT_COUNT
	};
};

class IDisplacement : public virtual BaseSurfaceNode{
protected:
    IDisplacement(uint, NodeTree *);
    ~IDisplacement();
public:
    virtual void Evaluate(const void *) = 0;
    static IDisplacement * Create(uint, NodeTree *);
    enum INPUT{
        INPUT_DISTANCE,
		INPUT_MAXIMUM,
		INPUT_BILLOW,
        INPUT_SURFACE,
        INPUT_COUNT
    };
};

#ifdef BLCLOUD_DEPRECATED
class IfBmPerlinNoise : public virtual BaseSurfaceNode{
protected:
    IfBmPerlinNoise(uint, NodeTree *);
    ~IfBmPerlinNoise();
public:
    virtual void Evaluate(const void *) = 0;
    static IfBmPerlinNoise * Create(uint, NodeTree *);
    enum INPUT{
        INPUT_OCTAVES,
        INPUT_FREQ,
        INPUT_AMP,
        INPUT_FJUMP,
        INPUT_GAIN,
        INPUT_BILLOW,
        INPUT_SURFACE,
        INPUT_COUNT
    };
};
#endif

class IVectorFieldSampler : public virtual BaseValueNode<dfloat3>{
protected:
	IVectorFieldSampler(uint, NodeTree *);
	~IVectorFieldSampler();
public:
	virtual void Evaluate(const void *) = 0;
	static IVectorFieldSampler * Create(uint, NodeTree *);
	enum INPUT{
		INPUT_FIELD,
		INPUT_POSITION
	};
};

//fake proxy nodes to split material output tree into subtrees
/*class SurfaceOutputNode : public BaseNode{
public:
    SurfaceOutputNode(NodeTree *);
    ~SurfaceOutputNode();
    void Evaluate(const void *);
    enum INPUT{
        INPUT_SURFACE,
        INPUT_COUNT
    };
};*/

class OutputNode : public BaseNode{
public:
    OutputNode(NodeTree *);
    ~OutputNode();
    //BaseNode * NewNode(void *, uint) const;
    void Evaluate(const void *);
    enum INPUT{
        //INPUT_SHADER,
        //INPUT_FIELD,
        INPUT_FOG,
		INPUT_FOGPOST,
        //INPUT_FIELD,
        INPUT_SURFACE,
        INPUT_COUNT,
    };
    //static OutputNode *proot; //global root node
};

class NodeTree{
public:
    NodeTree();
    ~NodeTree();
    void EvaluateNodes0(const void *, uint, uint);
    void EvaluateNodes1(const void *, uint, uint);
    //void ApplyBranchMask(BaseNode *, uint);
    void ApplyBranchMask();
    void SortNodes();
    BaseNode * GetRoot() const;
    static void DeleteAll();
    std::vector<BaseNode *> nodes0; //low
    std::vector<BaseNode *> nodes1; //high
    static std::vector<NodeTree *> ntrees;
    //OutputNode *proot; //hnodes.back()
};

//Instead of each class holding its own global node vector, one could use this to manage evaluation groups.
//Pass it to constructor or something.
/*class NodeEvaluationGroup{
public:
    NodeEvaluationGroup();
    ~NodeEvaluationGroup();
    //Collect()
    std::vector<BaseNode *> nodes;
};*/

BaseNode * CreateNodeByType(const char *, const void *, uint, NodeTree *);
BaseNode * CreateNodeBySocket(const char *, const void *, uint, NodeTree *);
//void EvaluateValueGroup(uint);
//void SortNodes();
//void DeleteNodes();

}

#endif
