#ifndef NODE_H
#define NODE_H

namespace Node{

class IValueNodeParams{
public:
	IValueNodeParams();
	~IValueNodeParams();
	virtual const dfloat3 * GetObjectPosW() const = 0;
	virtual const dfloat3 * GetVoxPosW() const = 0;
	virtual const dfloat3 * GetCptPosW() const = 0;
	virtual float GetLocalDistance() const = 0;
	virtual float GetLocalDensity() const = 0;
	virtual const dfloat3 * GetVoxPosWAdv() const = 0;
	virtual float GetAdvectionDistance() const = 0;
	virtual float GetAdvectionDensity() const = 0;
	virtual float SampleGlobalDistance(const dfloat3 &, bool) const = 0;
	virtual float SampleGlobalDensity(const dfloat3 &) const = 0;
	virtual dfloat3 SampleGlobalVector(const dfloat3 &) const = 0;
	virtual dfloat3 SampleGlobalGradient(const dfloat3 &) const = 0;
};

class NodeTree;

class BaseNode{
protected:
	BaseNode(uint, NodeTree *);
public:
	virtual ~BaseNode();
	virtual void Evaluate(const void *) = 0;
	virtual void Clear();
	BaseNode *pnodes[12];
	NodeTree *pntree; //can be null
	uint indices[12]; //index to the input node (pnodes[x]) output socket, per-socket-type basis
	uint imask; //input mask - a bitmask to indicate if any (real) input nodes are connected (not using default automatic base node)
	uint omask; //output mask to help optimize storage and performance in some cases
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

class FloatInput : public BaseValueNode<float>{
public:
	FloatInput(uint, NodeTree *);
	~FloatInput();
	void Evaluate(const void *);
};

class ScalarMath : public BaseValueNode<float>{
public:
	ScalarMath(uint, NodeTree *, char);
	~ScalarMath();
	void Evaluate(const void *);
	char opch;
};

class VectorInput : public BaseValueNode<dfloat3>{
public:
	VectorInput(uint, NodeTree *);
	~VectorInput();
	void Evaluate(const void *);
};

class VectorMath : public BaseValueNode<dfloat3>{
public:
	VectorMath(uint, NodeTree *, char);
	~VectorMath();
	void Evaluate(const void *);
	char opch;
};

class VectorMix : public BaseValueNode<dfloat3>{
public:
	VectorMix(uint, NodeTree *);
	~VectorMix();
	void Evaluate(const void *);
	enum OUTPUT_FLOAT{
		OUTPUT_FLOAT_X,
		OUTPUT_FLOAT_Y,
		OUTPUT_FLOAT_Z,
		OUTPUT_FLOAT_COUNT
	};
};

class VectorXYZ : public BaseValueNode<float>{
public:
	VectorXYZ(uint, NodeTree *);
	~VectorXYZ();
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

class IVoronoiLayers : public virtual BaseValueNode<float>{
protected:
	IVoronoiLayers(uint, NodeTree *);
	~IVoronoiLayers();
public:
	virtual void Evaluate(const void *);
	static IVoronoiLayers * Create(uint, NodeTree *);
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
};

class BaseFogNode : public virtual BaseNode{
protected:
	BaseFogNode(uint, NodeTree *);
	virtual ~BaseFogNode();
public:
	virtual void Evaluate(const void *);
	virtual void Clear();
	static BaseFogNode * Create(uint, NodeTree *);
};

class BaseSurfaceNode : public virtual BaseNode{
protected:
	//BaseSurfaceNode();
	BaseSurfaceNode(uint, NodeTree *);
	virtual ~BaseSurfaceNode();
public:
	virtual void Evaluate(const void *);
	virtual void Clear();
	static BaseSurfaceNode * Create(uint, NodeTree *);
};

class BaseVectorFieldNode : public virtual BaseNode{
protected:
	BaseVectorFieldNode(uint, NodeTree *);
	virtual ~BaseVectorFieldNode();
public:
	virtual void Evaluate(const void *);
	virtual void Clear();
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

class AdvectionInfo : public BaseValueNode<float>, public BaseValueNode<dfloat3>{
public:
	AdvectionInfo(uint, NodeTree *);
	~AdvectionInfo();
	void Evaluate(const void *);
	enum OUTPUT_VECTOR{
		OUTPUT_VECTOR_VOXPOSW,
		OUTPUT_VECTOR_COUNT
	};
	enum OUTPUT_FLOAT{
		OUTPUT_FLOAT_ADVDISTANCE,
		OUTPUT_FLOAT_DENSITY,
		OUTPUT_FLOAT_COUNT
	};
};

/*struct ObjectInfoData{
	dfloat3 location;
};*/

class ObjectInfo : public BaseValueNode<dfloat3>{
public:
	ObjectInfo(uint, NodeTree *);
	~ObjectInfo();
	void Evaluate(const void *);
	enum OUTPUT_VECTOR{
		OUTPUT_VECTOR_LOCATION,
		OUTPUT_VECTOR_COUNT
	};
};

class SceneInfo : public BaseValueNode<float>, public BaseValueNode<dfloat3>{
public:
	SceneInfo(uint, NodeTree *);
	~SceneInfo();
	void Evaluate(const void *);
	enum INPUT{
		INPUT_POSITION
	};
	enum OUTPUT_FLOAT{
		OUTPUT_FLOAT_DISTANCE,
		OUTPUT_FLOAT_SURFACE,
		OUTPUT_FLOAT_DENSITY,
		OUTPUT_FLOAT_FINAL,
		OUTPUT_FLOAT_COUNT
	};
	enum OUTPUT_VECTOR{
		OUTPUT_VECTOR_VECTOR,
		OUTPUT_VECTOR_GRADIENT,
		OUTPUT_VECTOR_COUNT
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
		INPUT_SIZE,
		INPUT_CUTOFF,
		INPUT_COUNT
	};
};

class IFieldInput : public virtual BaseFogNode, public virtual BaseVectorFieldNode{
protected:
	IFieldInput(uint, NodeTree *);
	~IFieldInput();
public:
	virtual void Evaluate(const void *) = 0;
	virtual void Clear() = 0;
	static IFieldInput * Create(uint, NodeTree *);
	enum INPUT{
		INPUT_RASTERIZATIONRES,
		INPUT_WEIGHT
	};
	enum OUTPUT_FOG{
		OUTPUT_FOG_DENSITY,
	};
	enum OUTPUT_VECTORFIELD{
		OUTPUT_VECTORFIELD_VELOCITY
	};
};

class ISmokeCache : public virtual BaseFogNode, public virtual BaseVectorFieldNode{
protected:
	ISmokeCache(uint, NodeTree *);
	~ISmokeCache();
public:
	virtual void Evaluate(const void *) = 0;
	virtual void Clear() = 0;
	static ISmokeCache * Create(uint, NodeTree *);
	enum OUTPUT_FOG{
		OUTPUT_FOG_DENSITY,
	};
	enum OUTPUT_VECTORFIELD{
		OUTPUT_VECTORFIELD_VELOCITY
	};
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

class ICombine : public virtual BaseFogNode{
protected:
	ICombine(uint, NodeTree *);
	~ICombine();
public:
	virtual void Evaluate(const void *) = 0;
	static ICombine * Create(uint, NodeTree *, char);
	enum INPUT{
		INPUT_FOGA,
		INPUT_FOGB,
		INPUT_COUNT
	};
};

class IAdvection : public virtual BaseFogNode{
protected:
	IAdvection(uint, NodeTree *);
	~IAdvection();
public:
	virtual void Evaluate(const void *) = 0;
	static IAdvection * Create(uint, NodeTree *, uint);
	enum INPUT{
		INPUT_THRESHOLD,
		INPUT_DISTANCE,
		INPUT_ITERATIONS,
		INPUT_DENSITY,
		INPUT_VELOCITY,
		INPUT_FOG,
		INPUT_COUNT
	};
	enum BOOL{
		BOOL_BREAK_ITERATION,
		BOOL_SAMPLE_LOCAL
	};
};

class ISurfaceToFog : public virtual BaseFogNode{
protected:
	ISurfaceToFog(uint, NodeTree *);
	~ISurfaceToFog();
public:
	virtual void Evaluate(const void *) = 0;
	static ISurfaceToFog * Create(uint, NodeTree *, float);
	enum INPUT{
		INPUT_SURFACE,
		INPUT_COUNT
	};
};

class IDisplacement : public virtual BaseSurfaceNode{
protected:
	IDisplacement(uint, NodeTree *);
	~IDisplacement();
public:
	virtual void Evaluate(const void *) = 0;
	static IDisplacement * Create(uint, NodeTree *, float);
	enum INPUT{
		INPUT_DISTANCE,
		INPUT_MAXIMUM,
		INPUT_BILLOW,
		INPUT_SURFACE,
		INPUT_COUNT
	};
};

class ICSG : public virtual BaseSurfaceNode{
protected:
	ICSG(uint, NodeTree *);
	~ICSG();
public:
	virtual void Evaluate(const void *) = 0;
	static ICSG * Create(uint, NodeTree *, char);
	enum INPUT{
		INPUT_SURFACEA,
		INPUT_SURFACEB,
		INPUT_COUNT
	};
};

class OutputNode : public BaseNode{
public:
	OutputNode(NodeTree *, char, bool);
	~OutputNode();
	void Evaluate(const void *);
	enum INPUT{
		INPUT_FOG,
		INPUT_FOGPOST,
		INPUT_VECTOR,
		INPUT_SURFACE,
		INPUT_COUNT,
	};
	char opch; //fog post processor blend operation
	bool qonly; //construct only query field
};

class NodeTree{
public:
	NodeTree(const char *);
	~NodeTree();
	void EvaluateNodes0(const void *, uint, uint);
	void EvaluateNodes1(const void *, uint, uint);
	//void Cleanup();
	void ApplyBranchMask();
	void SortNodes();
	BaseNode * GetRoot() const;
	static void DeleteAll();
	std::vector<BaseNode *> nodes0; //low-level nodes (math, info nodes, values etc)
	std::vector<BaseNode *> nodes1; //high-level nodes (surface and fog operations)
	char name[256];
	static std::vector<NodeTree *> ntrees;
};

BaseNode * CreateNodeByType(const char *, const void *, uint, NodeTree *);
BaseNode * CreateNodeBySocket(const char *, const void *, uint, NodeTree *);

}

#endif
