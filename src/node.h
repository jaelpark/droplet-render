#ifndef NODE_H
#define NODE_H

//class PyObject;

namespace Node{

/*template<class T>
struct PropertyDef{
    const char *pname;
    T *pp;
};*/

class NodeTree;

class BaseNode{
//protected: //should be protected
    //BaseNode();
protected:
    BaseNode(uint, NodeTree *);
public:
    virtual ~BaseNode();
    //virtual void SetConstants(const void *);
    virtual void Evaluate(const void *) = 0;
    //virtual BaseNode * NewNode(void *, uint) const = 0;
    BaseNode *pnodes[12];
    NodeTree *pntree; //can be null
    uint omask; //output mask to help optimize storage in some cases
    uint emask; //ouput root node branch mask (e.g. 0x1 if required for the 1st input, 0x2 for the second, 0x1|0x2 for both, etc.)
    uint level;
};

template<class T>
class BaseValueNode : public BaseNode{
public:
    BaseValueNode(T, uint, NodeTree *);
    BaseValueNode(uint, NodeTree *);
    virtual ~BaseValueNode();
    virtual void Evaluate(const void *);
    //BaseNode * NewNode(void *, uint) const;
    T result;
    /*static void EvaluateAll(const void *, uint);
    static void SortNodes();
    static std::vector<BaseValueNode<T> *> nodes;*/
};

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

class BaseFogNode : public BaseNode{
protected:
    BaseFogNode(uint, NodeTree *);
    virtual ~BaseFogNode();
public:
    virtual void Evaluate(const void *);
    static BaseFogNode * Create(uint, NodeTree *);
    static void SortNodes();
    static void EvaluateAll(const void *, uint);
    //static std::vector<BaseFogNode *> nodes;
};

class BaseSurfaceNode : public BaseNode{
protected:
    //BaseSurfaceNode();
    BaseSurfaceNode(uint, NodeTree *);
    virtual ~BaseSurfaceNode();
public:
    virtual void Evaluate(const void *);
    static BaseSurfaceNode * Create(uint, NodeTree *);
    //BaseNode * NewNode(void *, uint) const;
    /*std::vector<dfloat3> vl;
    std::vector<duint3> tl;
    std::vector<duint4> ql;*/
    /*static void SortNodes();
    static void EvaluateAll(const void *, uint);
    static std::vector<BaseSurfaceNode *> nodes;*/
};

class BaseVectorFieldNode : public BaseNode{
protected:
	BaseVectorFieldNode(uint, NodeTree *);
	virtual ~BaseVectorFieldNode();
public:
	virtual void Evaluate(const void *);
	static BaseVectorFieldNode * Create(uint, NodeTree *);
};

/*class fBmPerlinNoise : public BaseSurfaceNode{
public:
    fBmPerlinNoise(uint);
    ~fBmPerlinNoise();
    void Evaluate();
    //static Create() declaration in scene.h?
    //TODO: global static openvdb::FloatGrid member to save memory - currently only linear surface node trees are supported anyway
};*/

class ISurfaceInput : public virtual BaseSurfaceNode{
protected:
    ISurfaceInput(uint, NodeTree *);
    ~ISurfaceInput();
public:
    virtual void Evaluate(const void *) = 0;
    //virtual void CommitSurface(std::vector<dfloat3>, std::vector<uint>) = 0;
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
		INPUT_RASTERIZATIONRES,
		INPUT_WEIGHT,
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
		INPUT_DISTANCE,
		INPUT_ITERATIONS,
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
        INPUT_SURFACE,
        INPUT_COUNT
    };
};

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
        //INPUT_QSCALE,
        INPUT_SURFACE,
        //INPUT_GRID,
        INPUT_COUNT
    };
    /*enum OUTPUT{
        OUTPUT_SURFACE,
        OUTPUT_GRID
    };*/
    /*uint octaves;
    const PropertyDef<uint> cl[1] = {
        {"octaves",&octaves},
        {0,0}
    };*/
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
        INPUT_FIELD,
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

BaseNode * CreateNodeByType(const char *, uint, NodeTree *);
BaseNode * CreateNodeBySocket(const char *, const void *, uint, NodeTree *);
//void EvaluateValueGroup(uint);
//void SortNodes();
//void DeleteNodes();

}

#endif
