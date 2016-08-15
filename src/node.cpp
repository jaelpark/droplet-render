#include "main.h"
#include "node.h"
#include <python3.5m/Python.h>

#include <algorithm>
#include <functional>

namespace Node{

BaseNode::BaseNode(uint _level, NodeTree *_pntree) : omask(0), emask(0), level(_level), pntree(_pntree){
    memset(pnodes,0,sizeof(pnodes));
}

BaseNode::~BaseNode(){
    //
}

template<class T>
BaseValueResult<T>::BaseValueResult(const T &r){
	value[0] = r;
}

template<class T>
BaseValueResult<T>::BaseValueResult(){
	//
}

template<class T>
BaseValueNode<T>::BaseValueNode(T r, uint level, NodeTree *pnt) : BaseNode(level,pnt){
	for(uint i = 0; i < pnt->nodes0.size(); ++i) //HACK: prevent double listing caused by multiple inheritance
		if(pnt->nodes0[i] == this)
			return;
	//TODO: fix: in case of multiple inheritance (several different BaseValueNode types), do not add multiple instances of the same node
    result = r; //set the default value for every thread
    pnt->nodes0.push_back(this);
}

template<class T>
BaseValueNode<T>::BaseValueNode(uint level, NodeTree *pnt) : BaseNode(level,pnt){
	for(uint i = 0; i < pnt->nodes0.size(); ++i)
		if(pnt->nodes0[i] == this)
			return;
    pnt->nodes0.push_back(this);
}

template<class T>
BaseValueNode<T>::~BaseValueNode(){
    //
}

template<class T>
void BaseValueNode<T>::Evaluate(const void *pp){
    //
}

template<class T>
AddNode<T>::AddNode(uint level, NodeTree *pnt) : BaseValueNode<T>(level,pnt), BaseNode(level,pnt){
    //
}

template<class T>
AddNode<T>::~AddNode(){
    //
}

template<class T>
void AddNode<T>::Evaluate(const void *pp){
    T a = dynamic_cast<BaseValueNode<T>*>(this->pnodes[0])->BaseValueNode<T>::result.local().value[this->indices[0]];
    T b = dynamic_cast<BaseValueNode<T>*>(this->pnodes[1])->BaseValueNode<T>::result.local().value[this->indices[1]];
    this->result.local().value[0] = a+b;
    //DebugPrintf("---AddNode<float>::Evaluate(), result = %f\n",this->result);
}

template<class T>
SubNode<T>::SubNode(uint level, NodeTree *pnt) : BaseValueNode<T>(level,pnt), BaseNode(level,pnt){
    //
}

template<class T>
SubNode<T>::~SubNode(){
    //
}

template<class T>
void SubNode<T>::Evaluate(const void *pp){
    T a = dynamic_cast<BaseValueNode<T>*>(this->pnodes[0])->BaseValueNode<T>::result.local().value[this->indices[0]];
    T b = dynamic_cast<BaseValueNode<T>*>(this->pnodes[1])->BaseValueNode<T>::result.local().value[this->indices[1]];
    this->result.local().value[0] = a-b;
}

template<class T>
MulNode<T>::MulNode(uint level, NodeTree *pnt) : BaseValueNode<T>(level,pnt), BaseNode(level,pnt){
    //
}

template<class T>
MulNode<T>::~MulNode(){
    //
}

template<class T>
void MulNode<T>::Evaluate(const void *pp){
    T a = dynamic_cast<BaseValueNode<T>*>(this->pnodes[0])->BaseValueNode<T>::result.local().value[this->indices[0]];
    T b = dynamic_cast<BaseValueNode<T>*>(this->pnodes[1])->BaseValueNode<T>::result.local().value[this->indices[1]];
    this->result.local().value[0] = a*b;
}

template<class T>
DivNode<T>::DivNode(uint level, NodeTree *pnt) : BaseValueNode<T>(level,pnt), BaseNode(level,pnt){
    //
}

template<class T>
DivNode<T>::~DivNode(){
    //
}

template<class T>
void DivNode<T>::Evaluate(const void *pp){
    T a = dynamic_cast<BaseValueNode<T>*>(this->pnodes[0])->BaseValueNode<T>::result.local().value[this->indices[0]];
    T b = dynamic_cast<BaseValueNode<T>*>(this->pnodes[1])->BaseValueNode<T>::result.local().value[this->indices[1]];
    this->result.local().value[0] = a/b;
}

template<class T>
PowNode<T>::PowNode(uint level, NodeTree *pnt) : BaseValueNode<T>(level,pnt), BaseNode(level,pnt){
    //
}

template<class T>
PowNode<T>::~PowNode(){
    //
}

template<class T>
void PowNode<T>::Evaluate(const void *pp){
    T a = dynamic_cast<BaseValueNode<T>*>(this->pnodes[0])->BaseValueNode<T>::result.local().value[this->indices[0]];
    T b = dynamic_cast<BaseValueNode<T>*>(this->pnodes[1])->BaseValueNode<T>::result.local().value[this->indices[1]];
    this->result.local().value[0] = powf(a,b);
}

template<class T>
MinNode<T>::MinNode(uint level, NodeTree *pnt) : BaseValueNode<T>(level,pnt), BaseNode(level,pnt){
    //
}

template<class T>
MinNode<T>::~MinNode(){
    //
}

template<class T>
void MinNode<T>::Evaluate(const void *pp){
	T a = dynamic_cast<BaseValueNode<T>*>(this->pnodes[0])->BaseValueNode<T>::result.local().value[this->indices[0]];
    T b = dynamic_cast<BaseValueNode<T>*>(this->pnodes[1])->BaseValueNode<T>::result.local().value[this->indices[1]];
    this->result.local().value[0] = a < b?a:b;
}

template<class T>
MaxNode<T>::MaxNode(uint level, NodeTree *pnt) : BaseValueNode<T>(level,pnt), BaseNode(level,pnt){
    //
}

template<class T>
MaxNode<T>::~MaxNode(){
    //
}

template<class T>
void MaxNode<T>::Evaluate(const void *pp){
	T a = dynamic_cast<BaseValueNode<T>*>(this->pnodes[0])->BaseValueNode<T>::result.local().value[this->indices[0]];
    T b = dynamic_cast<BaseValueNode<T>*>(this->pnodes[1])->BaseValueNode<T>::result.local().value[this->indices[1]];
    this->result.local().value[0] = a > b?a:b;
}

/*FbmNoise::FbmNoise(uint _level, NodeTree *pnt) : BaseValueNode(_level,pnt){
	//
}

SalarFbmNoise::~FbmNoise(){
	//
}

void SalarFbmNoise::Evaluate(const void *){
	result = 0.0f; //fBm::noise
}*/

IFbmNoise::IFbmNoise(uint _level, NodeTree *pnt) : BaseValueNode<float>(_level,pnt), BaseValueNode<dfloat3>(_level,pnt), BaseNode(_level,pnt){
	//
}

IFbmNoise::~IFbmNoise(){
	//
}

void IFbmNoise::Evaluate(const void *pp){
	//
}

BaseFogNode::BaseFogNode(uint level, NodeTree *pnt) : BaseNode(level,pnt){
    //DebugPrintf(">> BaseSurfaceNode(level)\n");
    //nodes.push_back(this);
    pnt->nodes1.push_back(this);
}

BaseFogNode::~BaseFogNode(){
    //
}

void BaseFogNode::Evaluate(const void *pp){
    //
}

/*void BaseFogNode::SortNodes(){
    std::sort(nodes.begin(),nodes.end(),[&](BaseFogNode *pa, BaseFogNode *pb)->bool{
        return pa->level > pb->level;
    });
}

void BaseFogNode::EvaluateAll(const void *pp, uint max){
    for(uint i = 0; i < nodes.size() && nodes[i]->level > max; ++i)
        nodes[i]->Evaluate(pp);
}*/

BaseSurfaceNode::BaseSurfaceNode(uint level, NodeTree *pnt) : BaseNode(level,pnt){
    //DebugPrintf(">> BaseSurfaceNode(level)\n");
    //nodes.push_back(this);
    pnt->nodes1.push_back(this);
}

BaseSurfaceNode::~BaseSurfaceNode(){
    //
}

void BaseSurfaceNode::Evaluate(const void *pp){
    //
}

BaseVectorFieldNode::BaseVectorFieldNode(uint level, NodeTree *pnt) : BaseNode(level,pnt){
	pnt->nodes1.push_back(this);
}

BaseVectorFieldNode::~BaseVectorFieldNode(){
	//
}

void BaseVectorFieldNode::Evaluate(const void *pp){
	//
}

VoxelInfo::VoxelInfo(uint _level, NodeTree *pnt) : BaseValueNode<float>(_level,pnt), BaseValueNode<dfloat3>(_level,pnt), BaseNode(_level,pnt){
	//
}

VoxelInfo::~VoxelInfo(){
	//
}

void VoxelInfo::Evaluate(const void *pp){
	ValueNodeParams *pd = (ValueNodeParams*)pp;

	BaseValueResult<dfloat3> &rv = this->BaseValueNode<dfloat3>::result.local();
	rv.value[OUTPUT_VECTOR_VOXPOSW] = *std::get<VNP_VOXPOSW>(*pd);
	rv.value[OUTPUT_VECTOR_CPTPOSW] = *std::get<VNP_CPTPOSW>(*pd);

	BaseValueResult<float> &rs = this->BaseValueNode<float>::result.local();
	rs.value[OUTPUT_FLOAT_DISTANCE] = std::get<VNP_DISTANCE>(*pd);
	rs.value[OUTPUT_FLOAT_DENSITY] = std::get<VNP_DENSITY>(*pd);
}

ISurfaceInput::ISurfaceInput(uint _level, NodeTree *pnt) : BaseSurfaceNode(_level,pnt){
    //
    //DebugPrintf(">> ISurfaceInput()\n");
}

ISurfaceInput::~ISurfaceInput(){
    //
}

IParticleInput::IParticleInput(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt){
	//
}

IParticleInput::~IParticleInput(){
	//
}

ISmokeCache::ISmokeCache(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt){
	//
}

ISmokeCache::~ISmokeCache(){
	//
}

IComposite::IComposite(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt){
	//
}

IComposite::~IComposite(){
	//
}

IAdvection::IAdvection(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt){
	//
}

IAdvection::~IAdvection(){
	//
}

IDisplacement::IDisplacement(uint _level, NodeTree *pnt) : BaseSurfaceNode(_level,pnt){
    //
}

IDisplacement::~IDisplacement(){
    //
}

#ifdef BLCLOUD_DEPRECATED
IfBmPerlinNoise::IfBmPerlinNoise(uint _level, NodeTree *pnt) : BaseSurfaceNode(_level,pnt){
    //
}

IfBmPerlinNoise::~IfBmPerlinNoise(){
    //
}
#endif

IVectorFieldSampler::IVectorFieldSampler(uint _level, NodeTree *pnt) : BaseValueNode<dfloat3>(_level,pnt), BaseNode(_level,pnt){
	//
}

IVectorFieldSampler::~IVectorFieldSampler(){
	//
}

OutputNode::OutputNode(NodeTree *pnt) : BaseNode(0,pnt){
    pnt->nodes1.push_back(this);
    //proot = this; //hack
}

OutputNode::~OutputNode(){
    //
}

void OutputNode::Evaluate(const void *pp){
    //never used, output node won't be listed anywhere
}

NodeTree::NodeTree(){
    ntrees.push_back(this);
}

NodeTree::~NodeTree(){
    for(uint i = 0; i < nodes0.size(); ++i)
        delete nodes0[i];
    for(uint i = 0; i < nodes1.size(); ++i)
        delete nodes1[i];
}

void NodeTree::EvaluateNodes0(const void *pp, uint max, uint mask){
    for(uint i = 0; i < nodes0.size() && nodes0[i]->level >= max; ++i)
        if(nodes0[i]->emask & mask)
            nodes0[i]->Evaluate(pp);
}

void NodeTree::EvaluateNodes1(const void *pp, uint max, uint mask){
    for(uint i = 0; i < nodes1.size() && nodes1[i]->level >= max; ++i)
        if(nodes1[i]->emask & mask)
            nodes1[i]->Evaluate(pp);
}

void NodeTree::ApplyBranchMask(){
    BaseNode *proot = GetRoot();
    proot->emask = ~0u;

    std::function<void (BaseNode *, uint)> rbmask = [&](BaseNode *pnode, uint mask)->void{
        pnode->emask |= mask;
        for(uint i = 0; i < sizeof(pnode->pnodes)/sizeof(pnode->pnodes[0]); ++i)
            if(pnode->pnodes[i])
                rbmask(pnode->pnodes[i],mask);
    };

    for(uint i = 0; i < sizeof(proot->pnodes)/sizeof(proot->pnodes[0]); ++i)
        if(proot->pnodes[i])
            rbmask(proot->pnodes[i],1<<i);
}

void NodeTree::SortNodes(){
    std::sort(nodes0.begin(),nodes0.end(),[&](BaseNode *pa, BaseNode *pb)->bool{
        return pa->level > pb->level;
    });
    std::sort(nodes1.begin(),nodes1.end(),[&](BaseNode *pa, BaseNode *pb)->bool{
        return pa->level > pb->level;
    });
}

BaseNode * NodeTree::GetRoot() const{
    return nodes1.back(); //assume already sorted
}

void NodeTree::DeleteAll(){
    for(uint i = 0; i < ntrees.size(); ++i)
        delete ntrees[i];
    ntrees.clear();
}

BaseNode * CreateNodeByType(const char *pname, uint level, NodeTree *pnt){
    if(strcmp(pname,"ClNodeFloatAdd") == 0){
        return new AddNode<float>(level,pnt);
    }else if(strcmp(pname,"ClNodeFloatSub") == 0){
        return new SubNode<float>(level,pnt);
    }else if(strcmp(pname,"ClNodeFloatMul") == 0){
        return new MulNode<float>(level,pnt);
    }else if(strcmp(pname,"ClNodeFloatDiv") == 0){
        return new DivNode<float>(level,pnt);
    }else if(strcmp(pname,"ClNodeFloatPow") == 0){
        return new PowNode<float>(level,pnt);
	}else if(strcmp(pname,"ClNodeFloatMin") == 0){
		return new MinNode<float>(level,pnt);
	}else if(strcmp(pname,"ClNodeFloatMax") == 0){
		return new MaxNode<float>(level,pnt);

	}else if(strcmp(pname,"ClNodeFbmNoise") == 0){
		return IFbmNoise::Create(level,pnt);

	}else if(strcmp(pname,"ClNodeVoxelInfo") == 0){
		return new VoxelInfo(level,pnt);
    }else if(strcmp(pname,"ClNodeSurfaceInput") == 0){
        return ISurfaceInput::Create(level,pnt);
	}else if(strcmp(pname,"ClNodeParticleInput") == 0){
		return IParticleInput::Create(level,pnt);
	}else if(strcmp(pname,"ClNodeSmokeCache") == 0){
		return ISmokeCache::Create(level,pnt);
	}else if(strcmp(pname,"ClNodeComposite") == 0){
		return IComposite::Create(level,pnt);
	}else if(strcmp(pname,"ClNodeAdvection") == 0){
		return IAdvection::Create(level,pnt);
    }else if(strcmp(pname,"ClNodeDisplacement") == 0){
        return IDisplacement::Create(level,pnt);
#ifdef BLCLOUD_DEPRECATED
    }else if(strcmp(pname,"ClNodefBmPerlinNoise") == 0){
        return IfBmPerlinNoise::Create(level,pnt);
#endif
	}else if(strcmp(pname,"ClNodeVectorFieldSampler") == 0){
		return IVectorFieldSampler::Create(level,pnt);
    }else if(strcmp(pname,"ClNodeSurfaceOutput") == 0){
        return new OutputNode(pnt);
    }
    return 0;
}

BaseNode * CreateNodeBySocket(const char *pname, const void *pvalue, uint level, NodeTree *pnt){
    if(strcmp(pname,"ClNodeFloatSocket") == 0){
        float v = PyFloat_AsDouble((PyObject*)pvalue);
        return new BaseValueNode<float>(v,level,pnt);
    }else if(strcmp(pname,"ClNodeIntSocket") == 0){
        int v = PyLong_AsLong((PyObject*)pvalue);
        return new BaseValueNode<int>(v,level,pnt);
	}else if(strcmp(pname,"ClNodeVectorSocket") == 0)
		return new BaseValueNode<dfloat3>(dfloat3(0.0f),level,pnt);
    else if(strcmp(pname,"ClNodeShaderSocket") == 0)
        return BaseSurfaceNode::Create(level,pnt);//return new BaseSurfaceNode(level);//BaseNode(level);
    //else if(strcmp(pname,"ClNodeGridSocket") == 0)
        //return BaseSurfaceNode::Create(level);//return new BaseSurfaceNode(level);
    else if(strcmp(pname,"ClNodeFogSocket") == 0)
        return BaseFogNode::Create(level,pnt);
    else if(strcmp(pname,"ClNodeSurfaceSocket") == 0)
        return BaseSurfaceNode::Create(level,pnt);//return new BaseSurfaceNode(level);
	else if(strcmp(pname,"ClNodeVectorFieldSocket") == 0)
		return BaseVectorFieldNode::Create(level,pnt);
    return 0;

}

//force compile
template class BaseValueNode<float>;
//template<> std::vector<BaseValueNode<float> *> BaseValueNode<float>::nodes = std::vector<BaseValueNode<float> *>();
template class AddNode<float>;
template class SubNode<float>;
template class MulNode<float>;
template class DivNode<float>;
template class PowNode<float>;
template class MinNode<float>;
template class MaxNode<float>;

template class BaseValueNode<int>;
//template<> std::vector<BaseValueNode<int> *> BaseValueNode<int>::nodes = std::vector<BaseValueNode<int> *>();
/*template class AddNode<int>;
template class SubNode<int>;
template class MulNode<int>;
template class DivNode<int>;
template class PowNode<int>;*/

template class BaseValueNode<dfloat3>;

std::vector<NodeTree *> NodeTree::ntrees;

}
