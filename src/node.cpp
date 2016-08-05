#include "main.h"
#include "node.h"
#include <python3.5m/Python.h>

#include <algorithm>
#include <functional>

namespace Node{

/*BaseNode::BaseNode(){
    //
}*/

BaseNode::BaseNode(uint _level, NodeTree *_pntree) : omask(0), emask(0), level(_level), pntree(_pntree){
    memset(pnodes,0,sizeof(pnodes));
}

BaseNode::~BaseNode(){
    //
}

/*void BaseNode::SetConstants(const void *_pps){
    //
}*/

template<class T>
BaseValueNode<T>::BaseValueNode(T r, uint level, NodeTree *pnt) : BaseNode(level,pnt){
    result = r;
    pnt->nodes0.push_back(this);
    //nodes.push_back(this);
}

template<class T>
BaseValueNode<T>::BaseValueNode(uint level, NodeTree *pnt) : BaseNode(level,pnt){
    //nodes.push_back(this);
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

/*template<class T>
void BaseValueNode<T>::SortNodes(){
    std::sort(nodes.begin(),nodes.end(),[&](BaseValueNode<T> *pa, BaseValueNode<T> *pb)->bool{
        return pa->level > pb->level;
    });
}

template<class T>
void BaseValueNode<T>::EvaluateAll(const void *pp, uint max){
    //TODO: value nodes of all types should be evaluated simultanously. Need per-level evaluation.
    for(uint i = 0; i < nodes.size() && nodes[i]->level > max; ++i)
        nodes[i]->Evaluate(pp);
}*/

template<class T>
AddNode<T>::AddNode(uint level, NodeTree *pnt) : BaseValueNode<T>(level,pnt){
    //
}

template<class T>
AddNode<T>::~AddNode(){
    //
}

template<class T>
void AddNode<T>::Evaluate(const void *pp){
    T a = ((BaseValueNode<T>*)this->pnodes[0])->result;
    T b = ((BaseValueNode<T>*)this->pnodes[1])->result;
    this->result = a+b;
    //DebugPrintf("---AddNode<float>::Evaluate(), result = %f\n",this->result);
}

template<class T>
SubNode<T>::SubNode(uint level, NodeTree *pnt) : BaseValueNode<T>(level,pnt){
    //
}

template<class T>
SubNode<T>::~SubNode(){
    //
}

template<class T>
void SubNode<T>::Evaluate(const void *pp){
    T a = ((BaseValueNode<T>*)this->pnodes[0])->result;
    T b = ((BaseValueNode<T>*)this->pnodes[1])->result;
    this->result = a-b;
}

template<class T>
MulNode<T>::MulNode(uint level, NodeTree *pnt) : BaseValueNode<T>(level,pnt){
    //
}

template<class T>
MulNode<T>::~MulNode(){
    //
}

template<class T>
void MulNode<T>::Evaluate(const void *pp){
    T a = ((BaseValueNode<T>*)this->pnodes[0])->result;
    T b = ((BaseValueNode<T>*)this->pnodes[1])->result;
    this->result = a*b;
}

template<class T>
DivNode<T>::DivNode(uint level, NodeTree *pnt) : BaseValueNode<T>(level,pnt){
    //
}

template<class T>
DivNode<T>::~DivNode(){
    //
}

template<class T>
void DivNode<T>::Evaluate(const void *pp){
    T a = ((BaseValueNode<T>*)this->pnodes[0])->result;
    T b = ((BaseValueNode<T>*)this->pnodes[1])->result;
    this->result = a/b;
}

template<class T>
PowNode<T>::PowNode(uint level, NodeTree *pnt) : BaseValueNode<T>(level,pnt){
    //
}

template<class T>
PowNode<T>::~PowNode(){
    //
}

template<class T>
void PowNode<T>::Evaluate(const void *pp){
    T a = ((BaseValueNode<T>*)this->pnodes[0])->result;
    T b = ((BaseValueNode<T>*)this->pnodes[1])->result;
    this->result = powf(a,b);
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

/*void BaseSurfaceNode::SortNodes(){
    std::sort(nodes.begin(),nodes.end(),[&](BaseSurfaceNode *pa, BaseSurfaceNode *pb)->bool{
        return pa->level > pb->level;
    });
}

void BaseSurfaceNode::EvaluateAll(const void *pp, uint max){
    for(uint i = 0; i < nodes.size() && nodes[i]->level > max; ++i)
        nodes[i]->Evaluate(pp);
}*/

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

IfBmPerlinNoise::IfBmPerlinNoise(uint _level, NodeTree *pnt) : BaseSurfaceNode(_level,pnt){
    //
}

IfBmPerlinNoise::~IfBmPerlinNoise(){
    //
}

/*void IfBmPerlinNoise::SetConstants(const void *_pps){
    PyObject *pps = (PyObject*)_pps;
    PyObject *poc = PyObject_GetAttrString(pps,"octaves");
    octaves = PyLong_AsLong(poc);
    printf("octaves = %u\n",octaves);
    Py_DECREF(poc);
}*/

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
    for(uint i = 0; i < nodes0.size() && nodes0[i]->level > max; ++i)
        if(nodes0[i]->emask & mask)
            nodes0[i]->Evaluate(pp);
}

void NodeTree::EvaluateNodes1(const void *pp, uint max, uint mask){
    for(uint i = 0; i < nodes1.size() && nodes1[i]->level > max; ++i)
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
        return new DivNode<float>(level,pnt);
    }else if(strcmp(pname,"ClNodeFloatMul") == 0){
        return new MulNode<float>(level,pnt);
    }else if(strcmp(pname,"ClNodeFloatDiv") == 0){
        return new DivNode<float>(level,pnt);
    }else if(strcmp(pname,"ClNodeFloatPow") == 0){
        return new PowNode<float>(level,pnt);

    }else if(strcmp(pname,"ClNodeSurfaceInput") == 0){
        return ISurfaceInput::Create(level,pnt);
	}else if(strcmp(pname,"ClNodeParticleInput") == 0){
		return IParticleInput::Create(level,pnt);
	}else if(strcmp(pname,"ClNodeAdvection") == 0){
		return IAdvection::Create(level,pnt);
    }else if(strcmp(pname,"ClNodeDisplacement") == 0){
        return IDisplacement::Create(level,pnt);
    }else if(strcmp(pname,"ClNodefBmPerlinNoise") == 0){
        return IfBmPerlinNoise::Create(level,pnt);
        //return new fBmPerlinNoise(level);
        //ClNodeSurface
        //
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
    }else if(strcmp(pname,"ClNodeShaderSocket") == 0)
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

/*void EvaluateValueGroup(uint max){
    //TODO: nodes should be mixed and grouped by their levels
    BaseValueNode<float>::EvaluateAll(0,max);
    BaseValueNode<int>::EvaluateAll(0,max);
}

void SortNodes(){
    BaseValueNode<float>::SortNodes();
    BaseSurfaceNode::SortNodes();
}

void DeleteNodes(){
    for(uint i = 0; i < BaseValueNode<float>::nodes.size(); ++i)
        delete BaseValueNode<float>::nodes[i];
    BaseValueNode<float>::nodes.clear();
    for(uint i = 0; i < BaseSurfaceNode::nodes.size(); ++i)
        delete BaseSurfaceNode::nodes[i];
    BaseSurfaceNode::nodes.clear();

    //delete OutputNode::proot;
}*/

//force compile
template class BaseValueNode<float>;
//template<> std::vector<BaseValueNode<float> *> BaseValueNode<float>::nodes = std::vector<BaseValueNode<float> *>();
template class AddNode<float>;
template class SubNode<float>;
template class MulNode<float>;
template class DivNode<float>;
template class PowNode<float>;

template class BaseValueNode<int>;
//template<> std::vector<BaseValueNode<int> *> BaseValueNode<int>::nodes = std::vector<BaseValueNode<int> *>();
template class AddNode<int>;
template class SubNode<int>;
template class MulNode<int>;
template class DivNode<int>;
template class PowNode<int>;

std::vector<NodeTree *> NodeTree::ntrees;

/*std::vector<BaseFogNode *> BaseFogNode::nodes;

std::vector<BaseSurfaceNode *> BaseSurfaceNode::nodes;*/

//OutputNode * OutputNode::proot = 0;

}
