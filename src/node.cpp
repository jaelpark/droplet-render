#include "main.h"
#include "node.h"
#include <Python.h>

#include <algorithm>
#include <functional>

namespace Node{

IValueNodeParams::IValueNodeParams(){
	//
}

IValueNodeParams::~IValueNodeParams(){
	//
}

BaseNode::BaseNode(uint _level, NodeTree *_pntree) : imask(0), omask(0), emask(0), level(_level), pntree(_pntree){
	//printf("BaseNode()\n");
	memset(pnodes,0,sizeof(pnodes));
}

BaseNode::~BaseNode(){
	//
}

void BaseNode::Clear(){
	//
}

//HACK: prevent double listing caused by multiple inheritance
#define MIHACK(l) for(uint i = 0; i < l.size(); ++i)\
	if(l[i] == this)\
		return;

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
	//printf("BaseValueNode(T)\n");
	MIHACK(pnt->nodes0);
	result = r; //set the default value for every thread
	pnt->nodes0.push_back(this);
}

template<class T>
BaseValueNode<T>::BaseValueNode(uint level, NodeTree *pnt) : BaseNode(level,pnt){
	//printf("BaseValueNode()\n");
	MIHACK(pnt->nodes0);
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

FloatInput::FloatInput(uint level, NodeTree *pnt) : BaseValueNode<float>(level,pnt), BaseNode(level,pnt){
	//
}

FloatInput::~FloatInput(){
	//
}

void FloatInput::Evaluate(const void *pp){
	this->result.local().value[0] = dynamic_cast<BaseValueNode<float>*>(this->pnodes[0])->BaseValueNode<float>::result.local().value[this->indices[0]];
}

ScalarMath::ScalarMath(uint level, NodeTree *pnt, char _opch) : BaseValueNode<float>(level,pnt), BaseNode(level,pnt), opch(_opch){
	//
}

ScalarMath::~ScalarMath(){
	//
}

void ScalarMath::Evaluate(const void *pp){
	float a = dynamic_cast<BaseValueNode<float>*>(this->pnodes[0])->BaseValueNode<float>::result.local().value[this->indices[0]];
	float b = dynamic_cast<BaseValueNode<float>*>(this->pnodes[1])->BaseValueNode<float>::result.local().value[this->indices[1]];
	float r;
	switch(opch){
	case '+': r = a+b; break;
	case '-': r = a-b; break;
	case '*': r = a*b; break;
	case '/': r = a/b; break;
	case 'a': r = fabs(a); break;
	case 'm': r = std::min(a,b); break;
	case 'M': r = std::max(a,b); break;
	case 'q': r = sqrtf(a); break;
	case 'p': r = powf(a,b); break;
	case '0': r = floorf(a); break;
	case '1': r = ceilf(a); break;
	case 'e': r = expf(a); break;
	case 's': r = sinf(a); break;
	case 'c': r = cosf(a); break;
	case 't': r = tanf(a); break;
	case 'S': r = asinf(a); break;
	case 'C': r = acosf(a); break;
	case 'T': r = atan2f(a,b); break;
	case 'G': r = (float)(a > b); break;
	case 'g': r = (float)(a >= b); break;
	case 'L': r = (float)(a < b); break;
	case 'l': r = (float)(a <= b); break;
	default:
		r = 0.0f;
	}
	this->result.local().value[0] = r;
}

VectorInput::VectorInput(uint level, NodeTree *pnt) : BaseValueNode<dfloat3>(level,pnt), BaseNode(level,pnt){
	//
}

VectorInput::~VectorInput(){
	//
}

void VectorInput::Evaluate(const void *pp){
	this->result.local().value[0] = dfloat3(
		dynamic_cast<BaseValueNode<float>*>(this->pnodes[0])->BaseValueNode<float>::result.local().value[this->indices[0]],
		dynamic_cast<BaseValueNode<float>*>(this->pnodes[1])->BaseValueNode<float>::result.local().value[this->indices[1]],
		dynamic_cast<BaseValueNode<float>*>(this->pnodes[2])->BaseValueNode<float>::result.local().value[this->indices[2]]);
}

VectorMath::VectorMath(uint level, NodeTree *pnt, char _opch) : BaseValueNode<dfloat3>(level,pnt), BaseNode(level,pnt), opch(_opch){
	//
}

VectorMath::~VectorMath(){
	//
}

void VectorMath::Evaluate(const void *pp){
	dfloat3 sa = dynamic_cast<BaseValueNode<dfloat3>*>(this->pnodes[0])->BaseValueNode<dfloat3>::result.local().value[this->indices[0]];
	dfloat3 sb = dynamic_cast<BaseValueNode<dfloat3>*>(this->pnodes[1])->BaseValueNode<dfloat3>::result.local().value[this->indices[1]];
	float4 a = float4::load(&sa);
	float4 b = float4::load(&sb);
	float4 r;
	switch(opch){
	case '+': r = a+b; break;
	case '-': r = a-b; break;
	case '*': r = a*b; break;
	case '/': r = a/b; break;
	case 'X': r = float4::cross(a,b); break;
	case 'n': r = float4::normalize3(a); break;
	case '|':
		r = float4::dot3(a,b);
		break;
	default:
		r = float4::zero();
	}
	this->result.local().value[0] = dfloat3(r);
}

VectorMix::VectorMix(uint level, NodeTree *pnt) : BaseValueNode<dfloat3>(level,pnt), BaseNode(level,pnt){
	//
}

VectorMix::~VectorMix(){
	//
}

void VectorMix::Evaluate(const void *pp){
	dfloat3 sa = dynamic_cast<BaseValueNode<dfloat3>*>(this->pnodes[0])->BaseValueNode<dfloat3>::result.local().value[this->indices[0]];
	dfloat3 sb = dynamic_cast<BaseValueNode<dfloat3>*>(this->pnodes[1])->BaseValueNode<dfloat3>::result.local().value[this->indices[1]];
	float4 a = float4::load(&sa);
	float4 b = float4::load(&sb);
	float t = dynamic_cast<BaseValueNode<float>*>(this->pnodes[2])->BaseValueNode<float>::result.local().value[this->indices[2]];
	//
	this->result.local().value[0] = dfloat3((1.0f-t)*a+t*b);
}

VectorXYZ::VectorXYZ(uint level, NodeTree *pnt) : BaseValueNode<float>(level,pnt), BaseNode(level,pnt){
	//
}

VectorXYZ::~VectorXYZ(){
	//
}

void VectorXYZ::Evaluate(const void *pp){
	dfloat3 a = dynamic_cast<BaseValueNode<dfloat3>*>(this->pnodes[0])->BaseValueNode<dfloat3>::result.local().value[this->indices[0]];
	this->result.local().value[0] = a.x;
	this->result.local().value[1] = a.y;
	this->result.local().value[2] = a.z;
}

IFbmNoise::IFbmNoise(uint _level, NodeTree *pnt) : BaseValueNode<float>(_level,pnt), BaseValueNode<dfloat3>(_level,pnt), BaseNode(_level,pnt){
	//
}

IFbmNoise::~IFbmNoise(){
	//
}

void IFbmNoise::Evaluate(const void *pp){
	//unused, implementation at FbmNoise::Evaluate
}

IVoronoiLayers::IVoronoiLayers(uint _level, NodeTree *pnt) : BaseValueNode<float>(_level,pnt), BaseNode(_level,pnt){
	//
}

IVoronoiLayers::~IVoronoiLayers(){
	//
}

void IVoronoiLayers::Evaluate(const void *pp){
	//unused, implementation at VoronoiLayers::Evaluate
}

BaseFogNode::BaseFogNode(uint level, NodeTree *pnt) : BaseNode(level,pnt){
	//printf("BaseFogNode()\n");
	MIHACK(pnt->nodes1);
	pnt->nodes1.push_back(this);
}

BaseFogNode::~BaseFogNode(){
	//
}

void BaseFogNode::Evaluate(const void *pp){
	//
}

void BaseFogNode::Clear(){
	//
}

BaseSurfaceNode::BaseSurfaceNode(uint level, NodeTree *pnt) : BaseNode(level,pnt){
	//printf("BaseSurfaceNode()\n");
	MIHACK(pnt->nodes1);
	pnt->nodes1.push_back(this);
}

BaseSurfaceNode::~BaseSurfaceNode(){
	//
}

void BaseSurfaceNode::Clear(){
	//
}

void BaseSurfaceNode::Evaluate(const void *pp){
	//
}

BaseVectorFieldNode::BaseVectorFieldNode(uint level, NodeTree *pnt) : BaseNode(level,pnt){
	//printf("BaseVectorFieldNode()\n");
	MIHACK(pnt->nodes1);
	pnt->nodes1.push_back(this);
}

BaseVectorFieldNode::~BaseVectorFieldNode(){
	//
}

void BaseVectorFieldNode::Evaluate(const void *pp){
	//
}

void BaseVectorFieldNode::Clear(){
	//
}

VoxelInfo::VoxelInfo(uint _level, NodeTree *pnt) : BaseValueNode<float>(_level,pnt), BaseValueNode<dfloat3>(_level,pnt), BaseNode(_level,pnt){
	//
}

VoxelInfo::~VoxelInfo(){
	//
}

void VoxelInfo::Evaluate(const void *pp){
	IValueNodeParams *pd = (IValueNodeParams*)pp;

	BaseValueResult<dfloat3> &rv = this->BaseValueNode<dfloat3>::result.local();
	rv.value[OUTPUT_VECTOR_VOXPOSW] = *pd->GetVoxPosW();
	rv.value[OUTPUT_VECTOR_CPTPOSW] = *pd->GetCptPosW();

	BaseValueResult<float> &rs = this->BaseValueNode<float>::result.local();
	rs.value[OUTPUT_FLOAT_DISTANCE] = pd->GetLocalDistance();
	rs.value[OUTPUT_FLOAT_DENSITY] = pd->GetLocalDensity();
}

AdvectionInfo::AdvectionInfo(uint _level, NodeTree *pnt) : BaseValueNode<float>(_level,pnt), BaseValueNode<dfloat3>(_level,pnt), BaseNode(_level,pnt){
	//
}

AdvectionInfo::~AdvectionInfo(){
	//
}

void AdvectionInfo::Evaluate(const void *pp){
	IValueNodeParams *pd = (IValueNodeParams*)pp;

	BaseValueResult<dfloat3> &rv = this->BaseValueNode<dfloat3>::result.local();
	rv.value[OUTPUT_VECTOR_VOXPOSW] = *pd->GetVoxPosWAdv();

	BaseValueResult<float> &rs = this->BaseValueNode<float>::result.local();
	rs.value[OUTPUT_FLOAT_ADVDISTANCE] = pd->GetAdvectionDistance();
	rs.value[OUTPUT_FLOAT_DENSITY] = pd->GetAdvectionDensity();
}

ObjectInfo::ObjectInfo(uint _level, NodeTree *pnt) : BaseValueNode<dfloat3>(_level,pnt), BaseNode(_level,pnt){
	//
}

ObjectInfo::~ObjectInfo(){
	//
}

void ObjectInfo::Evaluate(const void *pp){
	IValueNodeParams *pd = (IValueNodeParams*)pp;

	BaseValueResult<dfloat3> &rv = this->BaseValueNode<dfloat3>::result.local();
	rv.value[OUTPUT_VECTOR_LOCATION] = *pd->GetObjectPosW();
}

SceneInfo::SceneInfo(uint _level, NodeTree *pnt) : BaseValueNode<float>(_level,pnt), BaseValueNode<dfloat3>(_level,pnt), BaseNode(_level,pnt){
	//
}

SceneInfo::~SceneInfo(){
	//
}

void SceneInfo::Evaluate(const void *pp){
	IValueNodeParams *pd = (IValueNodeParams*)pp;

	BaseValueNode<dfloat3> *pnode = dynamic_cast<BaseValueNode<dfloat3>*>(pnodes[INPUT_POSITION]);
	dfloat3 dposw = pnode->locr(indices[INPUT_POSITION]);

	BaseValueResult<float> &rs = this->BaseValueNode<float>::result.local();
	bool d = pd->SampleGlobalDistance(dposw,false) > 0.0f;
	rs.value[OUTPUT_FLOAT_DISTANCE] = pd->SampleGlobalDistance(dposw,true);
	rs.value[OUTPUT_FLOAT_SURFACE] = d?0.0f:1.0f;
	rs.value[OUTPUT_FLOAT_DENSITY] = pd->SampleGlobalDensity(dposw);
	rs.value[OUTPUT_FLOAT_FINAL] = d?rs.value[OUTPUT_FLOAT_DENSITY]:1.0f;

	BaseValueResult<dfloat3> &rv = this->BaseValueNode<dfloat3>::result.local();
	rv.value[OUTPUT_VECTOR_VECTOR] = pd->SampleGlobalVector(dposw);
	rv.value[OUTPUT_VECTOR_GRADIENT] = pd->SampleGlobalGradient(dposw);
}

ISurfaceInput::ISurfaceInput(uint _level, NodeTree *pnt) : BaseSurfaceNode(_level,pnt), BaseNode(_level,pnt){
	//
}

ISurfaceInput::~ISurfaceInput(){
	//
}

IParticleInput::IParticleInput(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseNode(_level,pnt){
	//
}

IParticleInput::~IParticleInput(){
	//
}

IFieldInput::IFieldInput(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseVectorFieldNode(_level,pnt), BaseNode(_level,pnt){
	//
}

IFieldInput::~IFieldInput(){
	//
}

ISmokeCache::ISmokeCache(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseVectorFieldNode(_level,pnt), BaseNode(_level,pnt){
	//
}

ISmokeCache::~ISmokeCache(){
	//
}

IFogPostInput::IFogPostInput(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseNode(_level,pnt){
	//
}

IFogPostInput::~IFogPostInput(){
	//
}

IComposite::IComposite(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseNode(_level,pnt){
	//
}

IComposite::~IComposite(){
	//
}

ICombine::ICombine(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseNode(_level,pnt){
	//
}

ICombine::~ICombine(){
	//
}

IAdvection::IAdvection(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseNode(_level,pnt){
	//
}

IAdvection::~IAdvection(){
	//
}

ISurfaceToFog::ISurfaceToFog(uint _level, NodeTree *pnt) : BaseFogNode(_level,pnt), BaseNode(_level,pnt){
	//
}

ISurfaceToFog::~ISurfaceToFog(){
	//
}

IDisplacement::IDisplacement(uint _level, NodeTree *pnt) : BaseSurfaceNode(_level,pnt), BaseNode(_level,pnt){
	//
}

IDisplacement::~IDisplacement(){
	//
}

ICSG::ICSG(uint _level, NodeTree *pnt) : BaseSurfaceNode(_level,pnt), BaseNode(_level,pnt){
	//
}

ICSG::~ICSG(){
	//
}

OutputNode::OutputNode(NodeTree *pnt, char _opch, bool _qonly) : BaseNode(0,pnt), opch(_opch), qonly(_qonly){
	//printf("OutputNode()\n");
	pnt->nodes1.push_back(this);
}

OutputNode::~OutputNode(){
	//
}

void OutputNode::Evaluate(const void *pp){
	//never used, output node won't be listed anywhere
}

NodeTree::NodeTree(const char *pn){
	strcpy(name,pn);
	ntrees.push_back(this);
}

NodeTree::~NodeTree(){
	for(uint i = 0; i < nodes0.size(); ++i)
		delete nodes0[i];
	for(uint i = 0; i < nodes1.size(); ++i)
		delete nodes1[i];
}

void NodeTree::EvaluateNodes0(const void *pp, uint max, uint mask){
	for(uint i = 0, n = nodes0.size(); i < n && nodes0[i]->level >= max; ++i)
		if(nodes0[i]->emask & mask)
			nodes0[i]->Evaluate(pp);
}

void NodeTree::EvaluateNodes1(const void *pp, uint max, uint mask){
	for(uint i = 0, n = nodes1.size(); i < n && nodes1[i]->level >= max; ++i)
		if(nodes1[i]->emask & mask)
			nodes1[i]->Evaluate(pp);
	//clear all the grids except for the root and the first level nodes connected to it (at the back of the list)
	if(nodes1.size() > 1){
		for(uint i = 0; nodes1[i]->level > 1; ++i)
			nodes1[i]->Clear();
	}
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

BaseNode * CreateNodeByType(const char *pname, const void *pnode, uint level, NodeTree *pnt){
	if(strcmp(pname,"ClNodeScalarMath") == 0){
		PyObject *pop = PyObject_GetAttrString((PyObject*)pnode,"op");
		const char opch = PyUnicode_AsUTF8(pop)[0];
		Py_DECREF(pop);

		return new ScalarMath(level,pnt,opch);

	}else if(strcmp(pname,"ClNodeVectorMath") == 0){
		PyObject *pop = PyObject_GetAttrString((PyObject*)pnode,"op");
		const char opch = PyUnicode_AsUTF8(pop)[0];
		Py_DECREF(pop);

		return new VectorMath(level,pnt,opch);

	}else if(strcmp(pname,"ClNodeVectorMix") == 0){
		return new VectorMix(level,pnt);
	}else if(strcmp(pname,"ClNodeVectorXYZ") == 0){
		return new VectorXYZ(level,pnt);
	}else if(strcmp(pname,"ClNodeFbmNoise") == 0){
		return IFbmNoise::Create(level,pnt);
	}else if(strcmp(pname,"ClNodeVoronoiLayers") == 0){
		return IVoronoiLayers::Create(level,pnt);

	}else if(strcmp(pname,"ClNodeFloatInput") == 0){
		return new FloatInput(level,pnt);
	}else if(strcmp(pname,"ClNodeVectorInput") == 0){
		return new VectorInput(level,pnt);
	}else if(strcmp(pname,"ClNodeVoxelInfo") == 0){
		return new VoxelInfo(level,pnt);
	}else if(strcmp(pname,"ClNodeAdvectionInfo") == 0){
		return new AdvectionInfo(level,pnt);
	}else if(strcmp(pname,"ClNodeObjectInfo") == 0){
		return new ObjectInfo(level,pnt);
	}else if(strcmp(pname,"ClNodeSceneInfo") == 0){
		return new SceneInfo(level,pnt);
	}else if(strcmp(pname,"ClNodeSurfaceInput") == 0){
		return ISurfaceInput::Create(level,pnt);
	}else if(strcmp(pname,"ClNodeParticleInput") == 0){
		return IParticleInput::Create(level,pnt);
	}else if(strcmp(pname,"ClNodeFieldInput") == 0){
		return IFieldInput::Create(level,pnt);
	}else if(strcmp(pname,"ClNodeSmokeCache") == 0){
		return ISmokeCache::Create(level,pnt);
	}else if(strcmp(pname,"ClNodeFogPostInput") == 0){
		return IFogPostInput::Create(level,pnt);
	}else if(strcmp(pname,"ClNodeComposite") == 0){
		return IComposite::Create(level,pnt);

	}else if(strcmp(pname,"ClNodeCombine") == 0){
		PyObject *pop = PyObject_GetAttrString((PyObject*)pnode,"op");
		const char opch = PyUnicode_AsUTF8(pop)[0];
		Py_DECREF(pop);

		return ICombine::Create(level,pnt,opch);
	}else if(strcmp(pname,"ClNodeAdvection") == 0){
		uint flags = 0;

		PyObject *pfi = PyObject_GetAttrString((PyObject*)pnode,"break_iter");
		flags |= PyObject_IsTrue(pfi)<<IAdvection::BOOL_BREAK_ITERATION;
		Py_DECREF(pfi);

		PyObject *psl = PyObject_GetAttrString((PyObject*)pnode,"sample_local");
		flags |= PyObject_IsTrue(psl)<<IAdvection::BOOL_SAMPLE_LOCAL;
		Py_DECREF(psl);

		return IAdvection::Create(level,pnt,flags);
	}else if(strcmp(pname,"ClNodeSurfaceToFog") == 0){
		PyObject *pcoff = PyObject_GetAttrString((PyObject*)pnode,"cutoff");
		float coff = (float)PyFloat_AsDouble(pcoff);
		Py_DECREF(pcoff);

		return ISurfaceToFog::Create(level,pnt,coff);
	}else if(strcmp(pname,"ClNodeDisplacement") == 0){
		PyObject *presf = PyObject_GetAttrString((PyObject*)pnode,"resf");
		float resf = (float)PyFloat_AsDouble(presf);
		Py_DECREF(presf);

		return IDisplacement::Create(level,pnt,resf);
	}else if(strcmp(pname,"ClNodeCSG") == 0){
		PyObject *pop = PyObject_GetAttrString((PyObject*)pnode,"op");
		const char opch = PyUnicode_AsUTF8(pop)[0];
		Py_DECREF(pop);

		return ICSG::Create(level,pnt,opch);
	}else if(strcmp(pname,"ClNodeSurfaceOutput") == 0){
		PyObject *pop = PyObject_GetAttrString((PyObject*)pnode,"op");
		const char opch = PyUnicode_AsUTF8(pop)[0];
		Py_DECREF(pop);

		PyObject *pbq = PyObject_GetAttrString((PyObject*)pnode,"bq");
		bool qonly = PyObject_IsTrue(pbq);
		Py_DECREF(pbq);

		return new OutputNode(pnt,opch,qonly);
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
	else if(strcmp(pname,"ClNodeFogSocket") == 0)
		return BaseFogNode::Create(level,pnt);
	else if(strcmp(pname,"ClNodeSurfaceSocket") == 0)
		return BaseSurfaceNode::Create(level,pnt);//return new BaseSurfaceNode(level);
	else if(strcmp(pname,"ClNodeVectorFieldSocket") == 0)
		return BaseVectorFieldNode::Create(level,pnt);
	return 0;

}

template class BaseValueNode<float>;
//template<> std::vector<BaseValueNode<float> *> BaseValueNode<float>::nodes = std::vector<BaseValueNode<float> *>();

template class BaseValueNode<int>;

template class BaseValueNode<dfloat3>;

std::vector<NodeTree *> NodeTree::ntrees;

}
