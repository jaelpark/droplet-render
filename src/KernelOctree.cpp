#include "main.h"
#include "scene.h"
#include "kernel.h"
#include "KernelOctree.h"

namespace KernelOctree{

//http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.29.987
static uint OctreeFirstNode(const dfloat3 &t0, const dfloat3 &tm){
	uint a = 0;
	if(t0.x > t0.y){
		if(t0.x > t0.z){
			if(tm.y < t0.x)
				a |= 2;
			if(tm.z < t0.x)
				a |= 1;
			return a;
		}
	}else{
		if(t0.y > t0.z){
			if(tm.x < t0.y)
				a |= 4;
			if(tm.z < t0.y)
				a |= 1;
			return a;
		}
	}

	if(tm.x < t0.z)
		a |= 4;
	if(tm.y < t0.z)
		a |= 2;
	return a;
}

static uint OctreeNextNode(const dfloat3 &tm, const duint3 &xyz){
	if(tm.x < tm.y){
		if(tm.x < tm.z)
			return xyz.x;
	}else{
		if(tm.y < tm.z)
			return xyz.y;
	}
	return xyz.z;
}

static void OctreeInitialize(const float4 &ro, const float4 &rd, const tbb::concurrent_vector<OctreeStructure> *pob, dfloat3 &t0s, dfloat3 &t1s, uint &a){
	float4 ce = float4::load(&(*pob)[0].ce);
	float4 ro1 = ro-ce+ce.splat<3>();
	float4 rd1 = rd;
	float4 scaabbmin = float4::zero();
	float4 scaabbmax = 2.0f*ce.splat<3>();

	dfloat3 ros = dfloat3(ro1);
	dfloat3 rds = dfloat3(rd1);

	a = 0;
	if(rds.x < 0.0f){
		ros.x = 2.0f*(*pob)[0].ce.w-ros.x;
		rds.x = -rds.x;
		a |= 4;
	}
	if(rds.y < 0.0f){
		ros.y = 2.0f*(*pob)[0].ce.w-ros.y;
		rds.y = -rds.y;
		a |= 2;
	}
	if(rds.z < 0.0f){
		ros.z = 2.0f*(*pob)[0].ce.w-ros.z;
		rds.z = -rds.z;
		a |= 1;
	}

	ro1 = float4::load(&ros);
	rd1 = float4::load(&rds);

	t0s = dfloat3((scaabbmin-ro1)/rd1);
	t1s = dfloat3((scaabbmax-ro1)/rd1);
}

BaseOctreeTraverser::BaseOctreeTraverser(){
	//
}

BaseOctreeTraverser::~BaseOctreeTraverser(){
	//
}

OctreeFullTraverser::OctreeFullTraverser(){
	//
}

OctreeFullTraverser::~OctreeFullTraverser(){
	//
}

void OctreeFullTraverser::Initialize(const sfloat4 &ro, const sfloat4 &rd, const dintN &gm, const dfloatN &maxd, const tbb::concurrent_vector<OctreeStructure> *_pob){
	pob = _pob;
	for(uint i = 0, a; i < BLCLOUD_VSIZE; ++i){
		ls[i].clear();
		if(gm.v[i] == 0)
			continue;

		dfloat3 t0s, t1s;
		OctreeInitialize(ro.get(i),rd.get(i),pob,t0s,t1s,a);
		if(std::max(std::max(t0s.x,t0s.y),t0s.z) < std::min(std::min(t1s.x,t1s.y),t1s.z))
			OctreeProcessSubtree(t0s,t1s,a,0,0,maxd.v[i],&ls[i]);
	}
}

//Use recursive algorithm for full traversals, which seems to be a bit faster
bool OctreeFullTraverser::OctreeProcessSubtree(const dfloat3 &t0, const dfloat3 &t1, uint a, uint n, uint l, float maxd, std::vector<Node> *pls){
	//n: node index, l: octree level
	if(t1.x < 0.0f || t1.y < 0.0f || t1.z < 0.0f || (n == 0 && l > 0))
		return false;

	if((*pob)[n].volx[VOLUME_BUFFER_SDF] != ~0u || (*pob)[n].volx[VOLUME_BUFFER_FOG] != ~0u){
		float tr0 = std::max(std::max(t0.x,t0.y),t0.z);
		if(tr0 > maxd)
			return true;
		float tr1 = std::min(std::min(t1.x,t1.y),t1.z);
		pls->push_back(Node(n,tr0,tr1));
		return false;
	}

	dfloat3 tm = dfloat3(
		0.5f*(t0.x+t1.x),
		0.5f*(t0.y+t1.y),
		0.5f*(t0.z+t1.z));
	uint tn = OctreeFirstNode(t0,tm); //current node

	do{
		//octree node convention:
		//buildIndex(x) = {0, 4, 2, 6, 1, 5, 3, 7}[x] = buildIndex^-1(x)
		static const uint nla[] = {0,4,2,6,1,5,3,7};
		switch(tn){
		case 0:
			if(OctreeProcessSubtree(t0,tm,a,(*pob)[n].chn[nla[a]],l+1,maxd,pls))
				return true;
			tn = OctreeNextNode(tm,duint3(4,2,1));
			break;
		case 1:
			if(OctreeProcessSubtree(dfloat3(t0.x,t0.y,tm.z),dfloat3(tm.x,tm.y,t1.z),a,(*pob)[n].chn[nla[1^a]],l+1,maxd,pls))
				return true;
			tn = OctreeNextNode(dfloat3(tm.x,tm.y,t1.z),duint3(5,3,8));
			break;
		case 2:
			if(OctreeProcessSubtree(dfloat3(t0.x,tm.y,t0.z),dfloat3(tm.x,t1.y,tm.z),a,(*pob)[n].chn[nla[2^a]],l+1,maxd,pls))
				return true;
			tn = OctreeNextNode(dfloat3(tm.x,t1.y,tm.z),duint3(6,8,3));
			break;
		case 3:
			if(OctreeProcessSubtree(dfloat3(t0.x,tm.y,tm.z),dfloat3(tm.x,t1.y,t1.z),a,(*pob)[n].chn[nla[3^a]],l+1,maxd,pls))
				return true;
			tn = OctreeNextNode(dfloat3(tm.x,t1.y,t1.z),duint3(7,8,8));
			break;
		case 4:
			if(OctreeProcessSubtree(dfloat3(tm.x,t0.y,t0.z),dfloat3(t1.x,tm.y,tm.z),a,(*pob)[n].chn[nla[4^a]],l+1,maxd,pls))
				return true;
			tn = OctreeNextNode(dfloat3(t1.x,tm.y,tm.z),duint3(8,6,5));
			break;
		case 5:
			if(OctreeProcessSubtree(dfloat3(tm.x,t0.y,tm.z),dfloat3(t1.x,tm.y,t1.z),a,(*pob)[n].chn[nla[5^a]],l+1,maxd,pls))
				return true;
			tn = OctreeNextNode(dfloat3(t1.x,tm.y,t1.z),duint3(8,7,8));
			break;
		case 6:
			if(OctreeProcessSubtree(dfloat3(tm.x,tm.y,t0.z),dfloat3(t1.x,t1.y,tm.z),a,(*pob)[n].chn[nla[6^a]],l+1,maxd,pls))
				return true;
			tn = OctreeNextNode(dfloat3(t1.x,t1.y,tm.z),duint3(8,8,7));
			break;
		case 7:
			if(OctreeProcessSubtree(dfloat3(tm.x,tm.y,tm.z),dfloat3(t1.x,t1.y,t1.z),a,(*pob)[n].chn[nla[7^a]],l+1,maxd,pls))
				return true;
			tn = 8;
			break;
		}
	}while(tn < 8);

	return false;
}

dintN OctreeFullTraverser::GetLeaf(uint i, duintN *pnodes, sfloat1 &tr0, sfloat1 &tr1){
	dintN mask;
	dfloatN TR0, TR1;
	for(uint j = 0; j < BLCLOUD_VSIZE; ++j){
		if(i >= ls[j].size()){
			mask.v[j] = 0;
			continue;
		}
		pnodes->v[j] = std::get<0>(ls[j][i]);
		TR0.v[j] = std::get<1>(ls[j][i]);
		TR1.v[j] = std::get<2>(ls[j][i]);
		mask.v[j] = -1;
	}
	tr0 = sfloat1::load(&TR0);
	tr1 = sfloat1::load(&TR1);

	//return sint1::load(&mask);
	return mask;
}

OctreeStepTraverser::OctreeStepTraverser(){
	//
}

OctreeStepTraverser::~OctreeStepTraverser(){
	//
}

void OctreeStepTraverser::Initialize(const sfloat4 &ro, const sfloat4 &rd, const dintN &gm, const dfloatN &_maxd, const tbb::concurrent_vector<OctreeStructure> *_pob){
	maxd = _maxd;
	pob = _pob;
	for(uint i = 0; i < BLCLOUD_VSIZE; ++i){
		mask.v[i] = 0;
		if(gm.v[i] == 0)
			continue;

		OctreeInitialize(ro.get(i),rd.get(i),pob,stack[i].t0[0],stack[i].t1[0],stack[i].a);
		if(std::max(std::max(stack[i].t0[0].x,stack[i].t0[0].y),stack[i].t0[0].z) < std::min(std::min(stack[i].t1[0].x,stack[i].t1[0].y),stack[i].t1[0].z)){
			mask.v[i] = -1;
			stack[i].l = 0;
			stack[i].n[0] = 0;
			stack[i].p[0] = false;
		}
	}
}

dintN OctreeStepTraverser::GetLeaf(uint i, duintN *pnodes, sfloat1 &tr0, sfloat1 &tr1){
	//Only sequential traversal is supported. Index is ignored here,
	dfloatN TR0, TR1;
	for(uint j = 0; j < BLCLOUD_VSIZE; ++j){
		if(mask.v[j] == 0)
			continue;
		Node node;
		if(!OctreeProcessSubtree(stack[j].t0,stack[j].t1,stack[j].tm,stack[j].tn,stack[j].n,stack[j].p,stack[j].l,stack[j].a,maxd.v[j],&node)){
			mask.v[j] = 0;
			continue;
		}
		pnodes->v[j] = std::get<0>(node);
		TR0.v[j] = std::get<1>(node);
		TR1.v[j] = std::get<2>(node);
	}
	tr0 = sfloat1::load(&TR0);
	tr1 = sfloat1::load(&TR1);

	return mask;
}

//Use a stack-saving iterative algorithm when only partial traversal is required
bool OctreeStepTraverser::OctreeProcessSubtree(dfloat3 *pt0, dfloat3 *pt1, dfloat3 *ptm, uint *ptn, uint *pn, bool *pp, uint &l, uint a, float maxd, Node *pnode){
	for(;;){
		if(pp[l]){
			static const uint nla[] = {0,4,2,6,1,5,3,7};
			switch(ptn[l]){
#define SUBNODE(t0n,t1n,q,x)\
	pt0[l+1] = t0n;\
	pt1[l+1] = t1n;\
	ptn[l] = OctreeNextNode(pt1[l+1],q);\
	pn[l+1] = (*pob)[pn[l]].chn[nla[x^a]];\
	pp[++l] = false;
			case 0:
				pt0[l+1] = pt0[l];
				pt1[l+1] = ptm[l];
				ptn[l] = OctreeNextNode(ptm[l],duint3(4,2,1));
				pn[l+1] = (*pob)[pn[l]].chn[nla[a]];
				pp[++l] = false;
				break;
			case 1:
				SUBNODE(dfloat3(pt0[l].x,pt0[l].y,ptm[l].z),
				dfloat3(ptm[l].x,ptm[l].y,pt1[l].z),duint3(5,3,8),1);
				break;
			case 2:
				SUBNODE(dfloat3(pt0[l].x,ptm[l].y,pt0[l].z),
				dfloat3(ptm[l].x,pt1[l].y,ptm[l].z),duint3(6,8,3),2);
				break;
			case 3:
				SUBNODE(dfloat3(pt0[l].x,ptm[l].y,ptm[l].z),
				dfloat3(ptm[l].x,pt1[l].y,pt1[l].z),duint3(7,8,8),3);
				break;
			case 4:
				SUBNODE(dfloat3(ptm[l].x,pt0[l].y,pt0[l].z),
				dfloat3(pt1[l].x,ptm[l].y,ptm[l].z),duint3(8,6,5),4);
				break;
			case 5:
				SUBNODE(dfloat3(ptm[l].x,pt0[l].y,ptm[l].z),
				dfloat3(pt1[l].x,ptm[l].y,pt1[l].z),duint3(8,7,8),5);
				break;
			case 6:
				SUBNODE(dfloat3(ptm[l].x,ptm[l].y,pt0[l].z),
				dfloat3(pt1[l].x,pt1[l].y,ptm[l].z),duint3(8,8,7),6);
				break;
			case 7:
				pt0[l+1] = dfloat3(ptm[l].x,ptm[l].y,ptm[l].z);
				pt1[l+1] = dfloat3(pt1[l].x,pt1[l].y,pt1[l].z);
				ptn[l] = 8;
				pn[l+1] = (*pob)[pn[l]].chn[nla[7^a]];
				pp[++l] = false;
				break;
			default:
				if(l == 0)
					return false;
				--l;
				break;
			}
		}else{
			if(pt1[l].x < 0.0f || pt1[l].y < 0.0f || pt1[l].z < 0.0f || (pn[l] == 0 && l > 0)){
				if(l == 0)
					return false;
				--l;
				continue;
			}

			if((*pob)[pn[l]].volx[VOLUME_BUFFER_SDF] != ~0u || (*pob)[pn[l]].volx[VOLUME_BUFFER_FOG] != ~0u){
				float tr0 = std::max(std::max(pt0[l].x,pt0[l].y),pt0[l].z);
				if(tr0 > maxd)
					return false;
				float tr1 = std::min(std::min(pt1[l].x,pt1[l].y),pt1[l].z);
				*pnode = Node(pn[l],tr0,tr1);
				--l;
				return true; //continue;
			}

			ptm[l] = dfloat3(
				0.5f*(pt0[l].x+pt1[l].x),
				0.5f*(pt0[l].y+pt1[l].y),
				0.5f*(pt0[l].z+pt1[l].z));
			ptn[l] = OctreeFirstNode(pt0[l],ptm[l]); //current node

			pp[l] = true;
		}
	}
}

}
