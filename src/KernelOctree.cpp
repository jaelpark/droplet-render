#include "main.h"
#include "scene.h"
#include "kernel.h"
#include "KernelOctree.h"

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

static bool OctreeProcessSubtree(const dfloat3 &t0, const dfloat3 &t1, uint a, const tbb::concurrent_vector<OctreeStructure> *pob, uint n, uint l, float maxd, std::vector<ParallelLeafList::Node> *pls){
	//n: node index, l: octree level
	if(t1.x < 0.0f || t1.y < 0.0f || t1.z < 0.0f || (n == 0 && l > 0))
		return false;
	//if(level == mlevel){... return;} //leaf, volume index != ~0
	if((*pob)[n].volx[VOLUME_BUFFER_SDF] != ~0u || (*pob)[n].volx[VOLUME_BUFFER_FOG] != ~0u){
		float tr0 = std::max(std::max(t0.x,t0.y),t0.z);
		if(tr0 > maxd)
			return true;
		float tr1 = std::min(std::min(t1.x,t1.y),t1.z);
		pls->push_back(ParallelLeafList::Node(n,tr0,tr1));
		return false;
	}
	//dfloat3 tm = 0.5f*(t0+t1);
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
			if(OctreeProcessSubtree(t0,tm,a,pob,(*pob)[n].chn[nla[a]],l+1,maxd,pls))
				return true;
			tn = OctreeNextNode(tm,duint3(4,2,1));
			break;
		case 1:
			if(OctreeProcessSubtree(dfloat3(t0.x,t0.y,tm.z),dfloat3(tm.x,tm.y,t1.z),a,pob,(*pob)[n].chn[nla[1^a]],l+1,maxd,pls))
				return true;
			tn = OctreeNextNode(dfloat3(tm.x,tm.y,t1.z),duint3(5,3,8));
			break;
		case 2:
			if(OctreeProcessSubtree(dfloat3(t0.x,tm.y,t0.z),dfloat3(tm.x,t1.y,tm.z),a,pob,(*pob)[n].chn[nla[2^a]],l+1,maxd,pls))
				return true;
			tn = OctreeNextNode(dfloat3(tm.x,t1.y,tm.z),duint3(6,8,3));
			break;
		case 3:
			if(OctreeProcessSubtree(dfloat3(t0.x,tm.y,tm.z),dfloat3(tm.x,t1.y,t1.z),a,pob,(*pob)[n].chn[nla[3^a]],l+1,maxd,pls))
				return true;
			tn = OctreeNextNode(dfloat3(tm.x,t1.y,t1.z),duint3(7,8,8));
			break;
		case 4:
			if(OctreeProcessSubtree(dfloat3(tm.x,t0.y,t0.z),dfloat3(t1.x,tm.y,tm.z),a,pob,(*pob)[n].chn[nla[4^a]],l+1,maxd,pls))
				return true;
			tn = OctreeNextNode(dfloat3(t1.x,tm.y,tm.z),duint3(8,6,5));
			break;
		case 5:
			if(OctreeProcessSubtree(dfloat3(tm.x,t0.y,tm.z),dfloat3(t1.x,tm.y,t1.z),a,pob,(*pob)[n].chn[nla[5^a]],l+1,maxd,pls))
				return true;
			tn = OctreeNextNode(dfloat3(t1.x,tm.y,t1.z),duint3(8,7,8));
			break;
		case 6:
			if(OctreeProcessSubtree(dfloat3(tm.x,tm.y,t0.z),dfloat3(t1.x,t1.y,tm.z),a,pob,(*pob)[n].chn[nla[6^a]],l+1,maxd,pls))
				return true;
			tn = OctreeNextNode(dfloat3(t1.x,t1.y,tm.z),duint3(8,8,7));
			break;
		case 7:
			if(OctreeProcessSubtree(dfloat3(tm.x,tm.y,tm.z),dfloat3(t1.x,t1.y,t1.z),a,pob,(*pob)[n].chn[nla[7^a]],l+1,maxd,pls))
				return true;
			tn = 8;
			break;
		}
	}while(tn < 8);

	return false;
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

void OctreeFullTraverser::Initialize(const sfloat4 &ro, const sfloat4 &rd, const dintN &gm, const dfloatN &maxd, const tbb::concurrent_vector<OctreeStructure> *pob){
	for(uint i = 0; i < BLCLOUD_VSIZE; ++i){
		if(gm.v[i] == 0)
			continue;
		float4 ce = float4::load(&(*pob)[0].ce);
		float4 ro1 = ro.get(i)-ce+ce.splat<3>();
		float4 rd1 = rd.get(i);
		float4 scaabbmin = float4::zero();
		float4 scaabbmax = 2.0f*ce.splat<3>();

		dfloat3 ros = dfloat3(ro1);
		dfloat3 rds = dfloat3(rd1);

		uint a = 0;
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

		float4 t0 = (scaabbmin-ro1)/rd1;
		float4 t1 = (scaabbmax-ro1)/rd1;

		dfloat3 t0s = dfloat3(t0);
		dfloat3 t1s = dfloat3(t1);

		if(std::max(std::max(t0s.x,t0s.y),t0s.z) < std::min(std::min(t1s.x,t1s.y),t1s.z))
			OctreeProcessSubtree(t0s,t1s,a,pob,0,0,maxd.v[i],&ls[i]);
	}
}

void OctreeFullTraverser::GetLeaf(uint x, uint i, uint *pnode, float *ptr0, float *ptr1){
	if(i >= ls[x].size())
		return;
	*pnode = std::get<0>(ls[x][i]);
	*ptr0 = std::get<1>(ls[x][i]);
	*ptr1 = std::get<2>(ls[x][i]);
}

OctreeStepTraverser::OctreeStepTraverser(){
	//
}

OctreeStepTraverser::~OctreeStepTraverser(){
	//
}

void OctreeStepTraverser::Initialize(const sfloat4 &ro, const sfloat4 &rd, const dintN &gm, const dfloatN &maxd, const tbb::concurrent_vector<OctreeStructure> *pob){
	//
}

void OctreeStepTraverser::GetLeaf(uint x, uint i, uint *pnode, float *ptr0, float *ptr1){
	//
}
