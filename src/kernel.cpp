#include "main.h"
#include "scene.h"
#include "kernel.h"
#include "KernelSampler.h"

#include <random>

#include <tbb/parallel_for.h>

#define USE_SSE2
#include "sse_mathfun.h"

#include "ArHosekSkyModel.h"

//http://developer.nvidia.com/GPUGems3/gpugems3_ch37.html
//vectorized version
inline __m128i RNG_TausStep(__m128i z, __m128i s1, __m128i s2, __m128i s3, __m128i m){
	__m128i a, b;
	a = _mm_sll_epi32(z,s1);
	a = _mm_xor_si128(a,z);
	a = _mm_srl_epi32(a,s2);
	b = _mm_and_si128(z,m);
	b = _mm_sll_epi32(b,s3);
	b = _mm_xor_si128(b,a);
	return b;
}

inline __m128i RNG_LCGStep(__m128i z, __m128i a, __m128i c){
	__m128i t1 = _mm_mul_epu32(a,z);
	__m128i t2 = _mm_mul_epu32(_mm_srli_si128(a,4),_mm_srli_si128(z,4));
	__m128i az = _mm_unpacklo_epi32(_mm_shuffle_epi32(t1,_MM_SHUFFLE(0,0,2,0)),_mm_shuffle_epi32(t2,_MM_SHUFFLE(0,0,2,0)));
	__m128i rr = _mm_add_epi32(az,c);
	return rr;
}

inline __m128 _mm_ctf_epu32(const __m128i v){
	const __m128 two16 = _mm_set1_ps(0x10000);
	//Avoid double rounding by doing two exact conversions
	//of high and low 16-bit segments
	const __m128i hi = _mm_srli_epi32(v,16);
	const __m128i lo = _mm_srli_epi32(_mm_slli_epi32(v,16),16);
	const __m128 fhi = _mm_mul_ps(_mm_cvtepi32_ps(hi),two16);
	const __m128 flo = _mm_cvtepi32_ps(lo);
	//do single rounding according to current rounding mode
	return _mm_add_ps(fhi,flo);
}

inline __m128 RNG_Sample(sint4 *prs){
	__m128i x;
    x = _mm_xor_si128(prs->v[0].v,prs->v[1].v);
    x = _mm_xor_si128(x,prs->v[2].v);
    x = _mm_xor_si128(x,prs->v[3].v);
	//__m128 q = _mm_castsi128_ps(x);
	__m128 q = _mm_ctf_epu32(x);
	__m128 f = _mm_set_ps1(2.3283064365387e-10f);
	__m128 r = _mm_mul_ps(f,q);
    prs->v[0].v = RNG_TausStep(prs->v[0].v,_mm_set1_epi32(13),_mm_set1_epi32(19),_mm_set1_epi32(12),_mm_set1_epi32(4294967294u));
    prs->v[1].v = RNG_TausStep(prs->v[1].v,_mm_set1_epi32(2),_mm_set1_epi32(25),_mm_set1_epi32(4),_mm_set1_epi32(4294967288u));
    prs->v[2].v = RNG_TausStep(prs->v[2].v,_mm_set1_epi32(3),_mm_set1_epi32(11),_mm_set1_epi32(17),_mm_set1_epi32(4294967280u));
    prs->v[3].v = RNG_LCGStep(prs->v[3].v,_mm_set1_epi32(1664525u),_mm_set1_epi32(1013904223u));
	return r;
}

inline sfloat4 RNG_SampleDir(sint4 *prs){
    //__m128 tt = _mm_set1_ps(2.0f);
    //__m128 tp = _mm_set1_ps(SM_PI);
    __m128 t = 2.0f*SM_PI*RNG_Sample(prs);
    __m128 p = SM_PI*RNG_Sample(prs);
    //__m128 r = sin_ps(p);//XMVectorSin(p);
	sfloat4 rd;
    __m128 sp, cp, st, ct;
    sincos_ps(p,&sp,&cp);
    sincos_ps(t,&st,&ct);
    rd.v[0].v = _mm_mul_ps(sp,ct);
    rd.v[1].v = _mm_mul_ps(sp,st);
    rd.v[2].v = cp;
    rd.v[3].v = _mm_setzero_ps();
	return rd;
}

inline void RNG_Init(sint4 *prs){
    //use mersenne twister to initialize the vectorized rng
    static std::mt19937 mt(1000);
    static std::uniform_int_distribution<int> rr(256,std::numeric_limits<int>::max());
    static tbb::spin_mutex sm;
    for(uint i = 0; i < 4; ++i){
        uint rs[BLCLOUD_VSIZE];
        sm.lock();
        for(uint j = 0; j < BLCLOUD_VSIZE; ++j)
            rs[j] = rr(mt);
        sm.unlock();
        prs->v[i].v = _mm_set_epi32(rs[0],rs[1],rs[2],rs[3]);
    }
}

//multiply four vectors simultaneously
inline sfloat4 mul(const sfloat4 &v, const matrix44 &m){
	sfloat4 r;
    r.v[0] = v.v[0]*m.r[0].splat<0>()+v.v[1]*m.r[1].splat<0>()+v.v[2]*m.r[2].splat<0>()+v.v[3]*m.r[3].splat<0>();
    r.v[1] = v.v[0]*m.r[0].splat<1>()+v.v[1]*m.r[1].splat<1>()+v.v[2]*m.r[2].splat<1>()+v.v[3]*m.r[3].splat<1>();
    r.v[2] = v.v[0]*m.r[0].splat<2>()+v.v[1]*m.r[1].splat<2>()+v.v[2]*m.r[2].splat<2>()+v.v[3]*m.r[3].splat<2>();
    r.v[3] = v.v[0]*m.r[0].splat<3>()+v.v[1]*m.r[1].splat<3>()+v.v[2]*m.r[2].splat<3>()+v.v[3]*m.r[3].splat<3>();
	return r;
}

inline sfloat1 IntersectSphere(const sfloat4 &pp, const sfloat4 &r, const sfloat1 &sp, float sr){
    sfloat4 p = pp-sfloat4(sp);
    sfloat1 b = sfloat4::dot3(r,p);
    sfloat1 c = sfloat4::dot3(p,p)-sfloat1(sr*sr);
    sfloat1 d2 = b*b-c;
    sfloat1 rm = sfloat1::GreaterOrEqual(d2,sfloat1::zero());
    sfloat1 sm = sfloat1::And(rm,sfloat1::one());
	//if(d2 < 0.0f)
		//return false;
    sfloat1 d = sfloat1::sqrt(sm*d2);
	/*tr0 = -b-d;
	tr1 = -b+d;
	return tr1 > 0.0f;*/
	//return (-b+d) > 0.0f;
    return sm*sfloat1::And(sfloat1::Greater(-b+d,sfloat1::zero()),sfloat1::one());
}

inline sfloat1 IntersectCube(const sfloat4 &p, const sfloat4 &r, const sfloat4 &bmin, const sfloat4 &bmax, sfloat1 &tr0, sfloat1 &tr1){
	sfloat4 t1 = (bmin-p)/r;
	sfloat4 t2 = (bmax-p)/r;
    sfloat4 tmin = sfloat4::min(t1,t2);
    sfloat4 tmax = sfloat4::max(t1,t2);
    tr0 = sfloat1::max(sfloat1::max(tmin.v[0],tmin.v[1]),tmin.v[2]);
    tr1 = sfloat1::min(sfloat1::min(tmax.v[0],tmax.v[1]),tmax.v[2]);
    return sfloat1::And(sfloat1::And(sfloat1::LessOrEqual(tr0,tr1),sfloat1::Greater(tr1,sfloat1::zero())),sfloat1::one());
}

inline float SampleVoxelSpace(const float4 &p, LeafVolume *pvol, const float4 &ce){
    float4 nv = float4((float)BLCLOUD_uN);
    float4 ni = -0.5f*(ce-ce.splat<3>()-p)*(nv-float4::one())/ce.splat<3>();
    float4 nf = float4::max(float4::floor(ni),float4::zero()); //max() shouldn't be necessary?
    float4 nc = float4::min(float4::ceil(ni),nv-float4::one());
    float4 nl = ni-nf;

    //Safe-guard to clamp indices. There's a rare case where indices go out bounds.
    nc = float4::max(float4::min(nc,float4(BLCLOUD_uN-1)),float4(0));
    nf = float4::max(float4::min(nf,float4(BLCLOUD_uN-1)),float4(0));

    duint3 a = duint3((uint)nf.get<0>(),(uint)nf.get<1>(),(uint)nf.get<2>());
    duint3 b = duint3((uint)nc.get<0>(),(uint)nc.get<1>(),(uint)nc.get<2>());

    float4 ua = float4(pvol->pvol[a.z*BLCLOUD_uN*BLCLOUD_uN+a.y*BLCLOUD_uN+a.x],pvol->pvol[a.z*BLCLOUD_uN*BLCLOUD_uN+b.y*BLCLOUD_uN+a.x],
        pvol->pvol[b.z*BLCLOUD_uN*BLCLOUD_uN+a.y*BLCLOUD_uN+a.x],pvol->pvol[b.z*BLCLOUD_uN*BLCLOUD_uN+b.y*BLCLOUD_uN+a.x]); //([x0,y0,z0],[x0,y1,z0],[x0,y0,z1],[x0,y1,z1])
    float4 ub = float4(pvol->pvol[a.z*BLCLOUD_uN*BLCLOUD_uN+a.y*BLCLOUD_uN+b.x],pvol->pvol[a.z*BLCLOUD_uN*BLCLOUD_uN+b.y*BLCLOUD_uN+b.x],
        pvol->pvol[b.z*BLCLOUD_uN*BLCLOUD_uN+a.y*BLCLOUD_uN+b.x],pvol->pvol[b.z*BLCLOUD_uN*BLCLOUD_uN+b.y*BLCLOUD_uN+b.x]); //([x1,y0,z0],[x1,y1,z0],[x1,y0,z1],[x1,y1,z1])
    float4 u = float4::lerp(ua,ub,nl.splat<0>()); //([x,y0,z0],[x,y1,z0],[x,y0,z1],[x,y1,z1])

    float4 va = u.swizzle<0,2,0,2>();
    float4 vb = u.swizzle<1,3,1,3>();
    float4 v = float4::lerp(va,vb,nl.splat<1>()); //([x,y,z0],[x,y,z1],-,-)

    float4 wa = v.splat<0>();
    float4 wb = v.splat<1>();
    float4 w = float4::lerp(wa,wb,nl.splat<2>()); //([x,y,z],-,-,-)

    return w.get<0>();
}

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

static void OctreeProcessSubtree(const dfloat3 &t0, const dfloat3 &t1, uint a, const tbb::concurrent_vector<OctreeStructure> *pob, uint n, uint l, std::vector<uint> *pls){
	//n: node index, l: octree level
	if(t1.x < 0.0f || t1.y < 0.0f || t1.z < 0.0f || (n == 0 && l > 0))
		return;
	//if(level == mlevel){... return;} //leaf, volume index != ~0
    if((*pob)[n].volx[VOLUME_BUFFER_SDF] != ~0u || (*pob)[n].volx[VOLUME_BUFFER_FOG] != ~0u){
		pls->push_back(n);
        return; //true: stop
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
		const uint nla[] = {0,4,2,6,1,5,3,7};
		switch(tn){
		case 0:
			OctreeProcessSubtree(t0,tm,a,pob,(*pob)[n].chn[nla[a]],l+1,pls);
            tn = OctreeNextNode(tm,duint3(4,2,1));
			break;
		case 1:
            OctreeProcessSubtree(dfloat3(t0.x,t0.y,tm.z),dfloat3(tm.x,tm.y,t1.z),a,pob,(*pob)[n].chn[nla[1^a]],l+1,pls);
            tn = OctreeNextNode(dfloat3(tm.x,tm.y,t1.z),duint3(5,3,8));
			break;
		case 2:
            OctreeProcessSubtree(dfloat3(t0.x,tm.y,t0.z),dfloat3(tm.x,t1.y,tm.z),a,pob,(*pob)[n].chn[nla[2^a]],l+1,pls);
            tn = OctreeNextNode(dfloat3(tm.x,t1.y,tm.z),duint3(6,8,3));
			break;
		case 3:
            OctreeProcessSubtree(dfloat3(t0.x,tm.y,tm.z),dfloat3(tm.x,t1.y,t1.z),a,pob,(*pob)[n].chn[nla[3^a]],l+1,pls);
            tn = OctreeNextNode(dfloat3(tm.x,t1.y,t1.z),duint3(7,8,8));
			break;
		case 4:
            OctreeProcessSubtree(dfloat3(tm.x,t0.y,t0.z),dfloat3(t1.x,tm.y,tm.z),a,pob,(*pob)[n].chn[nla[4^a]],l+1,pls);
            tn = OctreeNextNode(dfloat3(t1.x,tm.y,tm.z),duint3(8,6,5));
			break;
		case 5:
            OctreeProcessSubtree(dfloat3(tm.x,t0.y,tm.z),dfloat3(t1.x,tm.y,t1.z),a,pob,(*pob)[n].chn[nla[5^a]],l+1,pls);
            tn = OctreeNextNode(dfloat3(t1.x,tm.y,t1.z),duint3(8,7,8));
			break;
		case 6:
            OctreeProcessSubtree(dfloat3(tm.x,tm.y,t0.z),dfloat3(t1.x,t1.y,tm.z),a,pob,(*pob)[n].chn[nla[6^a]],l+1,pls);
            tn = OctreeNextNode(dfloat3(t1.x,t1.y,tm.z),duint3(8,8,7));
			break;
		case 7:
            OctreeProcessSubtree(dfloat3(tm.x,tm.y,tm.z),dfloat3(t1.x,t1.y,t1.z),a,pob,(*pob)[n].chn[nla[7^a]],l+1,pls);
			tn = 8;
			break;
		}
	}while(tn < 8);
}

inline void SamplingBasis(const sfloat4 &iv, sfloat4 *pb1, sfloat4 *pb2){
    //TODO: should probably handle zero-cases
    pb1->v[0] = -iv.v[2];
    pb1->v[1] = sfloat1::zero();
    pb1->v[2] = iv.v[0];
    pb1->v[3] = sfloat1::zero();
    *pb1 /= sfloat4::length3(iv.swizzle<0,3,2,3>()); //assume iv.w to be zero
    *pb2 = sfloat4::cross3(iv,*pb1);
}

/*#define PHASE_G 0.35f//0.88f //0.32f
inline sfloat1 HG_Phase(const sfloat1 &ct){
    sfloat1 g = sfloat1(PHASE_G);
    sfloat1 e = sfloat1(1.50f);
    return (1.0f-g*g)/(4.0f*SM_PI*sfloat1::pow(1.0f+g*g-2.0f*g*ct,e));
}

inline sfloat4 HG_Sample(const sfloat4 &iv, sint4 *prs){
    sfloat1 g = sfloat1(PHASE_G);
    sfloat1 sq = (1.0f-g*g)/(1.0f-g+2.0f*g*RNG_Sample(prs));
    sfloat1 ct = (1.0f+g*g-sq*sq)/(2.0f*g);
    sfloat1 st = sfloat1::sqrt(sfloat1::max(1.0f-ct*ct,sfloat1::zero()));
    sfloat1 ph = 2.0f*SM_PI*RNG_Sample(prs);

    sfloat4 b1, b2;
    SamplingBasis(iv,&b1,&b2);
    sfloat1 sph, cph;
    sincos_ps(ph.v,&sph.v,&cph.v);
    return b1*st*cph+b2*st*sph+iv*ct;
}*/

inline sfloat1 L_Pdf(const sfloat4 &iv, const sfloat1 &la){
    sfloat1 ctm = sfloat1::sqrt(1.0f-la*la);
    return sfloat1(1.0f/(2.0f*SM_PI*(1.0f-ctm)));
}

inline sfloat4 L_Sample(const sfloat4 &iv, const sfloat1 &la, sint4 *prs){
    //iv: light vector in this case
    //t = asin(theta = r/d)
    //cos(asin(t)) = sqrt(1-(r/d)^2)
    sfloat1 u1 = RNG_Sample(prs);
    sfloat1 ctm = sfloat1::sqrt(1.0f-la*la);
    sfloat1 ct = (1.0f-u1)+u1*ctm;
    sfloat1 st = sfloat1::sqrt(1.0f-ct*ct);
    sfloat1 ph = 2.0f*SM_PI*RNG_Sample(prs);

    sfloat4 b1, b2;
    SamplingBasis(iv,&b1,&b2);
    sfloat1 sph, cph;
    sincos_ps(ph.v,&sph.v,&cph.v);
    return b1*st*cph+b2*st*sph+iv*ct;
}

static sfloat4 SampleVolume(sfloat4 ro, sfloat4 rd, sfloat1 gm, RenderKernel *pkernel, sint4 *prs, ParallelLeafList &ls, uint r, uint samples, sfloat1 *prq){
	dintN sgm = dintN(gm);
	for(uint i = 0; i < BLCLOUD_VSIZE; ++i){
		if(sgm.v[i] == 0)
			continue;
        float4 ce = float4::load(&pkernel->pscene->ob[0].ce);
        float4 ro1 = ro.get(i)-ce+ce.splat<3>();
        float4 rd1 = rd.get(i);
        float4 scaabbmin = float4::zero();
        float4 scaabbmax = 2.0f*ce.splat<3>();

        dfloat3 ros = dfloat3(ro1);
        dfloat3 rds = dfloat3(rd1);

		uint a = 0;
		if(rds.x < 0.0f){
			ros.x = 2.0f*pkernel->pscene->ob[0].ce.w-ros.x;
			rds.x = -rds.x;
			a |= 4;
		}
		if(rds.y < 0.0f){
			ros.y = 2.0f*pkernel->pscene->ob[0].ce.w-ros.y;
			rds.y = -rds.y;
			a |= 2;
		}
		if(rds.z < 0.0f){
			ros.z = 2.0f*pkernel->pscene->ob[0].ce.w-ros.z;
			rds.z = -rds.z;
			a |= 1;
		}

        ro1 = float4::load(&ros);
        rd1 = float4::load(&rds);

		//float3 t0 = (make_float3(*(float4*)&ob[0].ce)-make_float3(ob[0].ce.w)-ro)*invrd;
		//float3 t1 = (make_float3(*(float4*)&ob[0].ce)+make_float3(ob[0].ce.w)-ro)*invrd;
        float4 t0 = (scaabbmin-ro1)/rd1;
        float4 t1 = (scaabbmax-ro1)/rd1;

        dfloat3 t0s = dfloat3(t0);
        dfloat3 t1s = dfloat3(t1);

        if(std::max(std::max(t0s.x,t0s.y),t0s.z) < std::min(std::min(t1s.x,t1s.y),t1s.z))
			OctreeProcessSubtree(t0s,t1s,a,&pkernel->pscene->ob,0,0,&ls.ls[r][i]);
	}

    sfloat4 c = sfloat1::zero();

	//rc-tr1[0],tr1[0]-tr1[1],...
	/*rc = ro;
	for(ll){
		if(tr0[i] > 0.0f){
			//if outside the next leaf...
			if(volx[SDF] != ~0u){
				d0 = distance(tr0[i]); //...get the distance at the leaf intersection...
				if(d0 > 0.0f)
					rc = tr0[i]; //...and skip the empty space if outside the surface (approaching the volume from outside)
			}else rc = tr0[i]; //sole fog leaf means we're completely outside surface (since the inside nodes were removed)
		}
		//if volx[SDF] != ~0u, use the max density 1.0 for the woodcock tracking
		for(woodcock){
	        rc += r*rd; //sample random position
			if(rc > tr1[i]){
	            //region after the current leaf (the region was skipped, move on to next one)
	            rc = tr1[i];
				break;
			}
	        if(rc < tr0[i]){
	            //region before the leaf (inside the volume, no empty space was skipped)
	            rho = 1; //uniform density
	            rm = 0; //ready: sample this location
				break;
	        }else{
	            //region inside the leaf
	            rho = woodcockSample(); //potentially heterogenous medium
				if(accept){
					rm = 0;
					break;
				}
			}
		}
	}

    //finally, sample with given density and position
	sample();*/

	//approximate very high orders of scattering by lowering the cross section as r gets higher
    /*sfloat1 msigmaa = sfloat1(2.3f)*expf(-2.0f*(float)r)+sfloat1(0.02f);
    //sfloat1 msigmas = sfloat1(8.0f)*expf(-2.0f*(float)r)+sfloat1(2.9f);//sfloat1(8.0f)*expf(-2.0f*(float)r)+sfloat1(2.9f);
	//0.5, 1.2, 1.9
	sfloat1 msigmas = sfloat1(8.0f)*expf(-2.0f*(float)r)+sfloat1(1.9f);//sfloat1(8.0f)*expf(-2.0f*(float)r)+sfloat1(2.9f);*/
#if 1
	sfloat1 msigmaa = sfloat1(pkernel->msigmaa);//sfloat1(0.05f);
	sfloat1 msigmas = sfloat1(pkernel->msigmas);//sfloat1(7.11f);
#else
	sfloat1 msigmaa = sfloat1(pkernel->msigmaa)*expf(-2.0f*(float)r)+sfloat1(0.02f);
	sfloat1 msigmas = sfloat1(pkernel->msigmas)*expf(-2.0f*(float)r)+sfloat1(2.9f);
#endif
    sfloat1 msigmae = msigmaa+msigmas;

    for(uint s = 0; s < samples; ++s){
        sfloat1 qm = gm;
        sfloat1 rm = sint1::trueI();
        sfloat1 zr = sint1::falseI();

		sfloat4 rc = ro;
		for(uint i = 0;; ++i){
			dintN leafcount;
			for(uint j = 0; j < BLCLOUD_VSIZE; ++j)
				leafcount.v[j] = ls.GetLeafCount(r,j);

            qm = sfloat1::And(qm,rm);
			qm = sfloat1::And(qm,sint1::Less(sint1(i),sint1::load(&leafcount)));
			if(qm.AllFalse())
				break;

            sfloat4 ce = sfloat1::zero();
			for(uint j = 0; j < BLCLOUD_VSIZE; ++j)
				if(i < ls.GetLeafCount(r,j))
                    ce.set(j,float4::load(&pkernel->pscene->ob[ls.GetLeaf(r,j,i)].ce));

			sfloat4 lo = rc; //local origin

            sfloat1 tr0, tr1;
            sfloat4 ee = ce.swizzle<3,3,3,3>(); //splat w (extent)
            IntersectCube(lo,rd,ce-ee,ce+ee,tr0,tr1);

            sfloat4 r0 = lo+rd*tr0;
			sfloat4 r1 = lo+rd*tr1;

            //sint1 sm = sfloat1::And(qm,sfloat1::Greater(tr0,zr));
			sint1 sm = sfloat1::Greater(tr0,zr);
			dintN SM = dintN(sm);

			dfloatN smax1, dist1, rho1;

			dintN QM = dintN(qm);
			dintN VM; //true: next leaf is a sole fog; false: sdf exists, but fog may not
			for(uint j = 0; j < BLCLOUD_VSIZE; ++j){
				if(QM.v[j] != 0){
					if(pkernel->pscene->ob[ls.GetLeaf(r,j,i)].volx[VOLUME_BUFFER_SDF] != ~0u){
						smax1.v[j] = 1.0f;
						VM.v[j] = 0;
					}else{
						smax1.v[j] = pkernel->pscene->ob[ls.GetLeaf(r,j,i)].qval[VOLUME_BUFFER_FOG];
						VM.v[j] = 1;
					}

					if(SM.v[j] != 0 && VM.v[j] == 0)
						dist1.v[j] = SampleVoxelSpace(r0.get(j),&pkernel->pscene->pbuf[VOLUME_BUFFER_SDF][pkernel->pscene->ob[ls.GetLeaf(r,j,i)].volx[VOLUME_BUFFER_SDF]],ce.get(j));
					else dist1.v[j] = 1.0f;
				}else VM.v[j] = 0;
			}
			sfloat1 d0 = sfloat1::load(&dist1);
			sint1 vm = sint1::load(&VM);

            sm = sfloat1::And(sm,sfloat1::Greater(d0,zr));
			sm = sfloat1::Or(sm,vm); //skip if the next leaf is a sole fog

			sm = sfloat1::And(sm,qm); //can't allow any further changes in 'rc' if qm == 0

            rc.v[0] = sfloat1::Or(sfloat1::And(sm,r0.v[0]),sfloat1::AndNot(sm,rc.v[0]));
            rc.v[1] = sfloat1::Or(sfloat1::And(sm,r0.v[1]),sfloat1::AndNot(sm,rc.v[1]));
            rc.v[2] = sfloat1::Or(sfloat1::And(sm,r0.v[2]),sfloat1::AndNot(sm,rc.v[2]));

            sfloat1 s0 = sfloat1::Or(sfloat1::And(sm,tr0),sfloat1::AndNot(sm,zr)); //positive distance skipped (moving to next leaf)

			sm = qm;

            sfloat1 smax = sfloat1::load(&smax1);//sfloat1(1.0f); //local max in this leaf
            for(sfloat1 sr = -log_ps(RNG_Sample(prs))/(msigmae*smax), sc, sh;; sr -= log_ps(RNG_Sample(prs))/(msigmae*smax)){
                sm = sfloat1::And(sm,rm);
                //if(sfloat1::AllTrue(sfloat1::EqualR(sm,zr)))
				if(sfloat1(sm).AllFalse())
					break;
				//rc = lo+rd*sr;
                sc = s0+sr;

                rc.v[0] = sfloat1::Or(sfloat1::And(sm,lo.v[0]+rd.v[0]*sc),sfloat1::AndNot(sm,rc.v[0]));
                rc.v[1] = sfloat1::Or(sfloat1::And(sm,lo.v[1]+rd.v[1]*sc),sfloat1::AndNot(sm,rc.v[1]));
                rc.v[2] = sfloat1::Or(sfloat1::And(sm,lo.v[2]+rd.v[2]*sc),sfloat1::AndNot(sm,rc.v[2]));

                sm = sfloat1::And(sm,sfloat1::Less(sc,tr1)); //check if out of extents
                sh = sfloat1::Or(sm,sfloat1::AndNot(rm,sint1::trueI())); //prevent modifications if rm == false

                rc.v[0] = sfloat1::Or(sfloat1::AndNot(sh,r1.v[0]),sfloat1::And(sh,rc.v[0])); //rc = tr1[i]
                rc.v[1] = sfloat1::Or(sfloat1::AndNot(sh,r1.v[1]),sfloat1::And(sh,rc.v[1]));
                rc.v[2] = sfloat1::Or(sfloat1::AndNot(sh,r1.v[2]),sfloat1::And(sh,rc.v[2]));

				//if sm == false, this is always true (so rm won't be changed)
				//-> unless the ray has run out of leafs
                rm = sfloat1::Or(sfloat1::And(sm,sfloat1::Greater(sc,tr0)),sfloat1::AndNot(sm,rm));

                /*sh = sfloat1::And(sm,rm);
				dintN SH = dintN(sh);
				for(uint j = 0; j < BLCLOUD_VSIZE; ++j)
					dist1.v[j] = SH.v[j] != 0?SampleVoxelSpace(rc.get(j),&pvol[ob[ls.GetLeaf(r,j,i)].volx[VOLUME_BUFFER_SDF]],ce.get(j)):1.0f;
				sfloat1 d = sfloat1::load(&dist1);

				rm = sfloat1::And(rm,sfloat1::Greater(d,zr));*/

				sh = sfloat1::And(sm,rm);
				dintN SH = dintN(sh);
				for(uint j = 0; j < BLCLOUD_VSIZE; ++j){
					if(SH.v[j] != 0){
						if(VM.v[j] == 0){
							dist1.v[j] = SampleVoxelSpace(rc.get(j),&pkernel->pscene->pbuf[VOLUME_BUFFER_SDF][pkernel->pscene->ob[ls.GetLeaf(r,j,i)].volx[VOLUME_BUFFER_SDF]],ce.get(j));
							//have to check if fog exists
							if(dist1.v[j] > 0.0f && pkernel->pscene->ob[ls.GetLeaf(r,j,i)].volx[VOLUME_BUFFER_FOG] != ~0u)
								rho1.v[j] = SampleVoxelSpace(rc.get(j),&pkernel->pscene->pbuf[VOLUME_BUFFER_FOG][pkernel->pscene->ob[ls.GetLeaf(r,j,i)].volx[VOLUME_BUFFER_FOG]],ce.get(j));
							else rho1.v[j] = -1.0f;
						}else{
							dist1.v[j] = 1.0f;
							rho1.v[j] = SampleVoxelSpace(rc.get(j),&pkernel->pscene->pbuf[VOLUME_BUFFER_FOG][pkernel->pscene->ob[ls.GetLeaf(r,j,i)].volx[VOLUME_BUFFER_FOG]],ce.get(j));
						}
					}else{
						dist1.v[j] = 1.0f;
						rho1.v[j] = -1.0f;
					}
				}
				sfloat1 d = sfloat1::load(&dist1);
				sfloat1 p = sfloat1::load(&rho1);
				sfloat1 q = RNG_Sample(prs);

				//if(p > q || d < 0) break;
				rm = sfloat1::And(rm,sfloat1::And(sfloat1::Less(p/smax,q),sfloat1::Greater(d,zr)));

                /*sh = sfloat1::And(sm,rm);
                sfloat1 p = sfloat1(
                    sh.get<0>() != 0.0f?SampleVoxelSpace<1>(rc.get(0),&pvol[ob[ls.GetLeaf(r,0,i)].volx],ce.get(0)):-1.0f,
                    sh.get<1>() != 0.0f?SampleVoxelSpace<1>(rc.get(1),&pvol[ob[ls.GetLeaf(r,1,i)].volx],ce.get(1)):-1.0f,
                    sh.get<2>() != 0.0f?SampleVoxelSpace<1>(rc.get(2),&pvol[ob[ls.GetLeaf(r,2,i)].volx],ce.get(2)):-1.0f,
                    sh.get<3>() != 0.0f?SampleVoxelSpace<1>(rc.get(3),&pvol[ob[ls.GetLeaf(r,3,i)].volx],ce.get(3)):-1.0f);
                sfloat1 q = RNG_Sample(prs);
                rm = sfloat1::And(rm,sfloat1::Less(p,q));*/
			}
		}

		//sample E(rc)/T here
		//T should be the value at rc; T(rc)

        sfloat4 lc = sfloat1::zero();
		sfloat4 ll;
		//if(!sfloat1::AnyTrue(sfloat1::EqualR(rm,zr))){ //skip (sky)lighting calculations if all incident rays scatter
		if(rm.AnyTrue()){ //skip (sky)lighting calculations if all incident rays scatter (at least one of rm != 0)
	        for(uint i = 0; i < pkernel->lightc; ++i){
	            //sfloat1 lt = sfloat1::Greater(sfloat4::dot3(rd,sfloat4(float4::load(&pkernel->plights[i].direction))),sfloat1(0.96f));
	            sfloat1 lt = sfloat1::Greater(sfloat4::dot3(rd,sfloat4(float4::load(&pkernel->plights[i].direction))),sfloat1(pkernel->plights[i].angle));
	            //lm = sfloat1::Or(lm,lt); //0.995f

				//!!!!!!!!! bug - not additive
				//replace with lc += BaseLight::Evaluate
	            float4 c = float4::load(&pkernel->plights[i].color);
	            lc.v[0] = sfloat1::Or(sfloat1::And(lt,c.splat<0>()),sfloat1::AndNot(lt,lc.v[0]));
	            lc.v[1] = sfloat1::Or(sfloat1::And(lt,c.splat<1>()),sfloat1::AndNot(lt,lc.v[1]));
	            lc.v[2] = sfloat1::Or(sfloat1::And(lt,c.splat<2>()),sfloat1::AndNot(lt,lc.v[2]));
	        }

	        //skylighting
	        //sfloat1 th = sfloat1::acos(rd.v[2]); //add/remove abs to remove/get ground
	        //sfloat1 ph = sfloat1::atan2(rd.v[0],rd.v[1]);
	        //sfloat1 sth, cth, szr, czr;
	        //sincos_ps(th,&sth.v,&cth.v);
	        sfloat1 rdz = sfloat1::abs(rd.v[2]);
	        sfloat1 sth = sfloat1::sqrt(1.0f-rd.v[2]*rd.v[2]);
	        sfloat1 cth = rdz;
	        sfloat1 lth = sfloat1::zero(); //sun theta
	        sfloat1 lph = sfloat1::zero();
	        sfloat1 slth, clth, slph, clph;
	        sincos_ps(lth,&slth.v,&clth.v);
	        sincos_ps(lph,&slph.v,&clph.v);

	        //sfloat1 cph = sth*slth*cos_ps(lph-ph)+cth*clth;
	        sfloat1 rdd = rd.v[0]/rd.v[1];
	        sfloat1 acr = 1.0f/sfloat1::sqrt(rdd*rdd+1.0f); //cos(ph = atan(rdd))
	        sfloat1 czp = acr*(rdd*slph+clph); //cos(zr-ph)
	        sfloat1 cph = sth*slth*czp+cth*clth;

	        //sfloat1 cc1 = sfloat1::Greater(cph,sfloat1::one());
	        //sfloat1 cc2 = sfloat1::Less(cph,-sfloat1::one());
	        //sfloat1 ga = sfloat1::acos(cph);
	        /*ga = sfloat1::Or(sfloat1::And(cc1,zr),sfloat1::AndNot(cc1,ga));
	        ga = sfloat1::Or(sfloat1::And(cc2,sfloat1(SM_PI)),sfloat1::AndNot(cc2,ga));*/

	        //th = sfloat1::min(th,sfloat1(XM_PIDIV2-0.001f));

	        //sfloat1 gacs = cos_ps(ga);
	        //sfloat1 thcs = cos_ps(th);
	        //sfloat1 gmma = 1.570796808f+(-0.7107710905f+(-0.9654782056e-5f+
	             //(-2.721495228f+(0.2766071913e-4f+(5.591013086f+(-0.1860396098e-4f-3.690226730f*cph)*cph)*cph)*cph)*cph)*cph)*cph; //acos minimax
			sfloat1 gmma = sfloat1::acos(cph);
	        sfloat1 gacs = sfloat1::saturate(cph);
	        sfloat1 thcs = sfloat1::max(rdz,sfloat1::zero());
	        sfloat1 raym = gacs*gacs;
			sfloat4 ca;
			for(uint i = 0; i < 3; ++i){
	            sfloat1 caf[9];
				for(uint j = 0; j < 9; ++j)
	                caf[j] = sfloat1(pkernel->pskyms->configs[i][j]);
	            //arhosek_tristim_skymodel_radiance(pkernel->pskyms,sth,sga,0)
	            sfloat1 expm = exp_ps(caf[4]*gmma);
	            sfloat1 miem = (sfloat1::one()+raym)/sfloat1::pow(sfloat1::one()+caf[8]*caf[8]-sfloat1(2.0f)*caf[8]*gacs,sfloat1(1.5f));
	            sfloat1 zenh = sfloat1::sqrt(thcs);
	            ca.v[i] = (sfloat1::one()+caf[0]*exp_ps(caf[1]/(thcs+sfloat1(0.01f))))*(caf[2]+caf[3]*expm+caf[5]*raym+caf[6]*miem+caf[7]*zenh);
				//ca.v[i] *= 0.028f*pkernel->pskyms->radiances[i];
				ca.v[i] *= pkernel->pskyms->radiances[i];
				ca.v[i] = 0.00035f*sfloat1::pow(ca.v[i],2.2f); //convert to linear and adjust exposure
				//TODO: use CIE version?
			}
	        ca.v[3] = sfloat1::zero();

	        ll = lc+ca;
		}

        //if(r < pkernel->scattevs && sfloat1::AnyTrue(sfloat1::EqualR(rm,zr))){
		if(r < pkernel->scattevs && rm.AnyFalse()){
#define BLCLOUD_MULTIPLE_IMPORTANCE
#ifdef BLCLOUD_MULTIPLE_IMPORTANCE
            //sample the phase function and lights
			sfloat1 la = sfloat1(pkernel->plights[0].angle); //TODO: choose randomly one the lights. Multiply the final estimate (s2) with the total number of lights (ref776).
			sfloat1 u1 = RNG_Sample(prs), u2 = RNG_Sample(prs);
            sfloat4 srd = KernelSampler::HGPhase::ghg.Sample(rd,u1,u2);//HG_Sample(rd,prs);
            sfloat4 lrd = L_Sample(sfloat4(float4::load(&pkernel->plights[0].direction)),la,prs);

            //pdfs for the balance heuristic w_x = p_x/sum(p_i,i=0..N)
            sfloat1 p1 = KernelSampler::HGPhase::ghg.Evaluate(sfloat4::dot3(srd,rd));//HG_Phase(sfloat4::dot3(srd,rd));
            sfloat1 p2 = L_Pdf(lrd,la);

            //need two samples - with the phase sampling keep on the recursion while for the light do only single scattering
            sfloat1 gm1 = sfloat1::AndNot(rm,sint1::trueI());
			sfloat1 rq;
            //sfloat4 s1 = SampleVolume(rc,srd,gm1,pkernel,prs,ls,r+1,1,&rq); //p1 canceled by phase=pdf=p1
            //sfloat4 s2 = SampleVolume(rc,lrd,gm1,pkernel,prs,ls,BLCLOUD_MAX_RECURSION-1,1,&rq)*HG_Phase(sfloat4::dot3(lrd,rd)); //p2 canceled by the MIS estimator

            //estimator S(1)*f1*w1/p1+S(2)*f2*w2/p2 /= woodcock pdf

			//TODO: should support additional lights
			sfloat4 s1 = SampleVolume(rc,srd,gm1,pkernel,prs,ls,r+1,1,&rq);
			sfloat4 s2 = SampleVolume(rc,lrd,gm1,pkernel,prs,ls,BLCLOUD_MAX_RECURSION-1,1,&rq);

			//(HG_Phase(X)*SampleVolume(X)*p1/(p1+L_Pdf(X)))/p1 => (HG_Phase(X)=p1)*SampleVolume(X)/(p1+L_Pdf(X)) = p1*SampleVolume(X)/(p1+L_Pdf(X))
			//(HG_Phase(Y)*SampleVolume(Y)*p2/(HG_Phase(Y)+p2))/p2 => HG_phase(Y)*SampleVolume(Y)/(HG_Phase(Y)+p2)
			sfloat1 p3 = KernelSampler::HGPhase::ghg.Evaluate(sfloat4::dot3(lrd,rd));//HG_Phase(sfloat4::dot3(lrd,rd));
			sfloat4 cm = s1*p1/(p1+L_Pdf(srd,la))+s2*p3/(p3+p2);
			cm *= msigmas/msigmae;

			//s1=HG_Phase(X)*SampleVolume(X)
			//s1*p1/(p1+L_Pdf(X))
#else
			sfloat1 rq;
            //phase function sampling
            sfloat4 srd = HG_Sample(rd,prs);
            sfloat4 cm = SampleVolume(rc,srd,sfloat1::AndNot(rm,sint1::trueI()),pkernel,prs,ls,r+1,1,&rq)*msigmas/msigmae;
            //phase/pdf(=phase)=1

            //light sampling, obviously won't alone result in any sky lighting
            /*sfloat1 la = sfloat1(pkernel->plights[0].angle);
            sfloat4 lrd = L_Sample(sfloat1(float4::load(&pkernel->plights[0].direction)),la,prs);
            sfloat4 cm = SampleVolume(rc,lrd,sfloat1::AndNot(rm,sint1::trueI()),pkernel,prs,ls,r+1,1)*HG_Phase(sfloat4::dot3(lrd,rd))*msigmas
                /(msigmae*L_Pdf(lrd,la));*/
#endif
			*prq = sint1::trueI();
            c.v[0] += sfloat1::Or(sfloat1::And(rm,ll.v[0]),sfloat1::AndNot(rm,cm.v[0]));
            c.v[1] += sfloat1::Or(sfloat1::And(rm,ll.v[1]),sfloat1::AndNot(rm,cm.v[1]));
            c.v[2] += sfloat1::Or(sfloat1::And(rm,ll.v[2]),sfloat1::AndNot(rm,cm.v[2]));
            c.v[3] += sfloat1::one();
        }else if(!(pkernel->flags & RENDER_TRANSPARENT) || r > 0){
			*prq = sint1::falseI();
            c.v[0] += sfloat1::Or(sfloat1::And(rm,ll.v[0]),sfloat1::AndNot(rm,zr));
            c.v[1] += sfloat1::Or(sfloat1::And(rm,ll.v[1]),sfloat1::AndNot(rm,zr));
            c.v[2] += sfloat1::Or(sfloat1::And(rm,ll.v[2]),sfloat1::AndNot(rm,zr));
            c.v[3] += sfloat1::one();
		}
	}

	for(uint i = 0; i < BLCLOUD_VSIZE; ++i)
		ls.ls[r][i].clear();

	return c;
}

static void K_Render(dmatrix44 *pviewi, dmatrix44 *pproji, RenderKernel *pkernel, uint x0, uint y0, uint rx, uint ry, uint w, uint h, uint samples, dfloat4 *pout){
    tbb::enumerable_thread_specific<ParallelLeafList> leafs; //list here to avoid memory allocations
#define BLCLOUD_MT
#ifdef BLCLOUD_MT
    tbb::parallel_for(tbb::blocked_range2d<size_t>(y0,y0+ry/BLCLOUD_VY,x0,x0+rx/BLCLOUD_VX),[&](const tbb::blocked_range2d<size_t> &nr){
		for(uint y = nr.rows().begin(); y < nr.rows().end(); ++y){
			for(uint x = nr.cols().begin(); x < nr.cols().end(); ++x){
#else
	{
		for(uint y = y0; y < y0+ry; ++y){
			for(uint x = x0; x < x0+rx/BLCLOUD_VSIZE; ++x){
#endif
				sfloat4 posh;
                posh.v[0] = -(2.0f*(sfloat1((float)(BLCLOUD_VX*(x-x0)+x0)+0.5f)+sfloat1(0,1,0,1))/sfloat1((float)w)-sfloat1::one());
                posh.v[1] = -(2.0f*(sfloat1((float)(BLCLOUD_VY*(y-y0)+y0)+0.5f)+sfloat1(0,0,1,1))/sfloat1((float)h)-sfloat1::one());
                posh.v[2] = sfloat1::zero();
                posh.v[3] = sfloat1::one();

                matrix44 viewi = matrix44::load(pviewi);
                matrix44 proji = matrix44::load(pproji);

				sfloat4 ro = sfloat4(viewi.r[3]);
				sfloat4 rv = mul(posh,proji);
				sfloat4 rd = mul(rv.xyz0(),viewi).xyz0();
                rd = sfloat4::normalize3(rd);

                sint4 rngs;
                RNG_Init(&rngs);

				ParallelLeafList &ls = leafs.local();

                sfloat1 gm = sfloat1::And(
                    sfloat1::And(sfloat1::Greater(posh.v[0],-sfloat1::one()),sfloat1::Less(posh.v[0],sfloat1::one())),
                    sfloat1::And(sfloat1::Greater(posh.v[1],-sfloat1::one()),sfloat1::Less(posh.v[1],sfloat1::one())));

				sfloat1 rq;
                sfloat4 c = SampleVolume(ro,rd,gm,pkernel,&rngs,ls,0,samples,&rq);

				dintN wmask = dintN(gm);
                if(wmask.v[0] != 0)
                    float4::store(&pout[(BLCLOUD_VY*(y-y0)+0)*rx+BLCLOUD_VX*(x-x0)+0],c.get(0));
                if(wmask.v[1] != 0)
                    float4::store(&pout[(BLCLOUD_VY*(y-y0)+0)*rx+BLCLOUD_VX*(x-x0)+1],c.get(1));
                if(wmask.v[2] != 0)
                    float4::store(&pout[(BLCLOUD_VY*(y-y0)+1)*rx+BLCLOUD_VX*(x-x0)+0],c.get(2));
                if(wmask.v[3] != 0)
                    float4::store(&pout[(BLCLOUD_VY*(y-y0)+1)*rx+BLCLOUD_VX*(x-x0)+1],c.get(3));
			}
		}
#ifdef BLCLOUD_MT
	},tbb::auto_partitioner());
#else
	}
#endif
}

RenderKernel::RenderKernel(){
    //
}

RenderKernel::~RenderKernel(){
	//
}

bool RenderKernel::Initialize(const Scene *pscene, const dmatrix44 *pviewi, const dmatrix44 *pproji, const std::vector<Light> *plights, uint scattevs,
    float msigmas, float msigmaa, uint rx, uint ry, uint w, uint h, uint flags){
    if(!(phb = (dfloat4*)_mm_malloc(rx*ry*16,16))){
        //DebugPrintf("Framebuffer allocation failure.\n");
		return false;
    }

	this->pscene = pscene;
    this->scattevs = scattevs;
	this->msigmas = msigmas;
	this->msigmaa = msigmaa;
	this->rx = rx;
	this->ry = ry;
	this->w = w;
	this->h = h;
    this->flags = flags;
    viewi = *pviewi;
	proji = *pproji;

	this->lightc = plights->size();
	this->plights = new Light[this->lightc];
	memcpy(this->plights,plights->data(),this->lightc*sizeof(Light));

    //dfloat3 sdir = dfloat3(0.0f,0.0f,1.0f);
    float th = 0.0f;//acosf(sdir.z);
    //float ph = atan2f(sdir.x,sdir.y);
    float se = XM_PIDIV2-th; //solar elevation

	pskyms = arhosek_rgb_skymodelstate_alloc_init(2.2,0.6,se);

	DebugPrintf("Initialized render kernel.\n");

	return true;
}

void RenderKernel::Render(uint x0, uint y0, uint samples){
    //Wrapper to launch gpu-kernels (there have been a few tests).
    /*dim3 db = dim3(8,8,1); //[numthreads(...)]
    dim3 dg = dim3(rx/8,ry/8,1); //Dispatch
	K_Render<<<db,dg>>>(*(matrix*)&viewi,*(matrix*)&proji,ob,x0,y0,w,h,prt);

	if(cudaMemcpy(phb,prt,rx*ry*16,cudaMemcpyDeviceToHost) != cudaSuccess)
		printf("cudaMemcpy() failure\n");*/

    K_Render(&viewi,&proji,this,x0,y0,rx,ry,w,h,samples,phb);
}

void RenderKernel::Destroy(){
	arhosekskymodelstate_free(pskyms);
	/*cudaFree(ob);
	cudaFree(prt);*/
	delete []plights;
    _mm_free(phb);
}
