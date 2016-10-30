#include "main.h"
#include "scene.h"
#include "SceneOcclusion.h"
#include "kernel.h"
#include "KernelOctree.h"
#include "KernelSampler.h"
#include "ArHosekSkyModel.h"

#include <random>
//#include <fenv.h>

#include <tbb/parallel_for.h>

#define USE_SSE2
#include "sse_mathfun.h"

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

//_mm_ctf_epu32 from the web somewhere
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
	r.v[0] = v.v[0]*m.r[0].splatN<0>()+v.v[1]*m.r[1].splatN<0>()+v.v[2]*m.r[2].splatN<0>()+v.v[3]*m.r[3].splatN<0>();
	r.v[1] = v.v[0]*m.r[0].splatN<1>()+v.v[1]*m.r[1].splatN<1>()+v.v[2]*m.r[2].splatN<1>()+v.v[3]*m.r[3].splatN<1>();
	r.v[2] = v.v[0]*m.r[0].splatN<2>()+v.v[1]*m.r[1].splatN<2>()+v.v[2]*m.r[2].splatN<2>()+v.v[3]*m.r[3].splatN<2>();
	r.v[3] = v.v[0]*m.r[0].splatN<3>()+v.v[1]*m.r[1].splatN<3>()+v.v[2]*m.r[2].splatN<3>()+v.v[3]*m.r[3].splatN<3>();
	return r;
}

inline float SampleVoxelSpace(const float4 &p, float *pvol, const float4 &ce, uint lvoxc){
	float4 nv = float4((float)lvoxc);
	float4 ni = -0.5f*(ce-ce.splat<3>()-p)*(nv-float4::one())/ce.splat<3>();
	float4 nf = float4::max(float4::floor(ni),float4::zero()); //max() shouldn't be necessary?
	float4 nc = float4::min(float4::ceil(ni),nv-float4::one());
	float4 nl = ni-nf;

	//Safe-guard to clamp indices. There's a rare case where indices go out bounds.
	nc = float4::max(float4::min(nc,float4(lvoxc-1)),float4(0));
	nf = float4::max(float4::min(nf,float4(lvoxc-1)),float4(0));

	duint3 a = duint3((uint)nf.get<0>(),(uint)nf.get<1>(),(uint)nf.get<2>());
	duint3 b = duint3((uint)nc.get<0>(),(uint)nc.get<1>(),(uint)nc.get<2>());

	uint lvoxc2 = lvoxc*lvoxc;
	float4 ua = float4(pvol[a.z*lvoxc2+a.y*lvoxc+a.x],pvol[a.z*lvoxc2+b.y*lvoxc+a.x],
		pvol[b.z*lvoxc2+a.y*lvoxc+a.x],pvol[b.z*lvoxc2+b.y*lvoxc+a.x]); //([x0,y0,z0],[x0,y1,z0],[x0,y0,z1],[x0,y1,z1])
	float4 ub = float4(pvol[a.z*lvoxc2+a.y*lvoxc+b.x],pvol[a.z*lvoxc2+b.y*lvoxc+b.x],
		pvol[b.z*lvoxc2+a.y*lvoxc+b.x],pvol[b.z*lvoxc2+b.y*lvoxc+b.x]); //([x1,y0,z0],[x1,y1,z0],[x1,y0,z1],[x1,y1,z1])
	float4 u = float4::lerp(ua,ub,nl.splat<0>()); //([x,y0,z0],[x,y1,z0],[x,y0,z1],[x,y1,z1])

	float4 va = u.swizzle<0,2,0,2>();
	float4 vb = u.swizzle<1,3,1,3>();
	float4 v = float4::lerp(va,vb,nl.splat<1>()); //([x,y,z0],[x,y,z1],-,-)

	float4 wa = v.splat<0>();
	float4 wb = v.splat<1>();
	float4 w = float4::lerp(wa,wb,nl.splat<2>()); //([x,y,z],-,-,-)

	return w.get<0>();
}

static sfloat4 SampleVolume(sfloat4 ro, sfloat4 rd, sfloat1 gm, RenderKernel *pkernel, sint4 *prs, ParallelLeafList &ls1, uint r, uint samples, sfloat1 *prq){
	dfloatN maxd = dfloatN(MAX_OCCLUSION_DIST); //max path
	dintN ogm = dintN(0);
	if(pkernel->psceneocc)
		pkernel->psceneocc->Intersect(ro,rd,gm,&ogm,&maxd);

	dintN sgm = dintN(gm);
	KernelOctree::OctreeFullTraverser traverser;
	traverser.Initialize(ro,rd,sgm,maxd,&pkernel->pscene->ob);

	sfloat4 c = sfloat4::zero();

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

	sfloat1 msigmaa = sfloat1(pkernel->msigmaa);
	sfloat1 msigmas = sfloat1(pkernel->msigmas);
	sfloat1 msigmae = msigmaa+msigmas;

	for(uint s = 0; s < samples; ++s){
		sfloat1 qm = gm;
		sfloat1 rm = sint1::trueI();
		sfloat1 zr = sint1::falseI();

		sfloat1 td = sfloat1::zero(); //distance travelled
		for(uint i = 0;; ++i){
			/*dintN leafcount;
			for(uint j = 0; j < BLCLOUD_VSIZE; ++j)
				leafcount.v[j] = ls.GetLeafCount(r,j);

			qm = sfloat1::And(qm,rm);
			qm = sfloat1::And(qm,sint1::Less(sint1(i),sint1::load(&leafcount)));
			if(qm.AllFalse())
				break;*/

			duintN nodes;
			sfloat1 tra, trb;

			dintN mask = traverser.GetLeaf(i,&nodes,tra,trb);
			sint1 lm = sint1::load(&mask);
			qm = sfloat1::And(qm,rm);
			qm = sfloat1::And(qm,lm);
			if(qm.AllFalse())
				break;

			sfloat4 ce = sfloat4::zero();
			sfloat1 lo = td; //local origin

			for(uint j = 0; j < BLCLOUD_VSIZE; ++j)
				if(mask.v[j] != 0)
					ce.set(j,float4::load(&pkernel->pscene->ob[nodes.v[j]].ce));

			sfloat1 tr0 = tra-td;
			sfloat1 tr1 = trb-td;

			/*dfloatN TRA, TRB;
			for(uint j = 0; j < BLCLOUD_VSIZE; ++j)
				if(i < ls.GetLeafCount(r,j)){
					ls.GetHit(r,j,i,&TRA.v[j],&TRB.v[j]);
					ce.set(j,float4::load(&pkernel->pscene->ob[ls.GetLeaf(r,j,i)].ce));
				}
			sfloat1 tra = sfloat1::load(&TRA); //leaf distance from the ray origin
			sfloat1 trb = sfloat1::load(&TRB);
			sfloat1 tr0 = tra-td; //leaf distance from the current pointer
			sfloat1 tr1 = trb-td;*/

			sfloat4 r0 = ro+rd*tra;

			sint1 sm = sfloat1::Greater(tra,td);
			dintN SM = dintN(sm);

			dfloatN smax1, dist1, rho1;

			dintN QM = dintN(qm);
			dintN VM; //true: next leaf is a sole fog; false: sdf exists, but fog may not
			for(uint j = 0; j < BLCLOUD_VSIZE; ++j){
				if(QM.v[j] != 0){
					if(pkernel->pscene->ob[nodes.v[j]].volx[VOLUME_BUFFER_SDF] != ~0u){
						smax1.v[j] = 1.0f;
						VM.v[j] = 0;
					}else{
						smax1.v[j] = pkernel->pscene->ob[nodes.v[j]].qval[VOLUME_BUFFER_FOG];
						VM.v[j] = 1;
					}

					if(SM.v[j] != 0 && VM.v[j] == 0)
						dist1.v[j] = SampleVoxelSpace(r0.get(j),pkernel->pscene->pvol[VOLUME_BUFFER_SDF]+pkernel->pscene->lvoxc3*pkernel->pscene->ob[nodes.v[j]].volx[VOLUME_BUFFER_SDF],ce.get(j),pkernel->pscene->lvoxc);
					else dist1.v[j] = 1.0f;
				}else VM.v[j] = 0;
			}
			sfloat1 d0 = sfloat1::load(&dist1);
			sint1 vm = sint1::load(&VM);

			sm = sfloat1::And(sm,sfloat1::Greater(d0,zr));
			sm = sfloat1::Or(sm,vm); //skip if the next leaf is a sole fog

			sm = sfloat1::And(sm,qm); //can't allow any further changes in 'td' if qm == 0
			td = sfloat1::Or(sfloat1::And(sm,tra),sfloat1::AndNot(sm,td));

			sfloat1 s0 = sfloat1::Or(sfloat1::And(sm,tr0),sfloat1::AndNot(sm,zr)); //positive distance skipped (moving to next leaf)

			sm = qm;

			sfloat1 smax = sfloat1::load(&smax1);//sfloat1(1.0f); //local max in this leaf
			for(sfloat1 sr = -log_ps(RNG_Sample(prs))/(msigmae*smax), sc, sh;; sr -= log_ps(RNG_Sample(prs))/(msigmae*smax)){
				sm = sfloat1::And(sm,rm);
				if(sfloat1(sm).AllFalse())
					break;
				sc = s0+sr;
				td = sfloat1::Or(sfloat1::And(sm,lo+sc),sfloat1::AndNot(sm,td));

				sm = sfloat1::And(sm,sfloat1::Less(sc,tr1)); //check if out of extents
				sh = sfloat1::Or(sm,sfloat1::AndNot(rm,sint1::trueI())); //prevent modifications if rm == false
				td = sfloat1::Or(sfloat1::AndNot(sh,trb),sfloat1::And(sh,td));

				//if sm == false, this is always true (so rm won't be changed)
				//-> unless the ray has run out of leafs
				rm = sfloat1::Or(sfloat1::And(sm,sfloat1::Greater(sc,tr0)),sfloat1::AndNot(sm,rm));

				sfloat4 rc = ro+rd*td;

				sh = sfloat1::And(sm,rm);
				dintN SH = dintN(sh);
				for(uint j = 0; j < BLCLOUD_VSIZE; ++j){
					if(SH.v[j] != 0){
						if(VM.v[j] == 0){
							dist1.v[j] = SampleVoxelSpace(rc.get(j),pkernel->pscene->pvol[VOLUME_BUFFER_SDF]+pkernel->pscene->lvoxc3*pkernel->pscene->ob[nodes.v[j]].volx[VOLUME_BUFFER_SDF],ce.get(j),pkernel->pscene->lvoxc);
							//have to check if fog exists
							if(dist1.v[j] > 0.0f && pkernel->pscene->ob[nodes.v[j]].volx[VOLUME_BUFFER_FOG] != ~0u)
								rho1.v[j] = SampleVoxelSpace(rc.get(j),pkernel->pscene->pvol[VOLUME_BUFFER_FOG]+pkernel->pscene->lvoxc3*pkernel->pscene->ob[nodes.v[j]].volx[VOLUME_BUFFER_FOG],ce.get(j),pkernel->pscene->lvoxc);
							else rho1.v[j] = -1.0f;
						}else{
							dist1.v[j] = 1.0f;
							rho1.v[j] = SampleVoxelSpace(rc.get(j),pkernel->pscene->pvol[VOLUME_BUFFER_FOG]+pkernel->pscene->lvoxc3*pkernel->pscene->ob[nodes.v[j]].volx[VOLUME_BUFFER_FOG],ce.get(j),pkernel->pscene->lvoxc);
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
			}
		}

		//sample E(rc)/T here
		//T should be the value at rc; T(rc)

		sfloat4 ll; //total lighting (directional+sky)
		//skip (sky)lighting calculations if all the incident rays scatter (don't reach sun or sky) (at least one of rm != 0)
		//if(rm.AnyTrue() && !(pkernel->flags & RENDER_TRANSPARENT && r == 0)){
		if(rm.AnyTrue() && r > 0){
			sfloat4 lc = sfloat4::zero(); //total directional lighting
			for(uint i = 0, n = KernelSampler::BaseLight::lights.size(); i < n; ++i)
				lc += KernelSampler::BaseLight::lights[i]->Evaluate(rd);

			//skylighting
			sfloat1 rdz = sfloat1::abs(rd.v[2]); //add/remove abs to remove/get ground
			//sfloat1 rdz = rd.v[2]; //add/remove abs to remove/get ground
			sfloat1 sth = sfloat1::sqrt(1.0f-rd.v[2]*rd.v[2]);
			sfloat1 cth = rdz;
			sfloat1 lth = sfloat1::zero(); //sun theta
			sfloat1 lph = sfloat1::zero();
			sfloat1 slth, clth, slph, clph;
			sincos_ps(lth,&slth.v,&clth.v);
			sincos_ps(lph,&slph.v,&clph.v);

			sfloat1 rdd = rd.v[0]/rd.v[1];
			sfloat1 acr = 1.0f/sfloat1::sqrt(rdd*rdd+1.0f); //cos(ph = atan(rdd))
			sfloat1 czp = acr*(rdd*slph+clph); //cos(zr-ph)
			sfloat1 cph = sth*slth*czp+cth*clth;

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
				sfloat1 miem = (sfloat1::one()+raym)/sfloat1::pow(sfloat1::one()+caf[8]*caf[8]-2.0f*caf[8]*gacs,sfloat1(1.5f));
				sfloat1 zenh = sfloat1::sqrt(thcs);
				ca.v[i] = (sfloat1::one()+caf[0]*exp_ps(caf[1]/(thcs+0.01f)))*(caf[2]+caf[3]*expm+caf[5]*raym+caf[6]*miem+caf[7]*zenh);
				//ca.v[i] *= 0.028f*pkernel->pskyms->radiances[i];
				ca.v[i] *= pkernel->pskyms->radiances[i];
				ca.v[i] = 0.00035f*sfloat1::pow(ca.v[i],2.2f); //convert to linear and adjust exposure
			}

			ll = lc+ca;
			ll.v[3] = sfloat1::one(); //alpha doesn't matter when r > 0
		}else{
			ll = sfloat4::zero();
			ll.v[3] = sfloat1::AndNot(rm,sfloat1::one()); //alpha = 1 when scattering
		}

		if(r < pkernel->scattevs && rm.AnyFalse()){
#define BLCLOUD_MULTIPLE_IMPORTANCE
#ifdef BLCLOUD_MULTIPLE_IMPORTANCE
			//Sample the phase function and lights
			//TODO: choose randomly one the lights. Multiply the final estimate (s2) with the total number of lights (ref776).
			sfloat1 u1 = RNG_Sample(prs), u2 = RNG_Sample(prs);
			sfloat4 srd = pkernel->ppf->Sample(rd,u1,u2);//HG_Sample(rd,prs);

			sfloat1 u3 = RNG_Sample(prs), u4 = RNG_Sample(prs);
			sfloat4 lrd = KernelSampler::BaseLight::lights[0]->Sample(rd,u3,u4);

			//pdfs for the balance heuristic w_x = p_x/sum(p_i,i=0..N)
			sfloat4 p1 = pkernel->ppf->EvaluateRGB(sfloat4::dot3(srd,rd));//HG_Phase(sfloat4::dot3(srd,rd));
			sfloat1 p2 = KernelSampler::BaseLight::lights[0]->Pdf(lrd);//L_Pdf(lrd,la);

			//need two samples - with the phase sampling keep on the recursion while for the light do only single scattering
			sfloat1 gm1 = sfloat1::AndNot(rm,sint1::trueI());
			sfloat1 rq;

			sfloat4 rc = ro+rd*td;

			//estimator S(1)*f1*w1/p1+S(2)*f2*w2/p2 /= woodcock pdf
			sfloat4 s1 = SampleVolume(rc,srd,gm1,pkernel,prs,ls1,r+1,1,&rq);
			sfloat4 s2 = SampleVolume(rc,lrd,gm1,pkernel,prs,ls1,BLCLOUD_MAX_RECURSION-1,1,&rq);

			//(HG_Phase(X)*SampleVolume(X)*p1/(p1+L_Pdf(X)))/p1 => (HG_Phase(X)=p1)*SampleVolume(X)/(p1+L_Pdf(X)) = p1*SampleVolume(X)/(p1+L_Pdf(X))
			//(HG_Phase(Y)*SampleVolume(Y)*p2/(HG_Phase(Y)+p2))/p2 => HG_phase(Y)*SampleVolume(Y)/(HG_Phase(Y)+p2)
			sfloat4 p3 = pkernel->ppf->EvaluateRGB(sfloat4::dot3(lrd,rd));
			sfloat4 cm = s1*p1/(p1+KernelSampler::BaseLight::lights[0]->Pdf(srd))+s2*p3/(p3+p2);//s1*p1/(p1+L_Pdf(srd,la))+s2*p3/(p3+p2);
			cm *= msigmas/msigmae;

			//s1=HG_Phase(X)*SampleVolume(X)
			//s1*p1/(p1+L_Pdf(X))
#else
			sfloat1 rq;
			//phase function sampling
			sfloat1 u1 = RNG_Sample(prs), u2 = RNG_Sample(prs);
			sfloat4 srd = pkernel->ppf->Sample(rd,u1,u2);
			sfloat4 cm = SampleVolume(rc,srd,sfloat1::AndNot(rm,sint1::trueI()),pkernel,prs,ls,r+1,1,&rq)*msigmas/msigmae;
			//phase/pdf(=phase)=1

			//light sampling, obviously won't alone result in any sky lighting or proper multiple scattering
			/*sfloat1 la = sfloat1(pkernel->plights[0].angle);
			sfloat4 lrd = L_Sample(sfloat1(float4::load(&pkernel->plights[0].direction)),la,prs);
			sfloat4 cm = SampleVolume(rc,lrd,sfloat1::AndNot(rm,sint1::trueI()),pkernel,prs,ls,r+1,1)*HG_Phase(sfloat4::dot3(lrd,rd))*msigmas
				/(msigmae*L_Pdf(lrd,la));*/
#endif
			sfloat1 hm = sfloat1::Greater(cm.v[0],sfloat1(1e5f)); //temp hack to deal with some bug, will be properly fixed later
			cm.v[0] = sfloat1::Or(sfloat1::And(hm,c.v[0]/(float)(s+1)),sfloat1::AndNot(hm,cm.v[0]));
			cm.v[1] = sfloat1::Or(sfloat1::And(hm,c.v[1]/(float)(s+1)),sfloat1::AndNot(hm,cm.v[1]));
			cm.v[2] = sfloat1::Or(sfloat1::And(hm,c.v[2]/(float)(s+1)),sfloat1::AndNot(hm,cm.v[2]));

			*prq = sint1::trueI();
			c.v[0] += sfloat1::Or(sfloat1::And(rm,ll.v[0]),sfloat1::AndNot(rm,cm.v[0]));
			c.v[1] += sfloat1::Or(sfloat1::And(rm,ll.v[1]),sfloat1::AndNot(rm,cm.v[1]));
			c.v[2] += sfloat1::Or(sfloat1::And(rm,ll.v[2]),sfloat1::AndNot(rm,cm.v[2]));
			c.v[3] += ll.v[3];
		}else{
			*prq = sint1::falseI();
			c.v[0] += sfloat1::Or(sfloat1::And(rm,ll.v[0]),sfloat1::AndNot(rm,zr));
			c.v[1] += sfloat1::Or(sfloat1::And(rm,ll.v[1]),sfloat1::AndNot(rm,zr));
			c.v[2] += sfloat1::Or(sfloat1::And(rm,ll.v[2]),sfloat1::AndNot(rm,zr));
			c.v[3] += ll.v[3];
		}
	}

	for(uint i = 0; i < BLCLOUD_VSIZE; ++i)
		ls1.ls[r][i].clear();

	return c;
}

static void K_Render(dmatrix44 *pviewi, dmatrix44 *pproji, RenderKernel *pkernel, uint x0, uint y0, uint rx, uint ry, uint w, uint h, uint samples, dfloat4 *pout){
	//feenableexcept(FE_ALL_EXCEPT&~FE_INEXACT);
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
	//fedisableexcept(FE_ALL_EXCEPT&~FE_INEXACT);
}

RenderKernel::RenderKernel(){
	//
}

RenderKernel::~RenderKernel(){
	//
}

bool RenderKernel::Initialize(const Scene *pscene, const SceneOcclusion *psceneocc, const dmatrix44 *pviewi, const dmatrix44 *pproji, KernelSampler::PhaseFunction *ppf, uint scattevs,
	float msigmas, float msigmaa, uint tilex, uint tiley, uint w, uint h, uint flags){
	if(!(phb = (dfloat4*)_mm_malloc(tilex*tiley*16,16)))
		return false;

	this->pscene = pscene;
	this->psceneocc = psceneocc;
	this->scattevs = scattevs;
	this->msigmas = msigmas;
	this->msigmaa = msigmaa;
	this->w = w;
	this->h = h;
	this->flags = flags;
	viewi = *pviewi;
	proji = *pproji;

	this->ppf = ppf;

	//TODO: set direction according to primary light
	float th = 0.0f;//acosf(sdir.z);
	float se = XM_PIDIV2-th; //solar elevation

	pskyms = arhosek_rgb_skymodelstate_alloc_init(2.2,0.6,se);
	//pskyms = arhosek_rgb_skymodelstate_alloc_init(7.4,0.2,se);

	if(KernelSampler::BaseLight::lights.size() != 1)
		DebugPrintf("Warning: only one directional light is currently properly supported.\n");

	DebugPrintf("Initialized render kernel.\n");

	return true;
}

void RenderKernel::Render(uint x0, uint y0, uint tilex, uint tiley, uint samples){
	//
	K_Render(&viewi,&proji,this,x0,y0,tilex,tiley,w,h,samples,phb);
}

void RenderKernel::Destroy(){
	arhosekskymodelstate_free(pskyms);
	_mm_free(phb);
}
