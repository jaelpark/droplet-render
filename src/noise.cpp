#include "main.h"
#include "node.h"
#include "noise.h"

#define USE_SSE2
#include "sse_mathfun.h"

//http://mrl.nyu.edu/~perlin/noise/
namespace PerlinNoise{

static const int p[] = {151,160,137,91,90,15,
    131,13,201,95,96,53,194,233,7,225,140,36,103,30,69,142,8,99,37,240,21,10,23,
    190, 6,148,247,120,234,75,0,26,197,62,94,252,219,203,117,35,11,32,57,177,33,
    88,237,149,56,87,174,20,125,136,171,168, 68,175,74,165,71,134,139,48,27,166,
    77,146,158,231,83,111,229,122,60,211,133,230,220,105,92,41,55,46,245,40,244,
    102,143,54, 65,25,63,161, 1,216,80,73,209,76,132,187,208, 89,18,169,200,196,
    135,130,116,188,159,86,164,100,109,198,173,186, 3,64,52,217,226,250,124,123,
    5,202,38,147,118,126,255,82,85,212,207,206,59,227,47,16,58,17,182,189,28,42,
    223,183,170,213,119,248,152, 2,44,154,163, 70,221,153,101,155,167, 43,172,9,
    129,22,39,253, 19,98,108,110,79,113,224,232,178,185, 112,104,218,246,97,228,
    251,34,242,193,238,210,144,12,191,179,162,241, 81,51,145,235,249,14,239,107,
    49,192,214, 31,181,199,106,157,184, 84,204,176,115,121,50,45,127, 4,150,254,
    138,236,205,93,222,114,67,29,24,72,243,141,128,195,78,66,215,61,156,180
};

inline sfloat4 fade(const sfloat4 &t){
    //return t * t * t * (t * (t * 6 - 15) + 10);
    return t*t*t*(t*(t*sfloat1(6.0f)-sfloat1(15.0f))+sfloat1(10.0f));
}

inline float4 fade(const float4 &t){
    return t*t*t*(t*(t*float4(6.0f)-float4(15.0f))+float4(10.0f));
}

inline sfloat1 grad(const sint1 &_h, const sfloat4 &pos){
    sint1 h = sint1::And(_h,sint1(15));
    sint1 mu = sint1::Less(h,sint1(8));
    sint1 mv = sint1::Less(h,sint1(4));
    sint1 mw = sint1::Or(sint1::Equal(h,sint1(12)),sint1::Equal(h,sint1(14)));
    sfloat1 u = sfloat1::Or(sfloat1::And(mu,pos.v[0]),sfloat1::AndNot(mu,pos.v[1]));
    sfloat1 v = sfloat1::Or(sfloat1::And(mv,pos.v[1]),sfloat1::AndNot(mv,sfloat1::Or(sfloat1::And(mw,pos.v[0]),sfloat1::AndNot(mw,pos.v[2]))));
    sint1 ms = sint1::Equal(sint1::And(h,sint1(1)),sint1(0));
    sint1 mt = sint1::Equal(sint1::And(h,sint1(2)),sint1(0));
    return sfloat1::Or(sint1::And(ms,u),sint1::AndNot(ms,-u))+sfloat1::Or(sfloat1::And(mt,v),sfloat1::AndNot(mt,-v));
}

sfloat1 noise(const sfloat4 &_pos){
    sfloat4 F = sfloat4::floor(_pos);
    sfloat4 pos = _pos-F; //frac
    sfloat4 fad = fade(pos);

    sint4 I = sint4::And(sint4::convert(F),sint4(sint1(255)));
	dintN AA, AB, BA, BB;
	dintN X, Y, Z;
    sint1::store(&X,I.v[0]);
    sint1::store(&Y,I.v[1]);
    sint1::store(&Z,I.v[2]);
    for(uint i = 0; i < BLCLOUD_VSIZE; ++i){
        int A = p[(X.v[i]+0)%256]+Y.v[i];
        AA.v[i] = p[(A+0)%256]+Z.v[i];
        AB.v[i] = p[(A+1)%256]+Z.v[i];
        int B = p[(X.v[i]+1)%256]+Y.v[i];
        BA.v[i] = p[(B+0)%256]+Z.v[i];
        BB.v[i] = p[(B+1)%256]+Z.v[i];
    }
    sint1 aa = sint1::load(&AA), ab = sint1::load(&AB), ba = sint1::load(&BA), bb = sint1::load(&BB);
    return sfloat1::lerp(
		sfloat1::lerp(sfloat1::lerp(grad(aa,pos),
	                                 grad(ba,pos+sfloat1(-1,0,0,0)),fad.v[0]),
	                   sfloat1::lerp(grad(ab,pos+sfloat1(0,-1,0,0)),
	                                 grad(bb,pos+sfloat1(-1,-1,0,0)),fad.v[0]),fad.v[1]),
	    sfloat1::lerp(sfloat1::lerp(grad(aa+sfloat1::one(),pos+sfloat1(0,0,-1,0)),
	                                 grad(ba+sfloat1::one(),pos+sfloat1(-1,0,-1,0)),fad.v[0]),
	                   sfloat1::lerp(grad(ab+sfloat1::one(),pos+sfloat1(0,-1,-1,0)),
	                                 grad(bb+sfloat1::one(),pos+sfloat1(-1,-1,-1,0)),fad.v[0]),fad.v[1]),fad.v[2]);
}

}

namespace HashNoise{

sfloat1 hash(sfloat1 n){
	sfloat1 t = sin_ps(n)*43758.5453f;
	return t-sfloat1::floor(t);
}

sfloat1 noise(const sfloat4 &_pos){
	sfloat4 F = sfloat4::floor(_pos);
    sfloat4 pos = _pos-F; //frac

	sfloat4 fad = pos*pos*(sfloat4(sfloat1(3.0f))-sfloat4(sfloat1(2.0f))*pos);
	sfloat1 n = F.v[0]+57.0f*F.v[1]+113.0f*F.v[2];
	return sfloat1::lerp(
		sfloat1::lerp(sfloat1::lerp(hash(n),
	                                 hash(n+sfloat1(1.0f)),fad.v[0]),
	                   sfloat1::lerp(hash(n+sfloat1(57.0f)),
	                                 hash(n+sfloat1(58.0f)),fad.v[0]),fad.v[1]),
	    sfloat1::lerp(sfloat1::lerp(hash(n+sfloat1(113.0f)),
	                                 hash(n+sfloat1(114.0f)),fad.v[0]),
	                   sfloat1::lerp(hash(n+sfloat1(170.0f)),
	                                 hash(n+sfloat1(171.0f)),fad.v[0]),fad.v[1]),fad.v[2]);
}

}

namespace fBm{

sfloat1 noise(const sfloat4 &_pos, uint octaves, float freq, float amp, float fjump, float gain){
    sfloat1 s = sfloat1::zero();
    for(uint i = 0; i < octaves; ++i){
		//warning: PerlinNoise range -1,1 while hash gives 0,1
		s += PerlinNoise::noise(_pos*freq)*amp;
        //s += HashNoise::noise(_pos*freq)*amp;
        freq *= fjump;
        amp *= gain;
    }
    return s;
}

float GetAmplitudeMax(uint octaves, float amp, float gain){
    float s = 0.0f;
    for(uint i = 0; i < octaves; ++i){
        s += amp;
        amp *= gain;
    }
    return s;
}

}

namespace Node{

FbmNoise::FbmNoise(uint _level, NodeTree *pnt) : BaseValueNode<float>(_level,pnt), BaseValueNode<dfloat3>(_level,pnt), BaseNode(_level,pnt), IFbmNoise(_level,pnt){
	//
}

FbmNoise::~FbmNoise(){
	//
}

void FbmNoise::Evaluate(const void *pp){
	BaseValueNode<int> *poctn = dynamic_cast<BaseValueNode<int>*>(pnodes[INPUT_OCTAVES]);
    BaseValueNode<float> *pfreqn = dynamic_cast<BaseValueNode<float>*>(pnodes[INPUT_FREQ]);
    BaseValueNode<float> *pampn = dynamic_cast<BaseValueNode<float>*>(pnodes[INPUT_AMP]);
    BaseValueNode<float> *pfjumpn = dynamic_cast<BaseValueNode<float>*>(pnodes[INPUT_FJUMP]);
    BaseValueNode<float> *pgainn = dynamic_cast<BaseValueNode<float>*>(pnodes[INPUT_GAIN]);
    BaseValueNode<dfloat3> *pnode = dynamic_cast<BaseValueNode<dfloat3>*>(pnodes[INPUT_POSITION]);

	dfloat3 dposw = pnode->locr(indices[INPUT_POSITION]);
	sfloat4 sposw = sfloat4(sfloat1(dposw.x,dposw.y,dposw.z,0.0f))+sfloat4(
		sfloat1(0.0f),
		sfloat1(1.0f,0.0f,0.0f,0.0f), //TODO: module the offset somehow
		sfloat1(0.0f,1.0f,0.0f,0.0f),
		sfloat1(0.0f,0.0f,1.0f,0.0f))*sfloat1(154.7f);

	sfloat1 f = fBm::noise(sposw,poctn->locr(indices[INPUT_OCTAVES]),pfreqn->locr(indices[INPUT_FREQ]),
		pampn->locr(indices[INPUT_AMP]),pfjumpn->locr(indices[INPUT_FJUMP]),pgainn->locr(indices[INPUT_GAIN]));

	BaseValueResult<float> &rs = this->BaseValueNode<float>::result.local();
	rs.value[OUTPUT_FLOAT_NOISE] = f.get<0>();
	rs.value[OUTPUT_FLOAT_MAXIMUM] = fBm::GetAmplitudeMax(poctn->locr(indices[INPUT_OCTAVES]),pampn->locr(indices[INPUT_AMP]),
		pgainn->locr(indices[INPUT_GAIN]));
	BaseValueResult<dfloat3> &rv = this->BaseValueNode<dfloat3>::result.local();
	rv.value[OUTPUT_VECTOR_NOISE] = dfloat3(0.57735f*f); //normalize by 1/sqrt(3) to have max length equal to amplitude
}

IFbmNoise * IFbmNoise::Create(uint level, NodeTree *pnt){
	return new FbmNoise(level,pnt);
}

}
