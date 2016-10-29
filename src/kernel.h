#ifndef KERNEL_H
#define KERNEL_H

//#define RENDER_TRANSPARENT 0x1

namespace KernelSampler{
class PhaseFunction;
}

#define BLCLOUD_MAX_RECURSION 32 //BLCLOUD_MAX_RECURSION
class ParallelLeafList{
public:
	ParallelLeafList(){}
	~ParallelLeafList(){}

	inline uint GetLeafCount(uint r, uint v) const{
		return ls[r][v].size();
	}

	inline uint GetLeaf(uint r, uint v, uint x) const{
		return std::get<0>(ls[r][v][x]);
	}

	inline void GetHit(uint r, uint v, uint x, float *ptr0, float *ptr1) const{
		*ptr0 = std::get<1>(ls[r][v][x]);
		*ptr1 = std::get<2>(ls[r][v][x]);
	};

	typedef std::tuple<uint, float, float> Node;
	std::vector<Node> ls[BLCLOUD_MAX_RECURSION][BLCLOUD_VSIZE];
};

class RenderKernel{
public:
	RenderKernel();
	~RenderKernel();
	bool Initialize(const Scene *, const SceneOcclusion *, const dmatrix44 *, const dmatrix44 *, KernelSampler::PhaseFunction *, uint, float, float, uint, uint, uint, uint, uint);
	void Render(uint, uint, uint, uint, uint);
	void Destroy();
	//
	dfloat4 *phb; //host buffer
	const Scene *pscene;
	const SceneOcclusion *psceneocc;
	struct ArHosekSkyModelState *pskyms;
	KernelSampler::PhaseFunction *ppf;
	dmatrix44 viewi;
	dmatrix44 proji;
	//uint samples;
	uint scattevs; //max number of scattering events
	float msigmas; //macroscopic scattering cross section
	float msigmaa; //-- absorption
	//uint tilex; //max tile size
	//uint tiley;
	uint w;
	uint h;
	uint flags;
};

#endif
