#ifndef KERNEL_H
#define KERNEL_H

//#define RENDER_TRANSPARENT 0x1

namespace KernelSampler{
class PhaseFunction;
class BaseEnv;
}

class RenderKernel{
public:
	RenderKernel();
	~RenderKernel();
	bool Initialize(const class Scene *, const class SceneOcclusion *, const dmatrix44 *, const dmatrix44 *, KernelSampler::PhaseFunction *, KernelSampler::BaseEnv *, float *, uint, float, float, uint, uint, uint, uint, uint);
	void Render(uint, uint, uint, uint, uint);
	void Shadow(uint, uint, uint, uint, uint);
	void Destroy();
	//
	enum BUFFER{
		BUFFER_LIGHTS,
		BUFFER_ENVIRONMENT,
		BUFFER_COUNT
	};

	dfloat4 *phb[BUFFER_COUNT]; //host buffer
	float *pdepth; //source depth for compositing shadow calculations

	const class Scene *pscene;
	const class SceneOcclusion *psceneocc;

#ifdef USE_ARHOSEK_SKYMODEL
	struct ArHosekSkyModelState *pskyms;
#endif

	KernelSampler::PhaseFunction *ppf;
	KernelSampler::BaseEnv *penv;

	dmatrix44 viewi;
	dmatrix44 proji;
	//
#ifdef USE_ARHOSEK_SKYMODEL
	dfloat3 skydir;
#endif
	//uint samples;
	uint scattevs; //max number of scattering events
	float msigmas; //macroscopic scattering cross section
	float msigmaa; //-- absorption
	//
	uint w;
	uint h;
	//tile size (from last Render())
	uint tilew;
	uint tileh;
	//
	uint flags;
	//
	static dint3 vpattern[BLCLOUD_VSIZE]; //vectorized pixel pattern
};

#endif
