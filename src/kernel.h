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
	bool Initialize(const class Scene *, const class SceneOcclusion *, const dmatrix44 *, const dmatrix44 *, KernelSampler::PhaseFunction *, KernelSampler::BaseEnv *, uint, float, float, uint, uint, uint, uint, uint);
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
	const class Scene *pscene;
	const class SceneOcclusion *psceneocc;
	struct ArHosekSkyModelState *pskyms;
	KernelSampler::PhaseFunction *ppf;
	KernelSampler::BaseEnv *penv;
	dmatrix44 viewi;
	dmatrix44 proji;
	//
	dfloat3 skydir;
	//uint samples;
	uint scattevs; //max number of scattering events
	float msigmas; //macroscopic scattering cross section
	float msigmaa; //-- absorption
	//
	uint w;
	uint h;
	//tile size
	uint tilew;
	uint tileh;
	//
	uint flags;
};

#endif
