#ifndef KERNEL_H
#define KERNEL_H

//#define RENDER_TRANSPARENT 0x1

namespace KernelSampler{
class PhaseFunction;
}

class RenderKernel{
public:
	RenderKernel();
	~RenderKernel();
	bool Initialize(const class Scene *, const class SceneOcclusion *, const dmatrix44 *, const dmatrix44 *, KernelSampler::PhaseFunction *, uint, float, float, uint, uint, uint, uint, uint);
	void Render(uint, uint, uint, uint, uint);
	void Destroy();
	//
	dfloat4 *phb; //host buffer
	const class Scene *pscene;
	const class SceneOcclusion *psceneocc;
	struct ArHosekSkyModelState *pskyms;
	KernelSampler::PhaseFunction *ppf;
	dmatrix44 viewi;
	dmatrix44 proji;
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
