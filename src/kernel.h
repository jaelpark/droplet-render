#ifndef KERNEL_H
#define KERNEL_H

#define RENDER_TRANSPARENT 0x1

class PhaseFunction{
public:
	PhaseFunction();
	~PhaseFunction();
	//TODO: color channel parameter for both below
	virtual sfloat1 Evaluate(const sfloat1 &) const = 0; //phase = pdf
	virtual sfloat4 Sample(const sfloat4 &, sint4 *) const = 0;
};

class HGPhase : public PhaseFunction{
public:
	HGPhase(float);
	~HGPhase();
	sfloat1 Evaluate(const sfloat1 &) const;
	sfloat4 Sample(const sfloat4 &, sint4 *) const;
	static HGPhase ghg;
private:
	float g1;
};

class MiePhase : public PhaseFunction{
public:
	MiePhase();
	~MiePhase();
	sfloat1 Evaluate(const sfloat1 &) const;
	sfloat4 Sample(const sfloat4 &, sint4 *) const; //sample from data relative to incident vector, then uniformly sample an azimuthal angle
	//Note: spectral rendering is not supported. To get MIE effects (for example), color channels need to be rendered separately with different
	//phase functions. Alternatively approximate the effect by assuming the different channels of the PDF to be close to each other, so that only
	//one can be used for sampling, while the full spectrum is evaluated for the color effects.
	//bool: spectral approximation
	static MiePhase gmie; //I suppose we have a singleton
};

struct Light{
	PhaseFunction *ppf;
    dfloat3 direction;
    dfloat3 color; //color*intensity
    float angle;
};

#define BLCLOUD_MAX_RECURSION 32 //BLCLOUD_MAX_RECURSION
class ParallelLeafList{
public:
	ParallelLeafList(){}
	~ParallelLeafList(){}

	inline uint GetLeafCount(uint r, uint v) const{
		return ls[r][v].size();
	}

	inline uint GetLeaf(uint r, uint v, uint x) const{
		return ls[r][v][x];
	}

    std::vector<uint> ls[BLCLOUD_MAX_RECURSION][BLCLOUD_VSIZE];
};

class RenderKernel{
public:
	RenderKernel();
	~RenderKernel();
    bool Initialize(const Scene *, const dmatrix44 *, const dmatrix44 *, const std::vector<Light> *, uint, float, float, uint, uint, uint, uint, uint);
    void Render(uint, uint, uint);
	void Destroy();
	//
    dfloat4 *phb; //host buffer
    const Scene *pscene;
	struct ArHosekSkyModelState *pskyms;
	Light *plights;
	uint lightc;
    dmatrix44 viewi;
    dmatrix44 proji;
    //uint samples;
    uint scattevs; //max number of scattering events
	float msigmas; //macroscopic scattering cross section
	float msigmaa; //-- absorption
	uint rx;
	uint ry;
	uint w;
	uint h;
    uint flags;
};

#endif
