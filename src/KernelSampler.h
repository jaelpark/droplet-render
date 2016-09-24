#ifndef KERNEL_SAMPLER_H
#define KERNEL_SAMPLER_H

namespace KernelSampler{

class PhaseFunction{
public:
	PhaseFunction();
	~PhaseFunction();
	//TODO: color channel parameter for both below
	virtual sfloat1 Evaluate(const sfloat1 &) const = 0; //phase = pdf
	virtual sfloat4 Sample(const sfloat4 &, const sfloat1 &, const sfloat1 &) const = 0;
};

class HGPhase : public PhaseFunction{
public:
	HGPhase(float);
	~HGPhase();
	sfloat1 Evaluate(const sfloat1 &) const;
	sfloat4 Sample(const sfloat4 &, const sfloat1 &, const sfloat1 &) const;
	static HGPhase ghg;
private:
	float g1;
};

class MiePhase : public PhaseFunction{
public:
	MiePhase();
	~MiePhase();
	sfloat1 Evaluate(const sfloat1 &) const;
	sfloat4 Sample(const sfloat4 &, const sfloat1 &, const sfloat1 &) const; //sample from data relative to incident vector, then uniformly sample an azimuthal angle
	//Note: spectral rendering is not supported. To get MIE effects (for example), color channels need to be rendered separately with different
	//phase functions. Alternatively approximate the effect by assuming the different channels of the PDF to be close to each other, so that only
	//one can be used for sampling, while the full spectrum is evaluated for the color effects.
	static MiePhase gmie;
};

class BaseLight{
public:
	BaseLight();
	~BaseLight();
	virtual sfloat4 Evaluate(const sfloat4 &) const = 0; //evaluate radiance for some direction
	virtual sfloat1 Pdf(const sfloat4 &) const = 0;
	virtual sfloat4 Sample(const sfloat4 &, const sfloat1 &, const sfloat1 &) const = 0;
};

class SunLight : public BaseLight{
public:
	SunLight(const dfloat3 *, const dfloat3 *, float);
	~SunLight();
	sfloat4 Evaluate(const sfloat4 &) const;
	sfloat1 Pdf(const sfloat4 &) const;
	sfloat4 Sample(const sfloat4 &, const sfloat1 &, const sfloat1 &) const;
	dfloat3 direction;
	dfloat3 color;
	float angle; //cross-section angle
};

/*class EnvironmentMap : public BaseLight{
public:
	EnvironmentMap();
	~EnvironmentMap();
	sfloat4 Evaluate(const sfloat4 &) const;
	sfloat1 Pdf(const sfloat4 &) const;
	sfloat4 Sample(const sfloat4 &, const sfloat1 &, const sfloat1 &) const;
};*/

}

#endif
