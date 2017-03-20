#ifndef KERNEL_SAMPLER_H
#define KERNEL_SAMPLER_H

namespace KernelSampler{

class PhaseFunction{
public:
	PhaseFunction();
	~PhaseFunction();
	//TODO: color channel parameter for both below
	virtual sfloat1 Evaluate(const sfloat1 &) const = 0; //phase = pdf
	virtual sfloat4 EvaluateRGB(const sfloat1 &) const = 0;
	virtual sfloat4 Sample(const sfloat4 &, const sfloat1 &, const sfloat1 &) const = 0;
};

class HGPhase : public PhaseFunction{
public:
	HGPhase(float);
	~HGPhase();
	sfloat1 Evaluate(const sfloat1 &) const;
	sfloat4 EvaluateRGB(const sfloat1 &) const;
	sfloat4 Sample(const sfloat4 &, const sfloat1 &, const sfloat1 &) const;
	static HGPhase ghg;
private:
	float g1; //anisotropy parameter
};

class MiePhase : public PhaseFunction{
public:
	MiePhase();
	~MiePhase();
	sfloat1 Evaluate(const sfloat1 &) const; //In future if spectral rendering is supported, an additional parameter indicating the color channel is needed.
	sfloat4 EvaluateRGB(const sfloat1 &) const;
	sfloat4 Sample(const sfloat4 &, const sfloat1 &, const sfloat1 &) const;
	static MiePhase gmie;
};

class BaseLight{
public:
	BaseLight();
	~BaseLight();
	virtual sfloat4 Evaluate(const sfloat4 &) const = 0; //evaluate radiance for some direction
	virtual sfloat1 Pdf(const sfloat4 &) const = 0;
	virtual sfloat4 Sample(const sfloat4 &, const sfloat1 &, const sfloat1 &) const = 0;
	static void DeleteAll();
	static std::vector<BaseLight *> lights;
};

class SunLight : public BaseLight{
public:
	SunLight(const dfloat3 *, const dfloat3 *, float);
	~SunLight();
	sfloat4 Evaluate(const sfloat4 &) const;
	sfloat1 Pdf(const sfloat4 &) const;
	sfloat4 Sample(const sfloat4 &, const sfloat1 &, const sfloat1 &) const;
	dfloat3 direction;
	dfloat3 color; //color*intensity
	float angle; //cross-section angle
	float cosAngle;
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
