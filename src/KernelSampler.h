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
	float g1; //anisotropy parameter
	static HGPhase ghg;
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

class BaseEnv{
public:
	BaseEnv();
	~BaseEnv();
	virtual sfloat4 Evaluate(const sfloat4 &) const = 0;
};

class NullEnv : public BaseEnv{
public:
	NullEnv();
	~NullEnv();
	virtual sfloat4 Evaluate(const sfloat4 &) const;
	static NullEnv nenv;
};

class MapEnv : public BaseEnv{
public:
	MapEnv();
	~MapEnv();
	dfloat4 * Initialize(uint, uint);
	sfloat4 Evaluate(const sfloat4 &) const;
	void Destroy();
/*private:
	float P(int, int, float) const;
	int fact(int);
	float K(int, int);
	float SH(int, int, float, float);
	sfloat4 Fetch(const sfloat4 &, uint, uint, const dfloat4 *) const;
public:
#define ENV_SHN 12
	dfloat3 c[ENV_SHN*ENV_SHN];*/
private:
	dfloat4 *ptex;
	uint w, h;
public:
	static MapEnv genv;
};

}

#endif
