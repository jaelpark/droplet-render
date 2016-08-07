#ifndef NOISE_H
#define NOISE_H

namespace PerlinNoise{
sfloat1 noise(const sfloat4 &); //range [-1,1]
//float noise(const float4 &);
}

namespace fBm{
sfloat1 noise(const sfloat4 &, uint, float, float, float, float);
//float noise(const float4 &, ...); //vectorized noise octaves
float GetAmplitudeMax(uint, float, float);
}

namespace Node{

class ScalarFbmNoise : public IScalarFbmNoise{
public:
	ScalarFbmNoise(uint, NodeTree *);
	~ScalarFbmNoise();
	void Evaluate(const void *);
	enum OUTPUT{
		OUTPUT_NOISE,
		OUTPUT_MAXIMUM
	};
};

}

#endif // NOISE_H
