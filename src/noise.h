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

#endif // NOISE_H
