#if defined(_MSC_VER)
#pragma once
#endif

#ifndef VOLUMEPHOTON_H
#define VOLUMEPHOTON_H

#include "volume.h"
#include "integrator.h"
#include "kdtree.h"

class VolumePhotonIntegrator : public VolumeIntegrator
{
public:
	VolumePhotonIntegrator(int ncaus, int nindir, int nLookup, int mdepth,
			 float maxdist, bool finalGather, int gatherSamples,
			 float rrt, float ga);
	~VolumePhotonIntegrator();

	void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    Spectrum Li(const Scene *scene, const Renderer *renderer,
            const RayDifferential &ray, const Sample *sample, RNG &rng,
            Spectrum *transmittance, MemoryArena &arena) const;
    
	Spectrum Transmittance(const Scene *scene, const Renderer *,
            const RayDifferential &ray, const Sample *sample, RNG &rng,
            MemoryArena &arena) const;

	void Preprocess(const Scene *scene, const Camera *camera,
            const Renderer *renderer);

private:

};





#endif //VOLUMEPHOTON_H
