#include "stdafx.h"
#include "integrators/volumephoton.h"
#include "paramset.h"

VolumePhotonIntegrator::VolumePhotonIntegrator(int ncaus, int nindir, 
		int nLookup, int mdepth, float maxdist, bool finalGather,
		int gatherSamples, float rrt, float ga)
{
}

VolumePhotonIntegrator::~VolumePhotonIntegrator()
{
}

void VolumePhotonIntegrator::RequestSamples(Sampler *sampler,
		Sample *sample, const Scene *scene)
{
}

Spectrum VolumePhotonIntegrator::Li(const Scene *scene, 
		const Renderer *renderer, const RayDifferential &ray, 
		const Sample *sample, RNG &rng, Spectrum *transmittance, 
		MemoryArena &arena) const
{
}

Spectrum VolumePhotonIntegrator::Transmittance(const Scene *scene, 
		const Renderer *, const RayDifferential &ray, 
		const Sample *sample, RNG &rng, MemoryArena &arena) const
{
}

void VolumePhotonIntegrator::Preprocess(const Scene *scene, 
		const Camera *camera, const Renderer *renderer)
{
}
