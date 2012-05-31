#include "stdafx.h"
#include "integrators/volumephoton.h"
#include "paramset.h"

VolumePhotonIntegrator::VolumePhotonIntegrator(int nvol, int nLookup, 
		int maxspecdepth, int maxphotondepth, float maxdist, 
		bool finalGather, int gatherSamples, float ga, float ss)
		: nVolumePhotonsWanted(nvol)
		, nLookup(nLookup)
		, maxSpecularDepth(maxspecdepth)
		, maxPhotonDepth(maxphotondepth)
		, maxDistSquared(maxdist * maxdist)
		, finalGather(finalGather)
		, gatherSamples(gatherSamples)
		, cosGatherAngle(cos(Radians(ga)))
		, stepSize(ss)
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
	return 0.f;
}

Spectrum VolumePhotonIntegrator::Transmittance(const Scene *scene, 
		const Renderer *, const RayDifferential &ray, 
		const Sample *sample, RNG &rng, MemoryArena &arena) const
{
	return 0.f;
}

void VolumePhotonIntegrator::Preprocess(const Scene *scene, 
		const Camera *camera, const Renderer *renderer)
{
}
