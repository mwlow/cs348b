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
	VolumePhotonIntegrator(int nvol, int nLookup, int maxspecdepth,
		int maxphotondepth, float maxdist, bool finalGather, 
		int gatherSamples, float ga, float ss);
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
	friend class PhotonShootingTask;	

	uint32_t nVolumePhotonsWanted, nLookup;
	int maxSpecularDepth, maxPhotonDepth;
	float maxDistSquared;
	bool finalGather;
	int gatherSamples;
	float cosGatherAngle;

	float stepSize;
	int tauSampleOffset, scatterSampleOffset;

	//KdTree<Photon> *volumeMap;
};

//VolumePhotonIntegrator *CreateVolumePhotonMapSurfaceIntegrator(const ParamSet &params);


#endif //VOLUMEPHOTON_H
