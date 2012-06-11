#if defined(_MSC_VER)
#pragma once
#endif

#ifndef VOLUMEPHOTON_H
#define VOLUMEPHOTON_H

#include "volume.h"
#include "integrator.h"
#include "kdtree.h"


#ifndef PHOTON_DEFINED
struct Photon;
struct ClosePhoton;
struct PhotonProcess;
#endif //PHOTON_DEFINED


class VolumePhotonIntegrator : public VolumeIntegrator
{
public:
	VolumePhotonIntegrator(int nvol, int nL, int mspecdepth, int mphotondepth,
		float mdist, bool fG, int gS, float ga, float ss, float cs)

	~VolumePhotonIntegrator();

	void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);

    Spectrum VolumePhotonIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Sample *sample, RNG &rng,
        Spectrum *T, MemoryArena &arena) const;

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

	float causticScale;

	int nVolumePaths;
	KdTree<Photon> *volumeMap;
};

//VolumePhotonIntegrator *CreateVolumePhotonMapSurfaceIntegrator(const ParamSet &params);


#endif //VOLUMEPHOTON_H
