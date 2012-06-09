
/*
    pbrt source code Copyright(c) 1998-2010 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    pbrt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.  Note that the text contents of
    the book "Physically Based Rendering" are *not* licensed under the
    GNU GPL.

    pbrt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

 */

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_VOLUMEPHOTONMAP_H
#define PBRT_INTEGRATORS_VOLUMEPHOTONMAP_H

#include "pbrt.h"
#include "integrator.h"
#include "volume.h"
#include "kdtree.h"

#ifndef USE_SURFACE_PHOTONS
struct Photon;
struct ClosePhoton;
struct PhotonProcess;
#endif


// VolumePhotonIntegrator Declarations
class VolumePhotonIntegrator : public VolumeIntegrator {
public:
    // VolumePhotonIntegrator Public Methods
    VolumePhotonIntegrator(int nvol, int nLookup, int maxspecdepth,
        int maxphotondepth, float maxdist, bool finalGather, int gatherSamples,
        float ga, float stepSize, float causticScale);
    ~VolumePhotonIntegrator();
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Sample *sample,
        RNG &rng, Spectrum *transmittance, MemoryArena &arena) const;
	Spectrum Transmittance(const Scene * scene, const Renderer * renderer, 
		const RayDifferential & ray, const Sample * sample, RNG & rng, 
		MemoryArena & arena) const;
	void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    void Preprocess(const Scene *scene, const Camera *camera, const Renderer *renderer);

	void WritePhotonsOutToMap(vector<Photon> * map);
private:
    // VolumePhotonIntegrator Private Methods
    friend class VolumePhotonShootingTask;

    // VolumePhotonIntegrator Private Data
	uint32_t nLookup;
    float maxDistSquared;
    int maxSpecularDepth, maxPhotonDepth;
    bool finalGather;
    int gatherSamples;
    float cosGatherAngle;

	float stepSize;

    int tauSampleOffset, scatterSampleOffset;

	//storing the volumetric photons
	KdTree<Photon> *volumeMap;
    uint32_t nVolumePhotonsWanted;
	int nVolumePaths;
	// region of water, fog etc.
	VolumeRegion * volumeRegion;

	float causticScale;
};


VolumePhotonIntegrator *CreateVolumePhotonMapSurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_VOLUMEPHOTONMAP_H
