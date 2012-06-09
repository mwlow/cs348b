
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

#include "stdafx.h"
#include "integrators/volumephotonmap.h"
#include "scene.h"
#include "montecarlo.h"
#include "sampler.h"
#include "progressreporter.h"
#include "intersection.h"
#include "paramset.h"
#include "camera.h"

#ifdef USE_SURFACE_PHOTONS
#include "integrators/photonmap.h"
#endif

#define VOLUME_PHOTONS	50000
#define VOLUME_LOOKUP	100
#define RAY_EPSILON		0.00001
#define CAUSTIC_SCALE	0.5

#include <fstream>

// PhotonIntegrator Local Declarations FIRST DEFINED IN PHOTONMAP.CPP/.H

#ifndef USE_SURFACE_PHOTONS
struct Photon {
    Photon(const Point &pp, const Spectrum &wt, const Vector &w)
        : p(pp), alpha(wt), wi(w) { }
    Photon() { }
    Point p;
    Spectrum alpha;
    Vector wi;
};

struct PhotonProcess {
    // PhotonProcess Public Methods
    PhotonProcess(uint32_t mp, ClosePhoton *buf);
    void operator()(const Point &p, const Photon &photon, float dist2,
         float &maxDistSquared);
    ClosePhoton *photons;
    uint32_t nLookup, nFound;
};


struct ClosePhoton {
    // ClosePhoton Public Methods
    ClosePhoton(const Photon *p = NULL, float md2 = INFINITY)
        : photon(p), distanceSquared(md2) { }
    bool operator<(const ClosePhoton &p2) const {
        return distanceSquared == p2.distanceSquared ?
            (photon < p2.photon) : (distanceSquared < p2.distanceSquared);
    }
    const Photon *photon;
    float distanceSquared;
};

/*
PhotonProcess::PhotonProcess(uint32_t mp, ClosePhoton *buf) {
    photons = buf;
    nLookup = mp;
    nFound = 0;
}

inline void PhotonProcess::operator()(const Point &p,
        const Photon &photon, float distSquared, float &maxDistSquared) {
    if (nFound < nLookup) {
        // Add photon to unordered array of photons
        photons[nFound++] = ClosePhoton(&photon, distSquared);
        if (nFound == nLookup) {
            std::make_heap(&photons[0], &photons[nLookup]);
            maxDistSquared = photons[0].distanceSquared;
        }
    }
    else {
        // Remove most distant photon from heap and add new photon
        std::pop_heap(&photons[0], &photons[nLookup]);
        photons[nLookup-1] = ClosePhoton(&photon, distSquared);
        std::push_heap(&photons[0], &photons[nLookup]);
        maxDistSquared = photons[0].distanceSquared;
    }
}
*/
#endif


/****************************** CODE COPIED FROM PHOTONMAP.CPP ******************************/
inline float kernel(const Photon *photon, const Point &p, float maxDist2);

static Spectrum LVolumePhoton(KdTree<Photon> *map, int nPaths, int nLookup,
    ClosePhoton *lookupBuf, VolumeRegion * vr, const Point intersectionPoint,
    const Vector &w, float maxDistSquared);

// VolumePhotonIntegrator Local Definitions
inline bool unsuccessful(uint32_t needed, uint32_t found, uint32_t shot) {
    return (found < needed && (found == 0 || found < shot / 1024));
}

inline float kernel(const Photon *photon, const Point &p,
                    float maxDist2) {
    float s = (1.f - DistanceSquared(photon->p, p) / maxDist2);
    return 3.f * INV_PI * s * s;
}

Spectrum LVolumePhoton(KdTree<Photon> *map, int nPaths, int nLookup,
    ClosePhoton *lookupBuf, VolumeRegion * vr, const Point intersectionPoint,
    const Vector &wo, float maxDist2){
	Spectrum L(0.0f);

	PhotonProcess proc(nLookup, lookupBuf);
	map->Lookup(intersectionPoint, proc, maxDist2);

	float maxDist3 = pow(sqrt(maxDist2), 3.0);

	if(proc.nFound == 0) return Spectrum(0.0f);

	float scale = 1.f / (float(nPaths) * maxDist3 * M_PI * 4.f / 3.f);
	
	ClosePhoton * photons = proc.photons;
	int nFound = proc.nFound;

	for(int i = 0; i < nFound; i++){
		const Photon * photon = photons[i].photon;
		L += vr->p(intersectionPoint, photon->wi, wo, 0.1f) * photon->alpha * scale;
	}

	float red[] = {1.0, 0.0, 0.0};
	RGBSpectrum s = RGBSpectrum::FromRGB(red);
	return L;
}

//shoot photons in the volume
class VolumePhotonShootingTask : public Task {
public:
    VolumePhotonShootingTask(int tn, float ti, Mutex &m, VolumePhotonIntegrator *in,
        ProgressReporter &prog, bool &at, int &ndp,
        vector<Photon> &volume, vector<Spectrum> &rpR, vector<Spectrum> &rpT,
        uint32_t &ns, Distribution1D *distrib, const Scene *sc,
        const Renderer *sr, float sz)
    : taskNum(tn), time(ti), mutex(m), integrator(in), progress(prog),
      abortTasks(at), nDirectPaths(ndp),
      volumePhotons(volume), rpReflectances(rpR), rpTransmittances(rpT),
      nshot(ns), lightDistribution(distrib), scene(sc), renderer (sr), stepSize(sz) { }
    void Run();

    int taskNum;
    float time;
    Mutex &mutex;
    VolumePhotonIntegrator *integrator;
    ProgressReporter &progress;
    bool &abortTasks;
    int &nDirectPaths;
    vector<Photon> &volumePhotons;
    vector<Spectrum> &rpReflectances, &rpTransmittances;
    uint32_t &nshot;
    const Distribution1D *lightDistribution;
    const Scene *scene;
    const Renderer *renderer;
	float stepSize;
};


// VolumePhotonIntegrator Method Definitions
VolumePhotonIntegrator::VolumePhotonIntegrator(int nvol,
        int nl, int mdepth, int mphodepth, float mdist, bool fg,
        int gs, float ga, float sz, float cs) {
    nLookup = nl;
    maxSpecularDepth = mdepth;
    maxPhotonDepth = mphodepth;
    maxDistSquared = mdist * mdist;
    finalGather = fg;
    cosGatherAngle = cos(Radians(ga));
    gatherSamples = gs;

	nVolumePhotonsWanted = nvol;
	volumeMap = NULL;
	nVolumePaths = 0;
	stepSize = sz;
	causticScale = cs;
}


VolumePhotonIntegrator::~VolumePhotonIntegrator() {
    delete volumeMap;
}

//Single scattering integrator code from pbrt
void VolumePhotonIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene) {
    tauSampleOffset = sample->Add1D(1);
    scatterSampleOffset = sample->Add1D(1);
}


void VolumePhotonIntegrator::Preprocess(const Scene *scene,
        const Camera *camera, const Renderer *renderer) {
    
	if (scene->lights.size() == 0) return;

	if (!scene->volumeRegion){
	}
	volumeRegion = scene->volumeRegion;
	
    Mutex *mutex = Mutex::Create();
    int nDirectPaths = 0;
	vector<Photon> volumetricPhotons;
	bool abortTasks = false;
    volumetricPhotons.reserve(nVolumePhotonsWanted);
    uint32_t nshot = 0;
    vector<Spectrum> rpReflectances, rpTransmittances;

    Distribution1D *lightDistribution = ComputeLightSamplingCDF(scene);

    ProgressReporter progress(nVolumePhotonsWanted, "Shooting photons");
    vector<Task *> photonShootingTasks;
    int nTasks = NumSystemCores();
    for (int i = 0; i < nTasks; ++i)
        photonShootingTasks.push_back(new VolumePhotonShootingTask(
            i, camera ? camera->shutterOpen : 0.f, *mutex, this, progress, abortTasks, nDirectPaths,
            volumetricPhotons, rpReflectances, rpTransmittances,
            nshot, lightDistribution, scene, renderer, stepSize));
    EnqueueTasks(photonShootingTasks);
    WaitForAllTasks();
    for (uint32_t i = 0; i < photonShootingTasks.size(); ++i)
        delete photonShootingTasks[i];
    Mutex::Destroy(mutex);
    progress.Done();

	if(volumetricPhotons.size() > 0){
		volumeMap = new KdTree<Photon>(volumetricPhotons);
	}
}


void VolumePhotonShootingTask::Run() {
    MemoryArena arena;
    RNG rng(31 * taskNum);
    vector<Photon> localVolumePhotons;
    uint32_t totalPaths = 0;
    bool volumeDone = (integrator->nVolumePhotonsWanted == 0);
    PermutedHalton halton(6, rng);
    vector<Spectrum> localRpReflectances, localRpTransmittances;

	VolumeRegion * volumeRegion = scene->volumeRegion;

    while (true) {
	
        const uint32_t blockSize = 1024;
        for (uint32_t i = 0; i < blockSize; ++i) {
	    float u[6];
            halton.Sample(++totalPaths, u);
            // Choose light to shoot photon from
            float lightPdf;
            int lightNum = lightDistribution->SampleDiscrete(u[0], &lightPdf);
            const Light *light = scene->lights[lightNum];

            // Generate _photonRay_ from light source and initialize _alpha_
            RayDifferential photonRay;
            float pdf;
            LightSample ls(u[1], u[2], u[3]);
            Normal Nl;
            Spectrum Le = light->Sample_L(scene, ls, u[4], u[5],
                                          time, &photonRay, &Nl, &pdf);

            if (pdf == 0.f || Le.IsBlack()) continue;
            Spectrum alpha = (AbsDot(Nl, photonRay.d) * Le) / (pdf * lightPdf);
		
			if (!alpha.IsBlack()) {

				float t0, t1;
				if(!volumeRegion->IntersectP(photonRay, &t0, &t1)) continue;

				Point xs = photonRay(t0);
				bool photonDead = false;
				bool scatteredPhoton = false;

				int numBounces = 0;

				if(t0 > photonRay.mint + RAY_EPSILON){
					photonRay = RayDifferential(photonRay(t0), photonRay.d, photonRay, RAY_EPSILON);
				}

				while(!photonDead && ( numBounces <= integrator->maxPhotonDepth ) && 
						volumeRegion->IntersectP(photonRay, &t0, &t1)){

					photonRay.mint = t0;
					photonRay.maxt = t0;

					photonRay.d = Normalize(photonRay.d);

					float volSigmaT = volumeRegion->sigma_t(photonRay(t0), photonRay.d, 0.1f).y();
					float volSigmaS = volumeRegion->sigma_s(photonRay(t0), photonRay.d, 0.1f).y();
					float albedo = volSigmaS / volSigmaT;
					float randNum = (float)rand()/(float)RAND_MAX;
					
					bool RAY_MARCHING = true;
					float avgDist = 0.0f;
					if(RAY_MARCHING){
						float rt0, rt1;
						float length = photonRay.d.Length();
						float u = randNum;
						float totalT = 0.0f;

						rt0 = t0;
						rt1 = t1;
						
						int count = 0;
						while(rt0 < rt1){
							totalT += volumeRegion->sigma_t(photonRay(rt0), -photonRay.d, 0.1f).y();
							float avgT = totalT/float(count+1);
							avgDist = -1 * log(u)/avgT;
							
							rt0 += stepSize;
							count++;
							
							if(avgDist > rt0)
								break;
						}
					    
						avgDist *= stepSize;
					}
					else{
						avgDist = - log(randNum) / volSigmaT;
					}
					Point interactionPoint = photonRay(t0 + avgDist);
					if(avgDist > (t1 - t0)){
						break;
					}

					if(avgDist < RAY_EPSILON){
						break;
					}
					photonRay.maxt += avgDist;
					Intersection surfaceIntersection;
					
					if(scene->Intersect(photonRay, &surfaceIntersection)){
						Vector diff = photonRay.o - surfaceIntersection.dg.p;
						float surfaceIntersectionTime = diff.Length();
						if(surfaceIntersectionTime < (t0 + avgDist)){
							BSDF * bsdf = surfaceIntersection.GetBSDF(photonRay, arena);							
							Vector wo = -photonRay.d;
					    		Vector wi;
					    		float pdf;
					    		BxDFType flags;
					    		Spectrum fr = bsdf->Sample_f(wo, &wi, BSDFSample(rng), &pdf, BSDF_ALL, &flags);
			
							if (fr.IsBlack() || pdf == 0.f) break;
							
							Spectrum anew = alpha * fr * AbsDot(wi, bsdf->dgShading.nn) / pdf;
							float continueProb = min(1.f, anew.y()/alpha.y());
							if(rng.RandomFloat() > continueProb){
								break;
							}
							alpha = anew / continueProb;
;
							photonRay = RayDifferential(surfaceIntersection.dg.p, wi, photonRay, surfaceIntersection.rayEpsilon);
							numBounces++;
							
							scatteredPhoton = true;
							continue;
						}
					}
					if(scatteredPhoton){
						Photon vp(interactionPoint, alpha, -photonRay.d);
						
						localVolumePhotons.push_back(vp);
					}

									
					randNum = ((float)rand()) / (float)RAND_MAX;
						
					if(randNum <= albedo){
						
						randNum = ((float)rand()) / (float)RAND_MAX;

						float g;
						if(volumeRegion){
							g = 20.0f;
						}else{
						}
						
						RayDifferential scatteredRay;
						float u1 = ((float)rand()) / (float)RAND_MAX;
						float u2 = ((float)rand()) / (float)RAND_MAX;

						Vector newDirection = SampleHG(photonRay.d, g, u1, u2);
						newDirection = Normalize(newDirection);
						scatteredRay.o = interactionPoint;
						scatteredRay.d = newDirection;
						photonRay = RayDifferential(interactionPoint, newDirection, photonRay, RAY_EPSILON);
					}
					else{
						photonDead = true;
						break;
					}
					numBounces++;
					scatteredPhoton = true;

				}			

            }
            arena.FreeAll();
			
        }	MutexLock lock(mutex);

		if(abortTasks)
			return;

		if(nshot > 500000 && unsuccessful(integrator->nVolumePhotonsWanted, volumePhotons.size(), blockSize)){
			Error("Unable to store enough volume photons. Giving up.\n");
			volumePhotons.erase(volumePhotons.begin(), volumePhotons.end());
			abortTasks = true;
			return;
		}
		
		progress.Update(localVolumePhotons.size());
		nshot += blockSize;

		if(!volumeDone){
			integrator->nVolumePaths += blockSize;
			for(uint32_t i = 0; i < localVolumePhotons.size(); ++i){
				volumePhotons.push_back(localVolumePhotons[i]);	
			}
			if(volumePhotons.size() >= integrator->nVolumePhotonsWanted){
				volumeDone = true;
			}
			localVolumePhotons.erase(localVolumePhotons.begin(), localVolumePhotons.end());
		}
        if (volumeDone)
            break;
    
	}
}


Spectrum VolumePhotonIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Sample *sample, RNG &rng, 
		Spectrum * T, MemoryArena &arena) const {

	VolumeRegion *vr = scene->volumeRegion;
    float t0, t1;
    if (!vr || !vr->IntersectP(ray, &t0, &t1) || (t1-t0) == 0.f) {
        *T = 1.f;
        return 0.f;
    }

    Spectrum Lv(0.);

    int nSamples = Ceil2Int((t1-t0) / stepSize);
    float step = (t1 - t0) / nSamples;
    Spectrum Tr(1.f);
    Point p = ray(t0), pPrev;
    Vector w = -ray.d;
    t0 += sample->oneD[scatterSampleOffset][0] * step;

    float *lightNum = arena.Alloc<float>(nSamples);
    LDShuffleScrambled1D(1, nSamples, lightNum, rng);
    float *lightComp = arena.Alloc<float>(nSamples);
    LDShuffleScrambled1D(1, nSamples, lightComp, rng);
    float *lightPos = arena.Alloc<float>(2*nSamples);
    LDShuffleScrambled2D(1, nSamples, lightPos, rng);
    uint32_t sampOffset = 0;
    for (int i = 0; i < nSamples; ++i, t0 += step) {
        pPrev = p;
        p = ray(t0);
        Ray tauRay(pPrev, p - pPrev, 0.f, 1.f, ray.time, ray.depth);
        Spectrum stepTau = vr->tau(tauRay,
                                   .5f * stepSize, rng.RandomFloat());
        Tr *= Exp(-stepTau);

        if (Tr.y() < 1e-3) {
            const float continueProb = .5f;
            if (rng.RandomFloat() > continueProb) break;
            Tr /= continueProb;
        }

        Lv += Tr * vr->Lve(p, w, ray.time);
        Spectrum ss = vr->sigma_s(p, w, ray.time);
        if (!ss.IsBlack() && scene->lights.size() > 0) {
            int nLights = scene->lights.size();
            int ln = min(Floor2Int(lightNum[sampOffset] * nLights),
                         nLights-1);
            Light *light = scene->lights[ln];
            float pdf;
            VisibilityTester vis;
            Vector wo;
            LightSample ls(lightComp[sampOffset], lightPos[2*sampOffset],
                           lightPos[2*sampOffset+1]);
            Spectrum L = light->Sample_L(p, 0.f, ls, ray.time, &wo, &pdf, &vis);
            
            if (!L.IsBlack() && pdf > 0.f && vis.Unoccluded(scene)) {
                Spectrum Ld = L * vis.Transmittance(scene, renderer, NULL, rng, arena);
                Lv += Tr * ss * vr->p(p, w, -wo, ray.time) * Ld * float(nLights) /
                        pdf;
            }

			ClosePhoton *lookupBuf = arena.Alloc<ClosePhoton>(nLookup);	
			Lv += LVolumePhoton(volumeMap, nVolumePaths, nLookup, lookupBuf, vr, p, w, maxDistSquared);
        }
        ++sampOffset;
    }

    *T = Tr;

    return Lv * step * causticScale;
}

Spectrum VolumePhotonIntegrator::Transmittance(const Scene * scene, const Renderer * renderer, 
		const RayDifferential & ray, const Sample * sample, RNG & rng, MemoryArena & arena) const{

    if (!scene->volumeRegion) return Spectrum(1.f);
    float step, offset;
    if (sample) {
        step = stepSize;
        offset = sample->oneD[tauSampleOffset][0];
    }
    else {
        step = 4.f * stepSize;
        offset = rng.RandomFloat();
    }
    Spectrum tau = scene->volumeRegion->tau(ray, step, offset);
    return Exp(-tau);
}

VolumePhotonIntegrator *CreateVolumePhotonMapSurfaceIntegrator(const ParamSet &params) {
    int nVolume = params.FindOneInt("volumephotons", VOLUME_PHOTONS);
    int nUsed = params.FindOneInt("nused", VOLUME_LOOKUP);
    if (PbrtOptions.quickRender) nVolume = nVolume / 10;
    if (PbrtOptions.quickRender) nUsed = max(1, nUsed / 10);
    int maxSpecularDepth = params.FindOneInt("maxspeculardepth", 5);
    int maxPhotonDepth = params.FindOneInt("maxphotondepth", 5);
    bool finalGather = params.FindOneBool("finalgather", true);
    int gatherSamples = params.FindOneInt("finalgathersamples", 32);
    if (PbrtOptions.quickRender) gatherSamples = max(1, gatherSamples / 4);
    float maxDist = params.FindOneFloat("maxdist", .5f);
    float gatherAngle = params.FindOneFloat("gatherangle", 10.f);
	float causticScale = params.FindOneFloat("scalefactor", CAUSTIC_SCALE);

	float stepSize = params.FindOneFloat("stepsize", 0.1f);

    return new VolumePhotonIntegrator(nVolume, nUsed, maxSpecularDepth, 
			maxPhotonDepth, maxDist, finalGather, gatherSamples, gatherAngle, stepSize, causticScale);
}


#define FILENAME "spotfog.pbrt"
void VolumePhotonIntegrator::WritePhotonsOutToMap(vector<Photon>* map)
{
	std::ofstream file;
	

	file.open("debug_file.pbrt", std::ios::trunc);

	file << "AttributeBegin\n";
	file << "\tMaterial \"matte\" \"color Kd\" [1.0 1.0 1.0]\n";

	Point position = Point(0,0,0);
	for(unsigned int a = 0; a < map->size(); a++)
	{
		Photon photon = (*map)[a];
		if(isinf(photon.p.x) || isinf(photon.p.y) || isinf(photon.p.z))
		{
			continue;
		}
		Vector difference = photon.p - position;
		file << "Translate " << difference.x << " " << difference.y << " " << difference.z << "\n";
		file << "Shape \"sphere\" \"float radius\" [0.002]\n";
		position = photon.p;
	}
	file << "AttributeEnd\n";
	file << "WorldEnd\n";
	file.close();
	printf("file written\n");
}

