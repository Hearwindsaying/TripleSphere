// STL
#include <random>
#include <vector>

// TripleSphere
#include <TripleSphere.h>

// Utils
#include <VecType.h>

using namespace ts;


std::vector<double> Sphere_MonteCarlo_Sl_Clipped(int l, float *SinSigma, float *CosSigma, float *SinW, float *CosW, float3 *L_world, float3 *outwDir_world)
{
    /*int iteration = 20000000;*/
    int iteration = 10000000;

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    float3 wDir_world = uniformSamplingHemisphere(dis(gen), dis(gen));
    assert(wDir_world.z >= 0.f && abs(length(wDir_world) - 1.f) <= 1e-6f);

    *outwDir_world = wDir_world;

    // Setup a local frame
    float3 localX, localY, localZ;
    localZ = wDir_world;
    BSDFMath::CoordinateSystem(localZ, &localX, &localY);

    // Cone direction in world frame
    float3 coneZ_world = uniformSamplingHemisphere(dis(gen), dis(gen));
    *L_world = coneZ_world;
    float3 coneX_world, coneY_world;
    BSDFMath::CoordinateSystem(coneZ_world, &coneX_world, &coneY_world);

    float coneSigma = dis(gen) * M_PI / 2.0;
    float coneCosMax = cos(coneSigma); *CosSigma = cos(coneSigma); *SinSigma = sin(coneSigma);

    // Cone direction in local frame
    float3 coneX_local = BSDFMath::WorldToLocal(coneX_world, localX, localY, localZ);
    float3 coneY_local = BSDFMath::WorldToLocal(coneY_world, localX, localY, localZ);
    float3 coneZ_local = BSDFMath::WorldToLocal(coneZ_world, localX, localY, localZ);

    *CosW = coneZ_local.z;
    *SinW = sqrt(1.f - *CosW * *CosW);

    std::vector<double> resultGenerals(l + 1, 0.0);
    for (int i = 0; i < iteration; ++i)
    {
        // for each band:
        for (int j = 0; j <= l; ++j)
        {
            float u = dis(gen), v = dis(gen);

            float3 wi_world = UniformSampleCone(u, v, coneCosMax, coneX_world, coneY_world, coneZ_world);
            if (wi_world.z >= 0.f)
                resultGenerals[j] += sqrtf((2.f * j + 1) / (4.f*M_PI)) * LegendreP(j, dot(wDir_world, wi_world));
        }
    }
    for (auto& resultGeneral : resultGenerals)
    {
        resultGeneral /= (iteration * UniformConePdf(coneCosMax));
    }

    return resultGenerals;
}

/**
 * @brief Generating random Light directions and Lobe directions to perform
 * Monte Carlo computing Sl, accounting for clipping case and
 * fully hemispherical cap case. Could be used to validate Sphere_computeSl()
 */
double MonteCarlo_Sl_Clipped(int l, float *SinSigma, float *CosSigma, float *SinW, float *CosW, float3 *L_world, float3 *outwDir_world)
{
    /*int iteration = 20000000;*/
    int iteration = 5000000;

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<> dis(0.0, 1.0);

    //float3 wDir_world = uniformSamplingHemisphere(dis(gen), dis(gen));
    float3 wDir_world = UniformSampleSphere(dis(gen), dis(gen));

    *outwDir_world = wDir_world;

    // Setup a local frame
    float3 localX, localY, localZ;
    localZ = wDir_world;
    BSDFMath::CoordinateSystem(localZ, &localX, &localY);

    // Cone direction in world frame
    float3 coneZ_world = UniformSampleSphere(dis(gen), dis(gen));
    *L_world = coneZ_world;
    float3 coneX_world, coneY_world;
    BSDFMath::CoordinateSystem(coneZ_world, &coneX_world, &coneY_world);

    float coneSigma = dis(gen) * M_PI / 2.0;
    float coneCosMax = cos(coneSigma); *CosSigma = cos(coneSigma); *SinSigma = sin(coneSigma);

    // Cone direction in local frame
    float3 coneX_local = BSDFMath::WorldToLocal(coneX_world, localX, localY, localZ);
    float3 coneY_local = BSDFMath::WorldToLocal(coneY_world, localX, localY, localZ);
    float3 coneZ_local = BSDFMath::WorldToLocal(coneZ_world, localX, localY, localZ);

    *CosW = coneZ_local.z;
    *SinW = sqrt(1.f - *CosW * *CosW);

    double resultGeneral = 0.f;
    for (int i = 0; i < iteration; ++i)
    {
        {
            float u = dis(gen), v = dis(gen);

            float3 wi_world = UniformSampleCone(u, v, coneCosMax, coneX_world, coneY_world, coneZ_world);
            if (wi_world.z >= 0.f)
                resultGeneral += LegendreP(l, dot(wDir_world, wi_world));
        }
    }
    resultGeneral /= (iteration * UniformConePdf(coneCosMax));
    return resultGeneral;
}

void TestSl(float epsilon = 1e-6f)
{
    int nfails = 0;
    int ntests = 0;
    constexpr int lmax = 9;

    // Test All Sl with clipping
#if 0
    {
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> disSigma(0.0, 2.0*M_PI);

        for (int subtest = 0; subtest < 50; ++subtest)
        {
            float sinSigma;
            float cosSigma;
            float sinW;
            float cosW;

            float3 L_world;
            float3 wDir_world;

            double mc = MonteCarlo_Sl_Clipped(lmax, &sinSigma, &cosSigma, &sinW, &cosW, &L_world, &wDir_world);
            float rec[lmax + 1][2 * lmax + 1];
            Sphere_computeLw<lmax>(sinSigma, cosSigma, L_world, wDir_world, 0, rec); // Only use 0th result to match 'Sphere_MonteCarlo_Sl_Clipped'

            if (std::abs(mc - rec[lmax][0]) >= epsilon)
            {
                ++nfails;
                printf("\033[0;31m");
                printf("[%d|FAILED!]Test failed at Line:%d mc:%f Sl Recurrence:%f\n", ntests, __LINE__, mc, rec[lmax][0]);
                printf("\033[0m");
            }
            else
            {
                printf("\033[0;32m");
                printf("[%d|PASSED!]\n", ntests);
                printf("\033[0m");
            }
            ++ntests;
        }
    }
#endif

    {
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> disSigma(0.0, 2.0*M_PI);
 
        for (int subtest = 0; subtest < 50; ++subtest)
        {
            float sinSigma;
            float cosSigma;
            float sinW;
            float cosW;

            float3 L_world;
            float3 wDir_world;

            std::vector<double> mc = Sphere_MonteCarlo_Sl_Clipped(lmax, &sinSigma, &cosSigma, &sinW, &cosW, &L_world, &wDir_world);
            float rec[lmax + 1][2 * lmax + 1];
            Sphere_computeLw<lmax>(sinSigma, cosSigma, L_world, wDir_world, 0, rec);

            bool testFailed = false;
            for (int i = 0; i <= lmax; ++i)
            {
                if (std::abs(mc[i] - rec[i][0]) >= epsilon)
                {
                    testFailed = true;
                    printf("\033[0;31m");
                    printf("[%d|FAILED!]bandNum=%d, mc:%f recurrence:%f\n", ntests, i, mc[i], rec[i][0]);
                    printf("\033[0m");
                }
            }
            if (testFailed)
            {
                ++nfails;
            }
            else
            {
                printf("\033[0;32m");
                printf("[%d|PASSED!]\n", ntests);
                printf("\033[0m");
            }
            ++ntests;
        }
    }

    printf("\nLegendre Sl Test coverage:%f%%(%d/%d) passed!\n", 100.f*static_cast<float>(ntests - nfails) / ntests, ntests - nfails, ntests);
}

int main()
{
    TestSl();

    return 0;
}