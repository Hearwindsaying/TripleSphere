// STL
#include <random>
#include <vector>

// TripleSphere
#include <TripleSphere.h>

// Utils
#include <VecType.h>

using namespace ts;

/**
 * @brief Generating random Light directions and Lobe directions to perform
 * Monte Carlo computing Sl, accounting for clipping case and
 * fully hemispherical cap case. Could be used to validate Sphere_computeSl()
 */
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
 * @brief Analytical boundary integral $B_l$ exported from Mathematica.
 */
float Bx_math(int x, float SinSigma, float SinW, float CosSigma, float CosW)
{
    float a1 = CosSigma * CosW;
    float b1 = -SinSigma * SinW;

    assert(x >= 0 && x <= 9);
    switch (x)
    {
    case 0: return 2 * (-1 + std::powf(CosSigma, 2))*CosW*M_PI;
    case 1: return (std::powf(b1, 2)*CosSigma + 2 * a1*(-1 + std::powf(CosSigma, 2))*CosW)*M_PI;
    case 2: return ((6 * a1*std::powf(b1, 2)*CosSigma + (-2 + 6 * std::powf(a1, 2) + 3 * std::powf(b1, 2))*(-1 + std::powf(CosSigma, 2))*CosW)*M_PI) / 2.;
    case 3: return ((3 * std::powf(b1, 2)*(-4 + 20 * std::powf(a1, 2) + 5 * std::powf(b1, 2))*CosSigma +
        4 * a1*(-6 + 10 * std::powf(a1, 2) + 15 * std::powf(b1, 2))*(-1 + std::powf(CosSigma, 2))*CosW)*M_PI) / 8.;
    case 4: return ((20 * a1*std::powf(b1, 2)*(-12 + 28 * std::powf(a1, 2) + 21 * std::powf(b1, 2))*CosSigma +
        (280 * std::powf(a1, 4) + 120 * std::powf(a1, 2)*(-2 + 7 * std::powf(b1, 2)) + 3 * (8 - 40 * std::powf(b1, 2) + 35 * std::powf(b1, 4)))*(-1 + std::powf(CosSigma, 2))*CosW)*M_PI) / 32.;
    case 5: return ((15 * std::powf(b1, 2)*(8 + 168 * std::powf(a1, 4) - 28 * std::powf(b1, 2) + 21 * std::powf(b1, 4) + 28 * std::powf(a1, 2)*(-4 + 9 * std::powf(b1, 2)))*CosSigma +
        2 * a1*(504 * std::powf(a1, 4) + 280 * std::powf(a1, 2)*(-2 + 9 * std::powf(b1, 2)) + 15 * (8 - 56 * std::powf(b1, 2) + 63 * std::powf(b1, 4)))*(-1 + std::powf(CosSigma, 2))*CosW)*
        M_PI) / 64.;
    case 6: return ((42 * a1*std::powf(b1, 2)*(264 * std::powf(a1, 4) + 60 * std::powf(a1, 2)*(-4 + 11 * std::powf(b1, 2)) + 5 * (8 - 36 * std::powf(b1, 2) + 33 * std::powf(b1, 4)))*CosSigma +
        (-80 + 840 * std::powf(b1, 2) + 21 * (176 * std::powf(a1, 6) - 90 * std::powf(b1, 4) + 55 * std::powf(b1, 6) + 120 * std::powf(a1, 4)*(-2 + 11 * std::powf(b1, 2)) +
            10 * std::powf(a1, 2)*(8 - 72 * std::powf(b1, 2) + 99 * std::powf(b1, 4))))*(-1 + std::powf(CosSigma, 2))*CosW)*M_PI) / 128.;
    case 7: return ((7 * std::powf(b1, 2)*(64 * (-5 + 135 * std::powf(a1, 2) - 495 * std::powf(a1, 4) + 429 * std::powf(a1, 6)) + 720 * (3 - 66 * std::powf(a1, 2) + 143 * std::powf(a1, 4))*std::powf(b1, 2) +
        3960 * (-1 + 13 * std::powf(a1, 2))*std::powf(b1, 4) + 2145 * std::powf(b1, 6))*CosSigma +
        8 * a1*(6864 * std::powf(a1, 6) + 5544 * std::powf(a1, 4)*(-2 + 13 * std::powf(b1, 2)) + 630 * std::powf(a1, 2)*(8 - 88 * std::powf(b1, 2) + 143 * std::powf(b1, 4)) +
            35 * (-16 + 216 * std::powf(b1, 2) - 594 * std::powf(b1, 4) + 429 * std::powf(b1, 6)))*(-1 + std::powf(CosSigma, 2))*CosW)*M_PI) / 1024.;
    case 8: return (51480 * std::powf(a1, 7)*std::powf(b1, 2)*CosSigma*M_PI + 18018 * std::powf(a1, 5)*std::powf(b1, 2)*(-4 + 15 * std::powf(b1, 2))*CosSigma*M_PI +
        3465 * std::powf(a1, 3)*std::powf(b1, 2)*(8 - 52 * std::powf(b1, 2) + 65 * std::powf(b1, 4))*CosSigma*M_PI +
        (315 * a1*std::powf(b1, 2)*(-64 + 528 * std::powf(b1, 2) - 1144 * std::powf(b1, 4) + 715 * std::powf(b1, 6))*CosSigma*M_PI) / 8. +
        12870 * std::powf(a1, 8)*(-1 + std::powf(CosSigma, 2))*CosW*M_PI + 12012 * std::powf(a1, 6)*(-2 + 15 * std::powf(b1, 2))*(-1 + std::powf(CosSigma, 2))*CosW*M_PI +
        (3465 * std::powf(a1, 4)*(8 - 104 * std::powf(b1, 2) + 195 * std::powf(b1, 4))*(-1 + std::powf(CosSigma, 2))*CosW*M_PI) / 2. +
        (315 * std::powf(a1, 2)*(-16 + 264 * std::powf(b1, 2) - 858 * std::powf(b1, 4) + 715 * std::powf(b1, 6))*(-1 + std::powf(CosSigma, 2))*CosW*M_PI) / 2. +
        (35 * (128 - 2304 * std::powf(b1, 2) + 9504 * std::powf(b1, 4) - 13728 * std::powf(b1, 6) + 6435 * std::powf(b1, 8))*(-1 + std::powf(CosSigma, 2))*CosW*M_PI) / 64.) / 128.;
    case 9: return ((45 * std::powf(b1, 2)*(896 + 1408 * std::powf(a1, 2)*(-28 + 182 * std::powf(a1, 2) - 364 * std::powf(a1, 4) + 221 * std::powf(a1, 6)) +
        9856 * (-1 + 39 * std::powf(a1, 2) - 195 * std::powf(a1, 4) + 221 * std::powf(a1, 6))*std::powf(b1, 2) +
        32032 * (1 - 30 * std::powf(a1, 2) + 85 * std::powf(a1, 4))*std::powf(b1, 4) + 40040 * (-1 + 17 * std::powf(a1, 2))*std::powf(b1, 6) + 17017 * std::powf(b1, 8))*CosSigma +
        2 * a1*(128 * (315 - 4620 * std::powf(a1, 2) + 143 * std::powf(a1, 4)*(126 - 180 * std::powf(a1, 2) + 85 * std::powf(a1, 4))) +
            126720 * (-7 + 91 * std::powf(a1, 2) - 273 * std::powf(a1, 4) + 221 * std::powf(a1, 6))*std::powf(b1, 2) +
            4324320 * (1 - 10 * std::powf(a1, 2) + 17 * std::powf(a1, 4))*std::powf(b1, 4) + 2402400 * (-3 + 17 * std::powf(a1, 2))*std::powf(b1, 6) + 3828825 * std::powf(b1, 8))*
            (-1 + std::powf(CosSigma, 2))*CosW)*M_PI) / 16384.;
    }

}

/**
 * @brief Recurrence formulae evaluation for $B_l$.
 * Note that 'sinSigma,sinW,cosSigma,cosW' are purely parameters for
 * the integral. We do not classify three cases as we do for $S_l$.
 */
std::vector<float> Bl(int lmax, float sinSigma, float sinW, float cosSigma, float cosW)
{
    //assert(abs(a - cosSigma * cosW) < 1e-6f);
    std::vector<float> Bls; Bls.reserve(lmax + 1);

    float a = cosSigma * cosW;
    float b = -sinSigma * sinW;

    float E0 = 0.f;
    float E1 = 2.f*M_PI;
    float F0 = 0.f;
    float F1 = 2.f*M_PI*a;
    float C0 = 2.f*M_PI;
    float C1 = 2.f*M_PI*a;
    float D0 = 2.f*M_PI*a;
    float D1 = (2.f*a*a + b * b)*M_PI;

    float CO; // undefined state
    float Df1, Df2;

    float B0 = cosSigma * D0 - cosW * C0;
    float B1 = cosSigma * D1 - cosW * C1;
    Bls.push_back(B0);
    Bls.push_back(B1);
    for (int l = 2; l <= lmax; ++l)
    {
        if (l > 2)
        {
            float l_1 = l - 1;
            C1 = (2.f*l_1 - 1.f) / l_1 * D0 - (l_1 - 1.f) / l_1 * CO;
        }
        float E2 = (l % 2 == 0 ? ((2.f*l - 1.f) * C1 + E0) : ((2.f*l - 1.f) * C1 + E1));

        float F2; // undefined state
        if (l > 3)
        {
            float l_2 = l - 2;
            F2 = (l % 2 == 0 ? (2.f*l_2 - 1.f)*Df1 + F0 : (2.f*l_2 - 1.f)*Df2 + F1);
        }

        float Dl;
        if (l == 2)
        {
            Dl = (a*(2.f*l - 1.f) / l + 2.f*a*(2.f*l - 1.f))*D1 - a * (l - 1.f) / l * C0 +/*LegendreP<2>()*/+(b*b - a * a - 1)*E2 + l * C1 + 2.f*a*F0;
        }
        else if (l == 3)
        {
            Dl = (a*(2.f*l - 1.f) / l + 2.f*a*(2.f*l - 1.f))*D1 - a * (l - 1.f) / l * C0 +/*LegendreP<2>()*/+(b*b - a * a - 1)*E2 + l * C1 + 2.f*a*F1;
        }
        else
        {
            Dl = (a*(2.f*l - 1.f) / l + 2.f*a*(2.f*l - 1.f))*D1 - a * (l - 1.f) / l * C0 +/*LegendreP<2>()*/+(b*b - a * a - 1)*E2 + l * C1 + 2.f*a*F2;
        }

        Dl /= l + 1;

        float Cl = (2.f*l - 1.f) / l * D1 - (l - 1.f) / l * C0;
        float Bl = cosSigma * Dl - cosW * Cl;
        Bls.push_back(Bl);

        if (l % 2 == 0)
        {
            Df1 = D1;
            E0 = E2;
        }
        else
        {
            Df2 = D1;
            E1 = E2;
        }
        D0 = D1; D1 = Dl;
        CO = C0; C0 = C1;
        if (l > 3)
        {
            if (l % 2 == 0)
            {
                F0 = F2;
            }
            else
            {
                F1 = F2;
            }
        }
    }
    return Bls;
}


/**
 * @brief Test Surface Integral $S_l$ up to order 10 in the paper.
 */
void TestSl(float epsilon = 1e-3f)
{
    int nfails = 0;
    int ntests = 0;
    constexpr int lmax = 9;

    // Test All Sl with clipping
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

/**
 * @brief Test Boundary Integral $B_l$ up to order 10 in the paper.
 */
void TestBl(float epsilon = 1e-6f)
{
    int nfails = 0;
    int ntests = 0;

    // Test All,[0,9]
#if 1
    {
        // Test D3
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> disSigma(0.0, 2.0*M_PI);

        printf("Test B0~B9.\n");
        for (int subtest = 0; subtest < 100; ++subtest)
        {
            float sigma = disSigma(gen), w = disSigma(gen);
            float sinSigma = sin(sigma);
            float cosSigma = cos(sigma);
            float sinW = sin(w);
            float cosW = cos(w);

            std::vector<float> mathematica;
            for (int dx = 0; dx < 10; ++dx)
            {
                mathematica.push_back(Bx_math(dx, sinSigma, sinW, cosSigma, cosW));
            }

            std::vector<float> Bx_rec = Bl(9, sinSigma, sinW, cosSigma, cosW);
            for (int dx = 0; dx < 10; ++dx)
            {
                if (std::abs(Bx_rec[dx] - mathematica[dx] >= epsilon))
                {
                    ++nfails;
                    printf("%d:mathematica=%f, Dx recurrence=%f\n", dx, mathematica[dx], Bx_rec[dx]);
                }
                ++ntests;
            }
        }
    }
#endif

    printf("\nLegendre Polynomial Test coverage:%f%%(%d/%d) passed!\n", 100.f*static_cast<float>(ntests - nfails) / ntests, ntests - nfails, ntests);
}

int main()
{
    TestBl();
    TestSl();

    return 0;
}