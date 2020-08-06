#pragma once

// STL
#include <cmath>
#include <cassert>
#include <stdexcept>

// Vector type and helper sampling functions
#include <VecType.h>
#include <SHUtils.h>

/* Provide a reference implementation of 
 * 'Spherical Light Integration over Spherical Caps via Spherical Harmonics'.
 * 
 * To use this header to compute light SH coefficients on the fly, you do not
 * need to access functions inside 'tsImpl', just passing your parms to
 * 'Sphere_computeCoeff'.
 * 
 * Note that for higher performance, we specialize our formulae for three cases
 * and thus have more simpler evaluation formulae. Interestingly, for the unclipped
 * case, we could even avoid the evaluation of Legendre polynomials, which involves
 * several 'powf' (computational expensive).
 *
 *    <'ts/tsImpl'>
 *    + 'Sphere_computeS1_clipped' returns initial condition for our recurrence
 *    formulae, which is used for the clipped case (phi0,phi1,phia,phib). For
 *    unclipped case (and clipped case S0), initial conditions are inside 'Sphere_Lw'
 *    as they have simpler forms.
 *
 *    + 'Sphere_Lw' returns sphere light coefficients up to $lmax$, designed for
 *    unclipped case.
 *    
 *    + 'Sphere_Lw_clipped' is the clipped case version of 'Sphere_Lw'.
 *    
 *    <'ts'>
 *    + 'Sphere_computeLw' identifies the case:
 *    unclipped/clipped/fully clipped and forward to corresponding functions
 *    to compute SH light coefficients. 
 *
 *    + 'Sphere_computeCoeff' is for client, applying Zonal Harmonics Factorization 
 *    and return the expected SH light coefficients for sphere lights.
 *
 * Feel free to contact hearwindsaying@gmail.com if you have any questions.
 */

namespace ts
{
    /************************************************************************/
    /*                         Inner Implmentation Code                     */
    /************************************************************************/
    namespace tsImpl
    {
        /**
          * @brief Compute several constants and give phia, phib, phi0, phi1.
          * Need to promise that it it's a clipped case! capCase=2
          * @param[in] wDir_world        Lobe direction in TBN space
          * @param[in] L_world           Light direction in TBN space
          * @param[in] sinSigma,cosSigma Spherical cap aperture angle
          * @param[in] sinW,cosW         direction between Lobe direction and Light direction.
          *
          * @param[out] outS1            Clipped case S1 constant
          * @param[out] outPhi0,outPhi1,outPhia,outPihb
          */
          //template<typename float3>
        void Sphere_computeS1_clipped(const float3 &wDir_world, const float3 &L_world,
            const float sinSigma, const float cosSigma,
            const float sinW, const float cosW,
            const float sinWprime, const float cosWprime,
            float *outS1, float *outPhi0, float *outPhi1, float *outPhia, float *outPhib)
        {
            auto phiToWi_Local = [](const float phi, const float sinSigma, const float cosSigma, const float sinW, const float cosW)->float3
            {
                return make_float3(cosSigma*sinW + sinSigma * (cosW*cosf(phi)),
                    sinSigma*sinf(phi),
                    cosSigma*cosW + sinSigma * cosf(phi)*(-sinW));
            };

            // Setup local frame
            float3 LocalZ = safe_normalize(wDir_world);
            float3 LocalY = safe_normalize(cross(LocalZ, L_world));
            float3 LocalX = safe_normalize(cross(LocalY, LocalZ));
            if (length(LocalY) < 0.1f)
                BSDFMath::CoordinateSystem(LocalZ, &LocalX, &LocalY);

            float aPrime = LocalX.z*sinSigma*cosW - LocalZ.z*sinSigma*sinW;
            float bPrime = LocalY.z*sinSigma;
            float cPrime = sign(aPrime)*std::sqrt(aPrime*aPrime + bPrime * bPrime);
            float phi = std::atanf(-bPrime / aPrime);
            float acosTerm = (-LocalX.z*cosSigma*sinW - LocalZ.z*cosSigma*cosW) / cPrime;
            acosTerm = clamp(acosTerm, -1.f, 1.f);

            // Part I
            float phi0 = std::acos(acosTerm) - phi;
            float phi1 = 2.f*M_PI - std::acos(acosTerm) - phi;

            // Sort
            if (!(phi0 < phi1))
                std::swap(phi0, phi1);

            // Updated wi_phi0
            float3 wi_phi0_local = phiToWi_Local(phi0, sinSigma, cosSigma, sinW, cosW);
            float3 wi_phi1_local = phiToWi_Local(phi1, sinSigma, cosSigma, sinW, cosW);

            float3 wi_phi0_world = BSDFMath::LocalToWorld(wi_phi0_local, LocalX, LocalY, LocalZ);
            float3 wi_phi1_world = BSDFMath::LocalToWorld(wi_phi1_local, LocalX, LocalY, LocalZ);

            if (std::abs(L_world.z) <= 1e-4f)
            {
                float3 crossValidate = safe_normalize(cross(wi_phi0_world, safe_normalize(wi_phi1_world - wi_phi0_world)));
                if (crossValidate.z >= 0)
                {
                    std::swap(phi0, phi1);
                }
            }
            else if (L_world.z < 0.f)
            {
                // We need inner arc in this case: <PI
                if (!(phi1 - phi0 < M_PI))
                {
                    phi1 = -(2.f*M_PI - phi1);
                    std::swap(phi0, phi1);
                }
            }
            else
            {
                // We need great arc in this case: >PI
                if (!((phi1 - phi0) > M_PI))
                {
                    phi0 = phi0 + 2.f*M_PI;
                    std::swap(phi0, phi1);
                }
            }

            // Part II
            // LocalZ = safe_normalize(wDir_world);
            LocalY = safe_normalize(cross(LocalZ, make_float3(0, 0, 1.0f)));
            LocalX = safe_normalize(cross(LocalY, LocalZ));

            if (std::abs(LocalZ.z - 1.f) < 1e-4f)
            {
                BSDFMath::CoordinateSystem(LocalZ, &LocalX, &LocalY);
            }

            float3 OL1 = BSDFMath::LocalToWorld(make_float3(cosWprime, 0, -sinWprime), LocalX, LocalY, LocalZ);
            float3 OL2 = BSDFMath::LocalToWorld(make_float3(0, 1, 0), LocalX, LocalY, LocalZ);
            const float3 &OA = wi_phi0_world;
            const float3 &OB = wi_phi1_world;

            float phia = std::acos(dot(OL1, OA));
            float phib = std::acos(dot(OL1, OB));

            if (dot(OA, OL2) < 0.f)
                phia = 2.f*M_PI - phia;

            if (dot(OB, OL2) < 0.f)
                phib = 2.f*M_PI - phib;

            // Sort
            if (!(phia < phib))
                std::swap(phia, phib);
            if (!(phib - phia < M_PI))
            {
                phib = -(2.f*M_PI - phib);
                std::swap(phia, phib);
            }

            *outS1 = 0.5*(sinSigma*cosSigma*sinW*(sinf(phi1) - sinf(phi0)) + sinSigma * sinSigma*cosW*(phi1 - phi0) +
                cosWprime * (phib - phia));
            *outPhi0 = phi0;
            *outPhi1 = phi1;
            *outPhia = phia;
            *outPhib = phib;
        }


        /**
          * @brief Compute Sl Coeffs up to $lmax$ for fully hemispherical cap case.
          * @param[in]  lmax              maximum band
          * @param[in]  sinSigma,cosSigma Spherical cap aperture angle
          * @param[in]  sinW,cosW         direction between Lobe direction and Light direction.
          * @param[in]  lobe_i            lobe direction index
          * @param[out] Lw
          */
        template<int lmax>
        void Sphere_Lw(float sinSigma, float sinW, float cosSigma, float cosW, int lobe_i, float Lw[lmax + 1][2 * lmax + 1])
        {
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

            float B1 = cosSigma * D1 - cosW * C1;

            float S0 = 2.f * M_PI * (1.f - cosSigma);
            float S1 = M_PI * sinSigma * sinSigma * cosW;

            Lw[0][lobe_i] = sqrtf(1.f / (4.f*M_PI)) * S0;
            Lw[1][lobe_i] = sqrtf(3.f / (4.f*M_PI)) * S1;
            // TODO: Redundant computation: I mistakenly compute Bl, where we only need Bl_1 actually.
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

                float Dl = ((a*(2.f*l - 1.f) / l + 2.f*a*(2.f*l - 1.f))*D1 - a * (l - 1.f) / l * C0 + (b*b - a * a - 1)*E2 + l * C1 + 2.f*a
                    *(l == 2 ? F0 : (l == 3 ? F1 : F2))) / (l + 1);

                float Cl = (2.f*l - 1.f) / l * D1 - (l - 1.f) / l * C0;
                float Bl = cosSigma * Dl - cosW * Cl;
                // Bug fixed: note that if we integrate from 0 to 2Pi, we are getting negative orientation of the boundary curve, which leads to negative result compared to Surface Integral Sl. We flipped the sign for Bl here in the final step to fix the issue.
                float Sl = (l - 2.f)*(l - 1.f) / (l*(l + 1.f))* (l % 2 == 0 ? S0 : S1) - (2.f*l - 1.f) / (l*(l + 1.f))*B1;

                Lw[l][lobe_i] = sqrtf((2.f * l + 1) / (4.f*M_PI)) * Sl;

                B1 = Bl;
                if (l % 2 == 0)
                {
                    S0 = Sl;
                }
                else
                {
                    S1 = Sl;
                }

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
        }

        /**
          * @brief Compute Sl Coeffs up to $lmax$ for clipped hemispherical cap case.
          * @param[in]  sinSigma,cosSigma   Spherical cap aperture angle
          * @param[in]  sinW,cosW           direction between Lobe direction and Light direction (both in the same space).
          * @param[in]  phi0,phi1,phia,phib Boundary constants output from Sphere_computeS1_clipped
          * @param[in]  _cosW,_sinW         direction between (0,0,1) and Lobe direction in TBN space.
          * @param[in]  inS0,inS1           clipped case constants S0,S1
          * @param[in]  lobe_i            lobe direction index
          *
          * @param[out] Lw
          */
        template<int lmax>
        void Sphere_Lw_clipped(
            float sinSigma, float sinW, float cosSigma, float cosW,
            float phi0, float phi1, float phia, float phib,
            // PartII params:cosw'
            float _cosW, float _sinW,
            // S0 S1 constants
            float inS0, float inS1,
            int lobe_i,
            float Lw[lmax + 1][2 * lmax + 1])
        {
            /************************************************************************/
            /*                                Part I                                */
            /************************************************************************/
            float a = cosSigma * cosW;
            float b = -sinSigma * sinW;

            /************************************************************************/
            /*                                Part II                               */
            /************************************************************************/
            float _b = -_sinW;

            /************************************************************************/
            /*                                Part I                                */
            /************************************************************************/
            float E0 = 0.f;
            float E1 = phi1 - phi0;
            float F0 = 0.f;
            float F1 = a * (phi1 - phi0) + b * (sinf(phi1) - sinf(phi0));
            float C0 = phi1 - phi0;
            float C1 = a * (phi1 - phi0) + b * (sinf(phi1) - sinf(phi0));
            float D0 = a * (phi1 - phi0) + b * (sinf(phi1) - sinf(phi0));
            float D1 = ((-2 * a*a - b * b)*(phi0 - phi1) - b * (4 * a + b * cosf(phi0))*sinf(phi0) + b * (4 * a + b * cosf(phi1))*sinf(phi1)) / 2;

            /************************************************************************/
            /*                                Part II                               */
            /************************************************************************/
            float _E0 = 0.f;
            float _E1 = phib - phia;
            float _C0 = phib - phia;
            float _C1 = _b * (sinf(phib) - sinf(phia));
            float _D0 = _b * (sinf(phib) - sinf(phia));
            float _D1 = ((-_b * _b)*(phia - phib) - _b * (_b * cosf(phia))*sinf(phia) + _b * (_b * cosf(phib))*sinf(phib)) / 2;

            /************************************************************************/
            /*                                Part I                                */
            /************************************************************************/
            float CO; // undefined state
            float Df1, Df2;

            /************************************************************************/
            /*                                Part II                               */
            /************************************************************************/
            float _CO;
            float _Df1, _Df2;

            /************************************************************************/
            /*                                Part I                                */
            /************************************************************************/
            float B1 = -cosSigma * D1 + cosW * C1;

            /************************************************************************/
            /*                                Part II                               */
            /************************************************************************/
            float _B1 = _cosW * _C1;

            Lw[0][lobe_i] = sqrtf(1.f / (4.f*M_PI)) * inS0;
            Lw[1][lobe_i] = sqrtf(3.f / (4.f*M_PI)) * inS1;

            float legendreP[(lmax - 1) * 2];
            for (int l = 2; l <= lmax; ++l)
            {
                legendreP[(l - 2) * 2] = LegendreP(l, a + b * cosf(phi1))*b*sinf(phi1) - LegendreP(l, a + b * cosf(phi0))*b*sinf(phi0);
                legendreP[(l - 2) * 2 + 1] = LegendreP(l, _b*cosf(phib))*_b*sinf(phib) - LegendreP(l, _b*cosf(phia))*_b*sinf(phia);
            }

            // TODO: Redundant computation: I mistakenly compute Bl, where we only need Bl_1 actually.
            for (int l = 2; l <= lmax; ++l)
            {
                if (l > 2)
                {
                    float l_1 = l - 1;
                    // Part I
                    C1 = (2.f*l_1 - 1.f) / l_1 * D0 - (l_1 - 1.f) / l_1 * CO;

                    // Part II
                    _C1 = (2.f*l_1 - 1.f) / l_1 * _D0 - (l_1 - 1.f) / l_1 * _CO;
                }
                // Part I
                float E2 = (l % 2 == 0 ? ((2.f*l - 1.f) * C1 + E0) : ((2.f*l - 1.f) * C1 + E1));

                // Part II
                float _E2 = (l % 2 == 0 ? ((2.f*l - 1.f) * _C1 + _E0) : ((2.f*l - 1.f) * _C1 + _E1));

                float F2; // undefined state
                if (l > 3)
                {
                    float l_2 = l - 2;
                    F2 = (l % 2 == 0 ? (2.f*l_2 - 1.f)*Df1 + F0 : (2.f*l_2 - 1.f)*Df2 + F1);
                }

                // Part I
                float Dl = ((a*(2.f*l - 1.f) / l + 2.f*a*(2.f*l - 1.f))*D1 - a * (l - 1.f) / l * C0
                    + legendreP[(l - 2) * 2]
                    + (b*b - a * a - 1)*E2 + l * C1
                    + 2.f*a*(l == 2 ? F0 : (l == 3 ? F1 : F2))) / (l + 1);

                // Part II (Same)
                float _Dl = (legendreP[(l - 2) * 2 + 1] + (_b*_b - 1)*_E2 + l * _C1) / (l + 1);

                // Part I
                float Cl = (2.f*l - 1.f) / l * D1 - (l - 1.f) / l * C0;
                // Part II
                float _Cl = (2.f*l - 1.f) / l * _D1 - (l - 1.f) / l * _C0;

                // Part I
                float Bl = -cosSigma * Dl + cosW * Cl;
                // Part II
                float _Bl = _cosW * _Cl;

                float Sl = (l - 2.f)*(l - 1.f) / (l*(l + 1.f))* (l % 2 == 0 ? inS0 : inS1) + (2.f*l - 1.f) / (l*(l + 1.f))*(B1 + _B1);
                Lw[l][lobe_i] = sqrtf((2.f * l + 1) / (4.f*M_PI)) * Sl;

                B1 = Bl;
                _B1 = _Bl;
                if (l % 2 == 0)
                {
                    inS0 = Sl;
                }
                else
                {
                    inS1 = Sl;
                }

                if (l % 2 == 0)
                {
                    // Part I
                    Df1 = D1;
                    E0 = E2;
                    // Part II
                    _E0 = _E2;
                }
                else
                {
                    // Part I
                    Df2 = D1;
                    E1 = E2;
                    // Part II
                    _E1 = _E2;
                }
                // Part I
                D0 = D1; D1 = Dl;
                CO = C0; C0 = C1;
                // Part II
                _D0 = _D1; _D1 = _Dl;
                _CO = _C0; _C0 = _C1;

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
        }
    }


    /************************************************************************/
    /*                               Client Code                            */
    /************************************************************************/
    /**
      * @brief Encapsulate Sl computation up to $lmax$ . Use this for client.
      *
      * @param[in] sinSigma,cosSigma aperture angle for spherical cap (light source)
      * @param[in] L_world    Light direction in TBN space
      * @param[in] wDir_world Lobe direction wi in TBN space
      * @param[in] lobe_i     lobe direction index
      *
      * @param[out] Lw     2D Array for all SH light coefficients (of all lobe directions).
      */
    template<int lmax>
    bool Sphere_computeLw(
        float sinSigma, float cosSigma,
        const float3 & L_world,
        const float3 & wDir_world,
        int lobe_i,
        float Lw[lmax + 1][2 * lmax + 1])
    {
        float cosW = dot(L_world, wDir_world);
        float sinW = std::sqrt(1.f - cosW * cosW);

        if (M_PI / 2.f - asin(sinSigma) >= std::acos(L_world.z))
        {
            // Full spherical cap case
            tsImpl::Sphere_Lw<lmax>(sinSigma, sinW, cosSigma, cosW, lobe_i, Lw);
            return true;
        }
        else if (asin(sinSigma) + M_PI / 2.f <= std::acos(L_world.z))
        {
            // Spherical cap are completely in lower hemisphere
            return false;
        }
        else
        {
            // Constant precomputation, S0:
            float psaiD = std::acos(L_world.z);
            float S0 = 2.f * M_PI - 2.f*M_PI*cosSigma - 2.f*std::acos(cosf(psaiD) / sinSigma) + 2.f*cosSigma*std::acos(cosSigma*cosf(psaiD) / (sinSigma*sinf(psaiD)));
            float S1;

            float cosWprime = wDir_world.z;
            float sinWprime = std::sqrt(1.f - cosWprime * cosWprime);

            // S1, phi:
            float phi0, phi1, phia, phib;
            tsImpl::Sphere_computeS1_clipped(wDir_world, L_world, sinSigma, cosSigma, sinW, cosW, sinWprime, cosWprime, &S1, &phi0, &phi1, &phia, &phib);

            tsImpl::Sphere_Lw_clipped<lmax>(sinSigma, sinW, cosSigma, cosW, phi0, phi1, phia, phib, cosWprime, sinWprime, S0, S1, lobe_i, Lw);

            return true;
        }
    }

    /**
      * @brief Compute spherical area light coefficients.
      * @note ylmCoeff has to be zero before passing in.
      *
      * @param[in]  L_TBN             Light direction in TBN
      * @param[in]  sinSigma,cosSigma Light spherical cap aperture angle
      * @param[in]  areaLightBasisVector basis directions for ZHF
      * @param[in]  areaLightAlphaCoeff  weighting coefficients for ZHF
      * @param[out] ylmCoeff          Array of length (lmax + 1)*(lmax + 1)
      * @param[out] lightIntensity    light intensity w.r.t this Llm.
      */
    template<int lmax>
    bool Sphere_computeCoeff(const float3 &L_TBN, float sinSigma, float cosSigma,
        float3 *ylmCoeff, const float3 &lightIntensity, float3 *areaLightBasisVector,
        float **areaLightAlphaCoeff)
    {
        float Lw[lmax + 1][2 * lmax + 1];

        // For each lobe direction wi
        for (int lobe_i = 0; lobe_i < 2 * lmax + 1; ++lobe_i)
        {
            const float3 &wi = areaLightBasisVector[lobe_i];
            if (!Sphere_computeLw<lmax>(sinSigma, cosSigma, L_TBN, wi, lobe_i, Lw))
                return false;
        }

        for (int j = 0; j <= lmax; ++j)
        {
            for (int i = 0; i < 2 * j + 1; ++i)
            {
                float coeff = 0.0f;
                for (int k = 0; k < 2 * j + 1; ++k)
                {
                    coeff += areaLightAlphaCoeff[j*j + i][k] * Lw[j][k];
                }
                ylmCoeff[j*j + i] += lightIntensity * coeff;
            }
        }

        return true;
    }
}