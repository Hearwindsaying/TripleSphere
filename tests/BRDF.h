#pragma once

#include <VecType.h>

namespace ts
{
    /**
     * @brief Sample implementation of BSDFClass class. The 
     * 'Eval_f(wo,wi)' function is requried.
     */
    struct RoughMetal
    {
        RoughMetal(float alpha, float3 eta, float3 kappa) :m_alpha(alpha), m_eta(eta), m_kappa(kappa) {}

        static float GGX_D(const float3 &m, float alphax)
        {
            /*ensure the orientation is consistent with the macrosurface normal,
             *enabling us to leave out the heaviside function in D(m)*/
            if (m.z <= 0.f)
                return 0.f;

            float CosThetaSqr = m.z*m.z;
            float alphaSqr = alphax * alphax;
            float value = 0.f;

            float denom = (1 + BSDFMath::Tan2Theta(m) / alphaSqr);
            value = 1.f / (M_PI * alphaSqr * CosThetaSqr * CosThetaSqr * denom * denom);

            /*potential numeric issues preventing*/
            if (value * BSDFMath::CosTheta(m) < 1e-20f)
                value = 0.f;
            return value;
        }

        static float Smith_G1(const float3 &v, const float3 &m, const float alphax)
        {
            /*ensure consistent orientation(inspired from mitsuba)*/
            if (dot(v, m) * BSDFMath::CosTheta(v) <= 0)
                return 0.0f;

            /* Perpendicular incidence -- no shadowing/masking */
            float tanTheta = fabsf(BSDFMath::TanTheta(v));
            if (tanTheta == 0.0f)
                return 1.0f;

            float root = alphax * tanTheta;
            return 2.0f / (1.0f + sqrtf(1.f + root * root));
        }

        /*Separable Masking and Shadowing Smith_G function*/
        static float Smith_G_Sep(const float3 &wo, const float3 &wi, const float3 &m, const float alpha)
        {
            /*if(print)
                rtPrintf("%f %f woloc:%f %f %f m:%f %f %f wiloc:%f %f %f\n", Smith_G1(wo, m, shaderParams), Smith_G1(wi, m, shaderParams),
                    wo.x,wo.y,wo.z,m.x,m.y,m.z,wi.x,wi.y,wi.z);*/
            return Smith_G1(wo, m, alpha) * Smith_G1(wi, m, alpha);
        }

        static float3 FresnelConductor(float cosi, const float3 & eta, const float3 & k)
        {
            float3 tmp = (eta*eta + k * k)*cosi*cosi;
            float3 Rparl2 = (tmp - (2.f * eta * cosi) + 1.0f) /
                (tmp + (2.f * eta * cosi) + 1.0f);
            float3 tmp_f = eta * eta + k * k;
            float3 Rperp2 =
                (tmp_f - (2.f * eta * cosi) + cosi * cosi) /
                (tmp_f + (2.f * eta * cosi) + cosi * cosi);
            return (Rparl2 + Rperp2) / 2.f;
        }

        static float3 FrCond_Evaluate(float cosi, const float3 & eta, const float3 & k)
        {
            return FresnelConductor(fabsf(cosi), eta, k);
        }

        float3 Eval_f(const float3 & wo_Local, const float3 & wi_Local) const
        {
            /*one-sided BSDF detection*/
            if (BSDFMath::CosTheta(wo_Local) <= 0.f || BSDFMath::CosTheta(wi_Local) <= 0.f)
                return make_float3(0.f, 0.f, 0.f);

            float cosThetaO = BSDFMath::AbsCosTheta(wo_Local);
            float cosThetaI = BSDFMath::AbsCosTheta(wi_Local);

            float3 wh = wi_Local + wo_Local;

            // Handle degenerate cases for microfacet reflection
            if (cosThetaI == 0 || cosThetaO == 0)
            {
                return make_float3(0.f, 0.f, 0.f);
            }
            if (wh.x == 0 && wh.y == 0 && wh.z == 0)
            {
                return make_float3(0.f, 0.f, 0.f);
            }
            wh = safe_normalize(wh);

            //note that when wh is downward,D(wh) could fail

            float3 F = FrCond_Evaluate(dot(wi_Local, wh), m_eta, m_kappa);

            return GGX_D(wh, m_alpha) * Smith_G_Sep(wo_Local, wi_Local, wh, m_alpha) * F / (4 * cosThetaI * cosThetaO);
        }

    private:
        float m_alpha;
        float3 m_eta;
        float3 m_kappa;
    };
}

