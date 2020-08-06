#pragma once

//#define _MATH_DEFINES_DEFINED // for M_PI
#include <corecrt_math_defines.h>
#include <cmath>


namespace ts
{
    /* Vector type used for TripleSphere. */
    struct float3
    {
        float x, y, z;

        float3()
        {
            x = y = z = 0.f;
        }
    };

    float3 make_float3(float x, float y, float z)
    {
        float3 t; t.x = x; t.y = y; t.z = z; return t;
    }

    float3 operator-(const float3& a, const float3& b)
    {
        return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
    }

    float3 operator*(const float3& a, const float3& b)
    {
        return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
    }
    float3 operator*(const float3& a, const float s)
    {
        return make_float3(a.x * s, a.y * s, a.z * s);
    }
    float3 operator*(const float s, const float3& a)
    {
        return make_float3(a.x * s, a.y * s, a.z * s);
    }
    void operator*=(float3& a, const float3& s)
    {
        a.x *= s.x; a.y *= s.y; a.z *= s.z;
    }
    void operator*=(float3& a, const float s)
    {
        a.x *= s; a.y *= s; a.z *= s;
    }


    float3 operator/(const float3& a, const float s)
    {
        float inv = 1.0f / s;
        return a * inv;
    }

    void operator/=(float3& a, const float s)
    {
        float inv = 1.0f / s;
        a *= inv;
    }

    float3 operator+(const float3& a, const float3& b)
    {
        return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
    }
    float3 operator+(const float3& a, const float b)
    {
        return make_float3(a.x + b, a.y + b, a.z + b);
    }
    float3 operator+(const float a, const float3& b)
    {
        return make_float3(a + b.x, a + b.y, a + b.z);
    }
    void operator+=(float3& a, const float3& b)
    {
        a.x += b.x; a.y += b.y; a.z += b.z;
    }

    float dot(const float3& a, const float3& b)
    {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    float3 safe_normalize(const float3& v)
    {
        float d = dot(v, v);
        if (d == 0.0f)
        {
            return v;
        }
        else
        {
            float invLen = 1.0f / sqrtf(d);
            return v * invLen;
        }
    }

    float clamp(const float f, const float a, const float b)
    {
        return fmaxf(a, fminf(f, b));
    }

    float length(const float3& v)
    {
        return sqrtf(dot(v, v));
    }

    float3 cross(const float3& a, const float3& b)
    {
        return make_float3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
    }

    int sign(float a)
    {
        if (a > 0.f)return 1;
        else return -1;
    }

    /* Sampling functions. */
    inline float Lerp(float t, float v1, float v2) { return (1 - t) * v1 + t * v2; }

    float3 UniformSampleSphere(const float u, const float v)
    {
        float z = 1 - 2 * u;
        float r = std::sqrt(fmax(0.f, 1.f - z * z));
        float phi = 2 * M_PI * v;
        return make_float3(r * std::cos(phi), r * std::sin(phi), z);
    }
    inline float pdfSphere() { return 1.f / (4.f * M_PI); };

    inline float3 uniformSamplingHemisphere(float x, float y)
    {
        float z = x;
        float r = std::sqrt(fmax(0.0, 1.0 - z * z));
        float phi = 2.0 * M_PI * y;
        return make_float3(r * std::cos(phi), r * std::sin(phi), z);
    };
    inline float pdfHemiSphere() { return 1.f / (2.f * M_PI); };


    float3 UniformSampleCone(const float u, const float v, float cosThetaMax)
    {
        float cosTheta = (1.f - u) + u * cosThetaMax;
        float sinTheta = std::sqrt(1.f - cosTheta * cosTheta);
        float phi = v * 2.f * M_PI;
        return make_float3(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta,
            cosTheta);
    }

    float3 UniformSampleCone(const float u, const float v, float cosThetaMax,
        const float3 &x, const float3 &y, const float3 &z) 
    {
        float cosTheta = Lerp(u, cosThetaMax, 1.f);
        float sinTheta = std::sqrt(1.f - cosTheta * cosTheta);
        float phi = v * 2 * M_PI;
        return std::cos(phi) * sinTheta * x + std::sin(phi) * sinTheta * y + cosTheta * z;
    }

    float UniformConePdf(float cosThetaMax) 
    {
        return 1.f / (2.f * M_PI * (1.f - cosThetaMax));
    }

    /* Local coordinates (TBN). */
    namespace BSDFMath
    {
        //////////////////////////////////////////////////////////////////////////
        //Forward declarations:
        //BSDFMath:
        float AbsCosTheta(const float3 & w_local);
        float CosTheta(const float3 & w_local);
        float Cos2Theta(const float3 & w_local);
        float Sin2Theta(const float3 & w_local);
        float SinTheta(const float3 &w_local);
        float Tan2Theta(const float3 &w_local);
        float CosPhi(const float3 &w_local);
        float SinPhi(const float3 &w_local);
        float TanTheta(const float3 &w_local);
        float Cos2Phi(const float3 &w_local);
        float Sin2Phi(const float3 &w_local);

        //BSDF Space conversion:
        float3 WorldToLocal(const float3 & v, const float3 & sn, const float3 & tn, const float3 & nn);
        float3 LocalToWorld(const float3 & v, const float3 & sn, const float3 & tn, const float3 & nn);
        inline void CoordinateSystem(const float3 &v1, float3 *v2, float3 *v3);

        //////////////////////////////////////////////////////////////////////////
        //BSDFMath:
        float AbsCosTheta(const float3 & w_local)
        {
            return fabsf(w_local.z);
        }

        float CosTheta(const float3 & w_local)
        {
            return w_local.z;
        }

        float Cos2Theta(const float3 & w_local)
        {
            return w_local.z * w_local.z;
        }

        float Sin2Theta(const float3 & w_local)
        {
            return fmaxf(0.f, 1 - Cos2Theta(w_local));
        }

        float SinTheta(const float3 &w_local)
        {
            return sqrtf(Sin2Theta(w_local));
        }

        float Tan2Theta(const float3 &w_local)
        {
            return Sin2Theta(w_local) / Cos2Theta(w_local);
        }

        float CosPhi(const float3 &w_local)
        {
            float sinTheta = SinTheta(w_local);
            return (sinTheta == 0) ? 1 : clamp(w_local.x / sinTheta, -1.f, 1.f);
        }

        float SinPhi(const float3 &w_local)
        {
            float sinTheta = SinTheta(w_local);
            return (sinTheta == 0) ? 0 : clamp(w_local.y / sinTheta, -1.f, 1.f);
        }

        float TanTheta(const float3 &w_local)
        {
            return SinTheta(w_local) / CosTheta(w_local);
        }

        float Cos2Phi(const float3 &w_local)
        {
            return CosPhi(w_local) * CosPhi(w_local);
        }

        float Sin2Phi(const float3 &w_local)
        {
            return SinPhi(w_local) * SinPhi(w_local);
        }


        //BSDF space conversion:
        float3 WorldToLocal(const float3 & v, const float3 & sn, const float3 & tn, const float3 & nn)
        {
            return make_float3(dot(v, sn), dot(v, tn), dot(v, nn));
        }
        float3 LocalToWorld(const float3 & v, const float3 & sn, const float3 & tn, const float3 & nn)
        {
            return make_float3(sn.x * v.x + tn.x * v.y + nn.x * v.z,
                sn.y * v.x + tn.y * v.y + nn.y * v.z,
                sn.z * v.x + tn.z * v.y + nn.z * v.z);
        }

        inline void CoordinateSystem(const float3 &v1, float3 *v2, float3 *v3)
        {
            if (std::abs(v1.x) > std::abs(v1.y))
                *v2 = make_float3(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
            else
                *v2 = make_float3(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
            *v3 = cross(v1, *v2);
        }
    }
}