#pragma once

// STL
#include <cmath>
#include <cassert>

namespace ts
{
    /************************************************************************/
    /*                              Legendre Polynomial                     */
    /************************************************************************/
#pragma region REGION_LEGENDREP

    // LegendreP
    template<int N>
    inline float LegendreP(float x);

    template<>
    inline float LegendreP<0>(float x)
    {
        return 1.0f;
    }
    template<>
    inline float LegendreP<1>(float x)
    {
        return x;
    }
    template<>
    inline float LegendreP<2>(float x)
    {
        return 0.5f*(3.f*x*x - 1);
    }
    template<>
    inline float LegendreP<3>(float x)
    {
        return 0.5f*(5.f*x*x*x - 3.f*x);
    }
    template<>
    inline float LegendreP<4>(float x)
    {
        return 0.125f*(35.f*x*x*x*x - 30.f*x*x + 3);
    }
    template<>
    inline float LegendreP<5>(float x)
    {
        return 0.125f*(63.f*powf(x, 5) - 70.f*x*x*x + 15.f*x);
    }
    template<>
    inline float LegendreP<6>(float x)
    {
        return (231.f*powf(x, 6) - 315.f*x*x*x*x + 105.f*x*x - 5.f) / 16.f;
    }
    template<>
    inline float LegendreP<7>(float x)
    {
        return (429.f*powf(x, 7) - 693.f*powf(x, 5) + 315.f*powf(x, 3) - 35.f*x) / 16.f;
    }
    template<>
    inline float LegendreP<8>(float x)
    {
        return (6435.f*powf(x, 8) - 12012.f*powf(x, 6) + 6930.f*powf(x, 4) - 1260.f*powf(x, 2) + 35.f) / 128.f;
    }
    template<>
    inline float LegendreP<9>(float x)
    {
        return (12155.f*powf(x, 9) - 25740.f*powf(x, 7) + 18018.f*powf(x, 5) - 4620.f*powf(x, 3) + 315.f*x) / 128.f;
    }
    template<>
    inline float LegendreP<10>(float x)
    {
        return (46189.f*powf(x, 10) - 109395.f*powf(x, 8) + 90090.f*powf(x, 6) - 30030.f*powf(x, 4) + 3465.f*powf(x, 2) - 63.f) / 256.f;
    }

    template<>
    inline float LegendreP<11>(float x)
    {
        return (-693 * x + 15015 * powf(x, 3) - 90090 * powf(x, 5) + 218790 * powf(x, 7) - 230945 * powf(x, 9) + 88179 * powf(x, 11)) / 256;
    }
    template<>
    inline float LegendreP<12>(float x)
    {
        return (231 - 18018 * powf(x, 2) + 225225 * powf(x, 4) - 1021020 * powf(x, 6) + 2078505 * powf(x, 8) - 1939938 * powf(x, 10) + 676039 * powf(x, 12)) / 1024;
    }
    template<>
    inline float LegendreP<13>(float x)
    {
        return (3003 * x - 90090 * powf(x, 3) + 765765 * powf(x, 5) - 2771340 * powf(x, 7) + 4849845 * powf(x, 9) - 4056234 * powf(x, 11) + 1300075 * powf(x, 13)) / 1024;
    }
    template<>
    inline float LegendreP<14>(float x)
    {
        return (-429 + 45045 * powf(x, 2) - 765765 * powf(x, 4) + 4849845 * powf(x, 6) - 14549535 * powf(x, 8) + 22309287 * powf(x, 10) - 16900975 * powf(x, 12) +
            5014575 * powf(x, 14)) / 2048;
    }
    template<>
    inline float LegendreP<15>(float x)
    {
        return (-6435 * x + 255255 * powf(x, 3) - 2909907 * powf(x, 5) + 14549535 * powf(x, 7) - 37182145 * powf(x, 9) + 50702925 * powf(x, 11) - 35102025 * powf(x, 13) +
            9694845 * powf(x, 15)) / 2048;
    }

    /**
     * @brief Legendre polynomial up to l=15.
     */
    inline float LegendreP(int l, float x)
    {
        assert(l >= 0 && l <= 15);
        switch (l)
        {
        case 0:return LegendreP<0>(x);
        case 1:return LegendreP<1>(x);
        case 2:return LegendreP<2>(x);
        case 3:return LegendreP<3>(x);
        case 4:return LegendreP<4>(x);
        case 5:return LegendreP<5>(x);
        case 6:return LegendreP<6>(x);
        case 7:return LegendreP<7>(x);
        case 8:return LegendreP<8>(x);
        case 9:return LegendreP<9>(x);
        case 10:return LegendreP<10>(x);
        case 11:return LegendreP<11>(x);
        case 12:return LegendreP<12>(x);
        case 13:return LegendreP<13>(x);
        case 14:return LegendreP<14>(x);
        case 15:return LegendreP<15>(x);
        default:
            throw std::out_of_range("LegendreP is out of range: 0<=l<=15");
        }
    }

#pragma endregion REGION_LEGENDREP

}