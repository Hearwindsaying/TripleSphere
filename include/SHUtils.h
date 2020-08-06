#pragma once

// STL
#include <cmath>
#include <cassert>
#include <thread>
#include <vector>
#include <random>
#include <iostream>
#include <fstream>

#include <Serialize.h>

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

    /************************************************************************/
    /*                              Spherical Harmonics                     */
    /************************************************************************/
    void SHEvalFast9(const float3 &w, float *pOut) 
    {
        const float fX = w.x;
        const float fY = w.y;
        const float fZ = w.z;

        float fC0, fS0, fC1, fS1, fPa, fPb, fPc;
        {
            float fZ2 = fZ * fZ;
            pOut[0] = 0.282094791774;
            pOut[2] = 0.488602511903f*fZ;
            pOut[6] = (0.946174695758f*fZ2) + -0.315391565253f;
            pOut[12] = fZ * ((fZ2* 1.86588166295f) + -1.11952899777f);
            pOut[20] = ((fZ* 1.9843134833f)*pOut[12]) + (-1.00623058987f*pOut[6]);
            pOut[30] = ((fZ* 1.98997487421f)*pOut[20]) + (-1.00285307284f*pOut[12]);
            pOut[42] = ((fZ* 1.99304345718f)*pOut[30]) + (-1.00154202096f*pOut[20]);
            pOut[56] = ((fZ* 1.99489143482f)*pOut[42]) + (-1.00092721392f*pOut[30]);
            pOut[72] = ((fZ* 1.99608992783f)*pOut[56]) + (-1.00060078107f*pOut[42]);
            pOut[90] = ((fZ* 1.99691119507f)*pOut[72]) + (-1.00041143799f*pOut[56]);
            fC0 = fX;
            fS0 = fY;
            fPa = -0.488602511903f;
            pOut[3] = fPa * fC0;
            pOut[1] = fPa * fS0;
            fPb = -1.09254843059f*fZ;
            pOut[7] = fPb * fC0;
            pOut[5] = fPb * fS0;
            fPc = (-2.28522899732f*fZ2) + 0.457045799464f;
            pOut[13] = fPc * fC0;
            pOut[11] = fPc * fS0;
            fPa = fZ * ((fZ2* -4.6833258049f) + 2.00713963067f);
            pOut[21] = fPa * fC0;
            pOut[19] = fPa * fS0;
            fPb = ((fZ* 2.03100960116f)*fPa) + (-0.991031208965f*fPc);
            pOut[31] = fPb * fC0;
            pOut[29] = fPb * fS0;
            fPc = ((fZ* 2.02131498924f)*fPb) + (-0.995226703056f*fPa);
            pOut[43] = fPc * fC0;
            pOut[41] = fPc * fS0;
            fPa = ((fZ* 2.01556443707f)*fPc) + (-0.997155044022f*fPb);
            pOut[57] = fPa * fC0;
            pOut[55] = fPa * fS0;
            fPb = ((fZ* 2.01186954041f)*fPa) + (-0.99816681789f*fPc);
            pOut[73] = fPb * fC0;
            pOut[71] = fPb * fS0;
            fPc = ((fZ* 2.00935312974f)*fPb) + (-0.998749217772f*fPa);
            pOut[91] = fPc * fC0;
            pOut[89] = fPc * fS0;
            fC1 = (fX*fC0) - (fY*fS0);
            fS1 = (fX*fS0) + (fY*fC0);
            fPa = 0.546274215296f;
            pOut[8] = fPa * fC1;
            pOut[4] = fPa * fS1;
            fPb = 1.44530572132f*fZ;
            pOut[14] = fPb * fC1;
            pOut[10] = fPb * fS1;
            fPc = (3.31161143515f*fZ2) + -0.473087347879f;
            pOut[22] = fPc * fC1;
            pOut[18] = fPc * fS1;
            fPa = fZ * ((fZ2* 7.19030517746f) + -2.39676839249f);
            pOut[32] = fPa * fC1;
            pOut[28] = fPa * fS1;
            fPb = ((fZ* 2.11394181566f)*fPa) + (-0.973610120462f*fPc);
            pOut[44] = fPb * fC1;
            pOut[40] = fPb * fS1;
            fPc = ((fZ* 2.08166599947f)*fPb) + (-0.984731927835f*fPa);
            pOut[58] = fPc * fC1;
            pOut[54] = fPc * fS1;
            fPa = ((fZ* 2.06155281281f)*fPc) + (-0.99033793766f*fPb);
            pOut[74] = fPa * fC1;
            pOut[70] = fPa * fS1;
            fPb = ((fZ* 2.04812235836f)*fPa) + (-0.99348527267f*fPc);
            pOut[92] = fPb * fC1;
            pOut[88] = fPb * fS1;
            fC0 = (fX*fC1) - (fY*fS1);
            fS0 = (fX*fS1) + (fY*fC1);
            fPa = -0.590043589927f;
            pOut[15] = fPa * fC0;
            pOut[9] = fPa * fS0;
            fPb = -1.77013076978f*fZ;
            pOut[23] = fPb * fC0;
            pOut[17] = fPb * fS0;
            fPc = (-4.40314469492f*fZ2) + 0.489238299435f;
            pOut[33] = fPc * fC0;
            pOut[27] = fPc * fS0;
            fPa = fZ * ((fZ2* -10.1332578547f) + 2.76361577854f);
            pOut[45] = fPa * fC0;
            pOut[39] = fPa * fS0;
            fPb = ((fZ* 2.20794021658f)*fPa) + (-0.9594032236f*fPc);
            pOut[59] = fPb * fC0;
            pOut[53] = fPb * fS0;
            fPc = ((fZ* 2.1532216877f)*fPb) + (-0.97521738656f*fPa);
            pOut[75] = fPc * fC0;
            pOut[69] = fPc * fS0;
            fPa = ((fZ* 2.11804417119f)*fPc) + (-0.983662844979f*fPb);
            pOut[93] = fPa * fC0;
            pOut[87] = fPa * fS0;
            fC1 = (fX*fC0) - (fY*fS0);
            fS1 = (fX*fS0) + (fY*fC0);
            fPa = 0.625835735449f;
            pOut[24] = fPa * fC1;
            pOut[16] = fPa * fS1;
            fPb = 2.07566231488f*fZ;
            pOut[34] = fPb * fC1;
            pOut[26] = fPb * fS1;
            fPc = (5.55021390802f*fZ2) + -0.504564900729f;
            pOut[46] = fPc * fC1;
            pOut[38] = fPc * fS1;
            fPa = fZ * ((fZ2* 13.4918050467f) + -3.11349347232f);
            pOut[60] = fPa * fC1;
            pOut[52] = fPa * fS1;
            fPb = ((fZ* 2.30488611432f)*fPa) + (-0.948176387355f*fPc);
            pOut[76] = fPb * fC1;
            pOut[68] = fPb * fS1;
            fPc = ((fZ* 2.22917715071f)*fPb) + (-0.967152839723f*fPa);
            pOut[94] = fPc * fC1;
            pOut[86] = fPc * fS1;
            fC0 = (fX*fC1) - (fY*fS1);
            fS0 = (fX*fS1) + (fY*fC1);
            fPa = -0.65638205684f;
            pOut[35] = fPa * fC0;
            pOut[25] = fPa * fS0;
            fPb = -2.36661916223f*fZ;
            pOut[47] = fPb * fC0;
            pOut[37] = fPb * fS0;
            fPc = (-6.74590252336f*fZ2) + 0.51891557872f;
            pOut[61] = fPc * fC0;
            pOut[51] = fPc * fS0;
            fPa = fZ * ((fZ2* -17.2495531105f) + 3.4499106221f);
            pOut[77] = fPa * fC0;
            pOut[67] = fPa * fS0;
            fPb = ((fZ* 2.40163634692f)*fPa) + (-0.939224604204f*fPc);
            pOut[95] = fPb * fC0;
            pOut[85] = fPb * fS0;
            fC1 = (fX*fC0) - (fY*fS0);
            fS1 = (fX*fS0) + (fY*fC0);
            fPa = 0.683184105192f;
            pOut[48] = fPa * fC1;
            pOut[36] = fPa * fS1;
            fPb = 2.6459606618f*fZ;
            pOut[62] = fPb * fC1;
            pOut[50] = fPb * fS1;
            fPc = (7.98499149089f*fZ2) + -0.53233276606f;
            pOut[78] = fPc * fC1;
            pOut[66] = fPc * fS1;
            fPa = fZ * ((fZ2* 21.3928901909f) + -3.77521591604f);
            pOut[96] = fPa * fC1;
            pOut[84] = fPa * fS1;
            fC0 = (fX*fC1) - (fY*fS1);
            fS0 = (fX*fS1) + (fY*fC1);
            fPa = -0.707162732525f;
            pOut[63] = fPa * fC0;
            pOut[49] = fPa * fS0;
            fPb = -2.9157066407f*fZ;
            pOut[79] = fPb * fC0;
            pOut[65] = fPb * fS0;
            fPc = (-9.26339318285f*fZ2) + 0.544905481344f;
            pOut[97] = fPc * fC0;
            pOut[83] = fPc * fS0;
            fC1 = (fX*fC0) - (fY*fS0);
            fS1 = (fX*fS0) + (fY*fC0);
            fPa = 0.728926660175f;
            pOut[80] = fPa * fC1;
            pOut[64] = fPa * fS1;
            fPb = 3.17731764895f*fZ;
            pOut[98] = fPb * fC1;
            pOut[82] = fPb * fS1;
            fC0 = (fX*fC1) - (fY*fS1);
            fS0 = (fX*fS1) + (fY*fC1);
            fPc = -0.748900951853f;
            pOut[99] = fPc * fC0;
            pOut[81] = fPc * fS0;
        }
    }

    void SHEvalFast15(const float3 &w, float *pOut)
    {
        const float fX = w.x;
        const float fY = w.y;
        const float fZ = w.z;
        float fC0, fS0, fC1, fS1, fPa, fPb, fPc;
        {
            float fZ2 = fZ * fZ;
            pOut[0] = 0.282094791774;
            pOut[2] = 0.488602511903f*fZ;
            pOut[6] = (0.946174695758f*fZ2) + -0.315391565253f;
            pOut[12] = fZ * ((fZ2* 1.86588166295f) + -1.11952899777f);
            pOut[20] = ((fZ* 1.9843134833f)*pOut[12]) + (-1.00623058987f*pOut[6]);
            pOut[30] = ((fZ* 1.98997487421f)*pOut[20]) + (-1.00285307284f*pOut[12]);
            pOut[42] = ((fZ* 1.99304345718f)*pOut[30]) + (-1.00154202096f*pOut[20]);
            pOut[56] = ((fZ* 1.99489143482f)*pOut[42]) + (-1.00092721392f*pOut[30]);
            pOut[72] = ((fZ* 1.99608992783f)*pOut[56]) + (-1.00060078107f*pOut[42]);
            pOut[90] = ((fZ* 1.99691119507f)*pOut[72]) + (-1.00041143799f*pOut[56]);
            pOut[110] = ((fZ* 1.99749843554f)*pOut[90]) + (-1.00029407441f*pOut[72]);
            pOut[132] = ((fZ* 1.99793281599f)*pOut[110]) + (-1.00021746222f*pOut[90]);
            pOut[156] = ((fZ* 1.99826313471f)*pOut[132]) + (-1.00016533025f*pOut[110]);
            pOut[182] = ((fZ* 1.99852016258f)*pOut[156]) + (-1.00012862564f*pOut[132]);
            pOut[210] = ((fZ* 1.9987240828f)*pOut[182]) + (-1.00010203561f*pOut[156]);
            pOut[240] = ((fZ* 1.99888858008f)*pOut[210]) + (-1.00008230114f*pOut[182]);
            fC0 = fX;
            fS0 = fY;
            fPa = -0.488602511903f;
            pOut[3] = fPa * fC0;
            pOut[1] = fPa * fS0;
            fPb = -1.09254843059f*fZ;
            pOut[7] = fPb * fC0;
            pOut[5] = fPb * fS0;
            fPc = (-2.28522899732f*fZ2) + 0.457045799464f;
            pOut[13] = fPc * fC0;
            pOut[11] = fPc * fS0;
            fPa = fZ * ((fZ2* -4.6833258049f) + 2.00713963067f);
            pOut[21] = fPa * fC0;
            pOut[19] = fPa * fS0;
            fPb = ((fZ* 2.03100960116f)*fPa) + (-0.991031208965f*fPc);
            pOut[31] = fPb * fC0;
            pOut[29] = fPb * fS0;
            fPc = ((fZ* 2.02131498924f)*fPb) + (-0.995226703056f*fPa);
            pOut[43] = fPc * fC0;
            pOut[41] = fPc * fS0;
            fPa = ((fZ* 2.01556443707f)*fPc) + (-0.997155044022f*fPb);
            pOut[57] = fPa * fC0;
            pOut[55] = fPa * fS0;
            fPb = ((fZ* 2.01186954041f)*fPa) + (-0.99816681789f*fPc);
            pOut[73] = fPb * fC0;
            pOut[71] = fPb * fS0;
            fPc = ((fZ* 2.00935312974f)*fPb) + (-0.998749217772f*fPa);
            pOut[91] = fPc * fC0;
            pOut[89] = fPc * fS0;
            fPa = ((fZ* 2.00756146364f)*fPc) + (-0.999108336871f*fPb);
            pOut[111] = fPa * fC0;
            pOut[109] = fPa * fS0;
            fPb = ((fZ* 2.00624026477f)*fPa) + (-0.999341888708f*fPc);
            pOut[133] = fPb * fC0;
            pOut[131] = fPb * fS0;
            fPc = ((fZ* 2.00523789636f)*fPb) + (-0.999500374688f*fPa);
            pOut[157] = fPc * fC0;
            pOut[155] = fPc * fS0;
            fPa = ((fZ* 2.00445931434f)*fPc) + (-0.999611725864f*fPb);
            pOut[183] = fPa * fC0;
            pOut[181] = fPa * fS0;
            fPb = ((fZ* 2.00384246272f)*fPa) + (-0.99969226034f*fPc);
            pOut[211] = fPb * fC0;
            pOut[209] = fPb * fS0;
            fPc = ((fZ* 2.00334541633f)*fPb) + (-0.999751953363f*fPa);
            pOut[241] = fPc * fC0;
            pOut[239] = fPc * fS0;
            fC1 = (fX*fC0) - (fY*fS0);
            fS1 = (fX*fS0) + (fY*fC0);
            fPa = 0.546274215296f;
            pOut[8] = fPa * fC1;
            pOut[4] = fPa * fS1;
            fPb = 1.44530572132f*fZ;
            pOut[14] = fPb * fC1;
            pOut[10] = fPb * fS1;
            fPc = (3.31161143515f*fZ2) + -0.473087347879f;
            pOut[22] = fPc * fC1;
            pOut[18] = fPc * fS1;
            fPa = fZ * ((fZ2* 7.19030517746f) + -2.39676839249f);
            pOut[32] = fPa * fC1;
            pOut[28] = fPa * fS1;
            fPb = ((fZ* 2.11394181566f)*fPa) + (-0.973610120462f*fPc);
            pOut[44] = fPb * fC1;
            pOut[40] = fPb * fS1;
            fPc = ((fZ* 2.08166599947f)*fPb) + (-0.984731927835f*fPa);
            pOut[58] = fPc * fC1;
            pOut[54] = fPc * fS1;
            fPa = ((fZ* 2.06155281281f)*fPc) + (-0.99033793766f*fPb);
            pOut[74] = fPa * fC1;
            pOut[70] = fPa * fS1;
            fPb = ((fZ* 2.04812235836f)*fPa) + (-0.99348527267f*fPc);
            pOut[92] = fPb * fC1;
            pOut[88] = fPb * fS1;
            fPc = ((fZ* 2.03868830379f)*fPb) + (-0.99539380324f*fPa);
            pOut[112] = fPc * fC1;
            pOut[108] = fPc * fS1;
            fPa = ((fZ* 2.03179849596f)*fPc) + (-0.996620470226f*fPb);
            pOut[134] = fPa * fC1;
            pOut[130] = fPa * fS1;
            fPb = ((fZ* 2.02660870844f)*fPa) + (-0.997445717412f*fPc);
            pOut[158] = fPb * fC1;
            pOut[154] = fPb * fS1;
            fPc = ((fZ* 2.02259958739f)*fPb) + (-0.998021758696f*fPa);
            pOut[184] = fPc * fC1;
            pOut[180] = fPc * fS1;
            fPa = ((fZ* 2.01943680268f)*fPc) + (-0.998436277386f*fPb);
            pOut[212] = fPa * fC1;
            pOut[208] = fPa * fS1;
            fPb = ((fZ* 2.01689694907f)*fPa) + (-0.998742296069f*fPc);
            pOut[242] = fPb * fC1;
            pOut[238] = fPb * fS1;
            fC0 = (fX*fC1) - (fY*fS1);
            fS0 = (fX*fS1) + (fY*fC1);
            fPa = -0.590043589927f;
            pOut[15] = fPa * fC0;
            pOut[9] = fPa * fS0;
            fPb = -1.77013076978f*fZ;
            pOut[23] = fPb * fC0;
            pOut[17] = fPb * fS0;
            fPc = (-4.40314469492f*fZ2) + 0.489238299435f;
            pOut[33] = fPc * fC0;
            pOut[27] = fPc * fS0;
            fPa = fZ * ((fZ2* -10.1332578547f) + 2.76361577854f);
            pOut[45] = fPa * fC0;
            pOut[39] = fPa * fS0;
            fPb = ((fZ* 2.20794021658f)*fPa) + (-0.9594032236f*fPc);
            pOut[59] = fPb * fC0;
            pOut[53] = fPb * fS0;
            fPc = ((fZ* 2.1532216877f)*fPb) + (-0.97521738656f*fPa);
            pOut[75] = fPc * fC0;
            pOut[69] = fPc * fS0;
            fPa = ((fZ* 2.11804417119f)*fPc) + (-0.983662844979f*fPb);
            pOut[93] = fPa * fC0;
            pOut[87] = fPa * fS0;
            fPb = ((fZ* 2.09394732136f)*fPa) + (-0.988623065486f*fPc);
            pOut[113] = fPb * fC0;
            pOut[107] = fPb * fS0;
            fPc = ((fZ* 2.07665596573f)*fPb) + (-0.991742220327f*fPa);
            pOut[135] = fPc * fC0;
            pOut[129] = fPc * fS0;
            fPa = ((fZ* 2.06379729122f)*fPc) + (-0.99380799f*fPb);
            pOut[159] = fPa * fC0;
            pOut[153] = fPa * fS0;
            fPb = ((fZ* 2.05395959064f)*fPa) + (-0.995233204046f*fPc);
            pOut[185] = fPb * fC0;
            pOut[179] = fPb * fS0;
            fPc = ((fZ* 2.04625652727f)*fPb) + (-0.996249651937f*fPa);
            pOut[213] = fPc * fC0;
            pOut[207] = fPc * fS0;
            fPa = ((fZ* 2.04010711411f)*fPc) + (-0.996994798511f*fPb);
            pOut[243] = fPa * fC0;
            pOut[237] = fPa * fS0;
            fC1 = (fX*fC0) - (fY*fS0);
            fS1 = (fX*fS0) + (fY*fC0);
            fPa = 0.625835735449f;
            pOut[24] = fPa * fC1;
            pOut[16] = fPa * fS1;
            fPb = 2.07566231488f*fZ;
            pOut[34] = fPb * fC1;
            pOut[26] = fPb * fS1;
            fPc = (5.55021390802f*fZ2) + -0.504564900729f;
            pOut[46] = fPc * fC1;
            pOut[38] = fPc * fS1;
            fPa = fZ * ((fZ2* 13.4918050467f) + -3.11349347232f);
            pOut[60] = fPa * fC1;
            pOut[52] = fPa * fS1;
            fPb = ((fZ* 2.30488611432f)*fPa) + (-0.948176387355f*fPc);
            pOut[76] = fPb * fC1;
            pOut[68] = fPb * fS1;
            fPc = ((fZ* 2.22917715071f)*fPb) + (-0.967152839723f*fPa);
            pOut[94] = fPc * fC1;
            pOut[86] = fPc * fS1;
            fPa = ((fZ* 2.17944947177f)*fPc) + (-0.977692361094f*fPb);
            pOut[114] = fPa * fC1;
            pOut[106] = fPa * fS1;
            fPb = ((fZ* 2.14476105895f)*fPa) + (-0.984083864633f*fPc);
            pOut[136] = fPb * fC1;
            pOut[128] = fPb * fS1;
            fPc = ((fZ* 2.11947811973f)*fPb) + (-0.988211768803f*fPa);
            pOut[160] = fPc * fC1;
            pOut[152] = fPc * fS1;
            fPa = ((fZ* 2.10042012604f)*fPc) + (-0.991008166818f*fPb);
            pOut[186] = fPa * fC1;
            pOut[178] = fPa * fS1;
            fPb = ((fZ* 2.08566536146f)*fPa) + (-0.992975326985f*fPc);
            pOut[214] = fPb * fC1;
            pOut[206] = fPb * fS1;
            fPc = ((fZ* 2.07399021374f)*fPb) + (-0.994402195129f*fPa);
            pOut[244] = fPc * fC1;
            pOut[236] = fPc * fS1;
            fC0 = (fX*fC1) - (fY*fS1);
            fS0 = (fX*fS1) + (fY*fC1);
            fPa = -0.65638205684f;
            pOut[35] = fPa * fC0;
            pOut[25] = fPa * fS0;
            fPb = -2.36661916223f*fZ;
            pOut[47] = fPb * fC0;
            pOut[37] = fPb * fS0;
            fPc = (-6.74590252336f*fZ2) + 0.51891557872f;
            pOut[61] = fPc * fC0;
            pOut[51] = fPc * fS0;
            fPa = fZ * ((fZ2* -17.2495531105f) + 3.4499106221f);
            pOut[77] = fPa * fC0;
            pOut[67] = fPa * fS0;
            fPb = ((fZ* 2.40163634692f)*fPa) + (-0.939224604204f*fPc);
            pOut[95] = fPb * fC0;
            pOut[85] = fPb * fS0;
            fPc = ((fZ* 2.30651251893f)*fPb) + (-0.960392076798f*fPa);
            pOut[115] = fPc * fC0;
            pOut[105] = fPc * fS0;
            fPa = ((fZ* 2.24304480562f)*fPc) + (-0.972483256519f*fPb);
            pOut[137] = fPa * fC0;
            pOut[127] = fPa * fS0;
            fPb = ((fZ* 2.19816577471f)*fPa) + (-0.9799919151f*fPc);
            pOut[161] = fPb * fC0;
            pOut[151] = fPb * fS0;
            fPc = ((fZ* 2.16506350946f)*fPb) + (-0.984940960491f*fPa);
            pOut[187] = fPc * fC0;
            pOut[177] = fPc * fS0;
            fPa = ((fZ* 2.13984751055f)*fPc) + (-0.988353228994f*fPb);
            pOut[215] = fPa * fC0;
            pOut[205] = fPa * fS0;
            fPb = ((fZ* 2.12014150471f)*fPa) + (-0.990790929847f*fPc);
            pOut[245] = fPb * fC0;
            pOut[235] = fPb * fS0;
            fC1 = (fX*fC0) - (fY*fS0);
            fS1 = (fX*fS0) + (fY*fC0);
            fPa = 0.683184105192f;
            pOut[48] = fPa * fC1;
            pOut[36] = fPa * fS1;
            fPb = 2.6459606618f*fZ;
            pOut[62] = fPb * fC1;
            pOut[50] = fPb * fS1;
            fPc = (7.98499149089f*fZ2) + -0.53233276606f;
            pOut[78] = fPc * fC1;
            pOut[66] = fPc * fS1;
            fPa = fZ * ((fZ2* 21.3928901909f) + -3.77521591604f);
            pOut[96] = fPa * fC1;
            pOut[84] = fPa * fS1;
            fPb = ((fZ* 2.49687304443f)*fPa) + (-0.931968978277f*fPc);
            pOut[116] = fPb * fC1;
            pOut[104] = fPb * fS1;
            fPc = ((fZ* 2.38376864254f)*fPb) + (-0.954701580788f*fPa);
            pOut[138] = fPc * fC1;
            pOut[126] = fPc * fS1;
            fPa = ((fZ* 2.30739551748f)*fPc) + (-0.967961183941f*fPb);
            pOut[162] = fPa * fC1;
            pOut[150] = fPa * fS1;
            fPb = ((fZ* 2.25281778445f)*fPa) + (-0.976346606979f*fPc);
            pOut[188] = fPb * fC1;
            pOut[176] = fPb * fS1;
            fPc = ((fZ* 2.21218218056f)*fPb) + (-0.981962321069f*fPa);
            pOut[216] = fPc * fC1;
            pOut[204] = fPc * fS1;
            fPa = ((fZ* 2.1809662438f)*fPc) + (-0.985889075035f*fPb);
            pOut[246] = fPa * fC1;
            pOut[234] = fPa * fS1;
            fC0 = (fX*fC1) - (fY*fS1);
            fS0 = (fX*fS1) + (fY*fC1);
            fPa = -0.707162732525f;
            pOut[63] = fPa * fC0;
            pOut[49] = fPa * fS0;
            fPb = -2.9157066407f*fZ;
            pOut[79] = fPb * fC0;
            pOut[65] = fPb * fS0;
            fPc = (-9.26339318285f*fZ2) + 0.544905481344f;
            pOut[97] = fPc * fC0;
            pOut[83] = fPc * fS0;
            fPa = fZ * ((fZ2* -25.9102413134f) + 4.09109073369f);
            pOut[117] = fPa * fC0;
            pOut[103] = fPa * fS0;
            fPb = ((fZ* 2.59004504465f)*fPa) + (-0.925989276585f*fPc);
            pOut[139] = fPb * fC0;
            pOut[125] = fPb * fS0;
            fPc = ((fZ* 2.46020966158f)*fPb) + (-0.949871380292f*fPa);
            pOut[163] = fPc * fC0;
            pOut[149] = fPc * fS0;
            fPa = ((fZ* 2.37170824513f)*fPc) + (-0.964026880376f*fPb);
            pOut[189] = fPa * fC0;
            pOut[175] = fPa * fS0;
            fPb = ((fZ* 2.30792777449f)*fPa) + (-0.973107792339f*fPc);
            pOut[217] = fPb * fC0;
            pOut[203] = fPb * fS0;
            fPc = ((fZ* 2.2600784379f)*fPb) + (-0.979267402942f*fPa);
            pOut[247] = fPc * fC0;
            pOut[233] = fPc * fS0;
            fC1 = (fX*fC0) - (fY*fS0);
            fS1 = (fX*fS0) + (fY*fC0);
            fPa = 0.728926660175f;
            pOut[80] = fPa * fC1;
            pOut[64] = fPa * fS1;
            fPb = 3.17731764895f*fZ;
            pOut[98] = fPb * fC1;
            pOut[82] = fPb * fS1;
            fPc = (10.5778117217f*fZ2) + -0.55672693272f;
            pOut[118] = fPc * fC1;
            pOut[102] = fPc * fS1;
            fPa = fZ * ((fZ2* 30.7915797034f) + -4.39879710048f);
            pOut[140] = fPa * fC1;
            pOut[124] = fPa * fS1;
            fPb = ((fZ* 2.68095132369f)*fPa) + (-0.920985497016f*fPc);
            pOut[164] = fPb * fC1;
            pOut[148] = fPb * fS1;
            fPc = ((fZ* 2.53546276419f)*fPb) + (-0.945732487487f*fPa);
            pOut[190] = fPc * fC1;
            pOut[174] = fPc * fS1;
            fPa = ((fZ* 2.43553242266f)*fPc) + (-0.960586941785f*fPb);
            pOut[218] = fPa * fC1;
            pOut[202] = fPa * fS1;
            fPb = ((fZ* 2.3630173363f)*fPa) + (-0.970226187228f*fPc);
            pOut[248] = fPb * fC1;
            pOut[232] = fPb * fS1;
            fC0 = (fX*fC1) - (fY*fS1);
            fS0 = (fX*fS1) + (fY*fC1);
            fPa = -0.748900951853f;
            pOut[99] = fPa * fC0;
            pOut[81] = fPa * fS0;
            fPb = -3.43189529989f*fZ;
            pOut[119] = fPb * fC0;
            pOut[101] = fPb * fS0;
            fPc = (-11.9255275395f*fZ2) + 0.567882263783f;
            pOut[141] = fPc * fC0;
            pOut[123] = fPc * fS0;
            fPa = fZ * ((fZ2* -36.0280906893f) + 4.69931617687f);
            pOut[165] = fPa * fC0;
            pOut[147] = fPa * fS0;
            fPb = ((fZ* 2.76955854703f)*fPa) + (-0.916741522875f*fPc);
            pOut[191] = fPb * fC0;
            pOut[173] = fPb * fS0;
            fPc = ((fZ* 2.60934774459f)*fPb) + (-0.942152946136f*fPa);
            pOut[219] = fPc * fC0;
            pOut[201] = fPc * fS0;
            fPa = ((fZ* 2.49861072509f)*fPc) + (-0.957561417515f*fPb);
            pOut[249] = fPa * fC0;
            pOut[231] = fPa * fS0;
            fC1 = (fX*fC0) - (fY*fS0);
            fS1 = (fX*fS0) + (fY*fC0);
            fPa = 0.767395118222f;
            pOut[120] = fPa * fC1;
            pOut[100] = fPa * fS1;
            fPb = 3.68029769881f*fZ;
            pOut[142] = fPb * fC1;
            pOut[122] = fPb * fS1;
            fPc = (13.3042542003f*fZ2) + -0.578445834794f;
            pOut[166] = fPc * fC1;
            pOut[146] = fPc * fS1;
            fPa = fZ * ((fZ2* 41.6119315355f) + -4.99343178426f);
            pOut[192] = fPa * fC1;
            pOut[172] = fPa * fS1;
            fPb = ((fZ* 2.8559149147f)*fPa) + (-0.913099118387f*fPc);
            pOut[220] = fPb * fC1;
            pOut[200] = fPb * fS1;
            fPc = ((fZ* 2.6817904467f)*fPb) + (-0.939030232622f*fPa);
            pOut[250] = fPc * fC1;
            pOut[230] = fPc * fS1;
            fC0 = (fX*fC1) - (fY*fS1);
            fS0 = (fX*fS1) + (fY*fC1);
            fPa = -0.784642105787f;
            pOut[143] = fPa * fC0;
            pOut[121] = fPa * fS0;
            fPb = -3.92321052894f*fZ;
            pOut[167] = fPb * fC0;
            pOut[145] = fPb * fS0;
            fPc = (-14.7120394835f*fZ2) + 0.58848157934f;
            pOut[193] = fPc * fC0;
            pOut[171] = fPc * fS0;
            fPa = fZ * ((fZ2* -47.5360543607f) + 5.28178381785f);
            pOut[221] = fPa * fC0;
            pOut[199] = fPa * fS0;
            fPb = ((fZ* 2.94010727172f)*fPa) + (-0.909940356832f*fPc);
            pOut[251] = fPb * fC0;
            pOut[229] = fPb * fS0;
            fC1 = (fX*fC0) - (fY*fS0);
            fS1 = (fX*fS0) + (fY*fC0);
            fPa = 0.800821995784f;
            pOut[168] = fPa * fC1;
            pOut[144] = fPa * fS1;
            fPb = 4.16119315355f*fZ;
            pOut[194] = fPb * fC1;
            pOut[170] = fPb * fS1;
            fPc = (16.1471947939f*fZ2) + -0.598044251627f;
            pOut[222] = fPc * fC1;
            pOut[198] = fPc * fS1;
            fPa = fZ * ((fZ2* 53.7940721231f) + -5.56490401273f);
            pOut[252] = fPa * fC1;
            pOut[228] = fPa * fS1;
            fC0 = (fX*fC1) - (fY*fS1);
            fS0 = (fX*fS1) + (fY*fC1);
            fPa = -0.816077118838f;
            pOut[195] = fPa * fC0;
            pOut[169] = fPa * fS0;
            fPb = -4.39470978027f*fZ;
            pOut[223] = fPb * fC0;
            pOut[197] = fPb * fS0;
            fPc = (-17.6082433888f*fZ2) + 0.607180806512f;
            pOut[253] = fPc * fC0;
            pOut[227] = fPc * fS0;
            fC1 = (fX*fC0) - (fY*fS0);
            fS1 = (fX*fS0) + (fY*fC0);
            fPa = 0.830522083065f;
            pOut[224] = fPa * fC1;
            pOut[196] = fPa * fS1;
            fPb = 4.62415125663f*fZ;
            pOut[254] = fPb * fC1;
            pOut[226] = fPb * fS1;
            fC0 = (fX*fC1) - (fY*fS1);
            fS0 = (fX*fS1) + (fY*fC1);
            fPc = -0.844250650857f;
            pOut[255] = fPc * fC0;
            pOut[225] = fPc * fS0;

        }
    }


    /************************************************************************/
    /*                              Project BRDF to SH                      */
    /************************************************************************/
    /* @ref Laurent Belcour's Code for parallel.
     * @param todo:remove inheritance from std::thread
     */
    template<typename BSDFClass>
    struct BSDFProjection : std::thread
    {
        std::vector<std::vector<float3>> Bvec;

        BSDFProjection(int matrixWidth, int matrixHeight, const float3 *w, float *Ylm, int skip, int nthread, int nSample, const BSDFClass& bsdf) :
            std::thread(&BSDFProjection::run, this, w, Ylm, skip, nthread, nSample, bsdf, matrixWidth, matrixHeight)
        {
            // Note don't store nSamples in BSDFProjection but pass nSamples directly to the constructor
            // of std:thread(). std::thread() constructs first and runs background. 
            // Alternative: initialize in run().
        }

        void run(const float3 *w, float *Ylm, int skip, int nthread, int nSamples, const BSDFClass& bsdf, int matrixWidth, int matrixHeight)
        {
            int &BSDFMATRIX_HEIGHT = matrixHeight;
            int &BSDFMATRIX_WIDTH = matrixWidth;
            // Init
            std::vector<std::vector<float3>> vec(BSDFMATRIX_HEIGHT, std::vector<float3>(BSDFMATRIX_WIDTH, float3()));
            Bvec = vec;

            // Compute double spherical integral for BSDF matrix
            for (int osamp = skip; osamp < nSamples; osamp += nthread)
            {
                const float3 &wo = w[osamp];

                if (skip == 0)
                {
                    std::cout << "MC Progress: " << osamp << " / " << nSamples << "     \r";
                    std::cout.flush();
                }
                for (int isamp = 0; isamp < nSamples; ++isamp)
                {
                    const float3 &wi = w[isamp];
                    // Update BSDF matrix elements for sampled directions
                    float3 f = bsdf.Eval_f(wo, wi);
                    if (f.x > 0.f && f.y > 0.f && f.z > 0.f)
                    {
                        for (int y = 0; y < BSDFMATRIX_HEIGHT; ++y)
                        {
                            float *ylmVector_wi = &Ylm[isamp * BSDFMATRIX_WIDTH];
                            float *ylmVector_wo = &Ylm[osamp * BSDFMATRIX_WIDTH];

                            float3 factorPerRow_x_ = f * fabsf(wi.z) / (nSamples*pdfHemiSphere())*ylmVector_wi[y];
                            for (int x = 0; x < BSDFMATRIX_WIDTH; ++x)
                            {
                                Bvec[y][x] += 1.f * factorPerRow_x_ / (nSamples*pdfHemiSphere())*ylmVector_wo[x];
                            }
                        }
                    }
                }
            }
        }
    };

    /**
     * @brief Project BRDF to SH, support order 16 and order 10.
     * @param[in] SHIOrder  Maximum order for integration
     * @param[in] BSDFOrder Maximum order for BRDF
     */
    template<typename BSDFClass>
    void projectBRDFtoSH(const std::string &outputFilename, int SHIOrder, int BSDFOrder, const BSDFClass& bsdf)
    {
        assert((SHIOrder == 9 || SHIOrder == 15) && (BSDFOrder == 9 || BSDFOrder == 15) && BSDFOrder >= SHIOrder);
        const int BSDFMATRIX_WIDTH = (BSDFOrder + 1)*(BSDFOrder + 1);
        const int BSDFMATRIX_HEIGHT = (SHIOrder + 1)*(SHIOrder + 1);

        int nSamples = 7000;
        std::vector<float3> BSDFMatrix;
        BSDFMatrix.resize(BSDFMATRIX_HEIGHT*BSDFMATRIX_WIDTH, float3());

        // Precompute directions $\w{}$ and SH values for directions
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<unsigned int> dis_uint(0);

        // BSDFMATRIX_WIDTH is always the largest.
        float *Ylm = new float[BSDFMATRIX_WIDTH * nSamples];
        float3 *w = new float3[nSamples];
        uint32_t scramble[2] = { /*dis_uint(gen),dis_uint(gen)*/0,0 };
        for (int i = 0; i < nSamples; ++i)
        {
            float u[2];
            Sample02(i, scramble, u);
            w[i] = uniformSamplingHemisphere(u[0], u[1]);

            if (BSDFOrder == 15)
                SHEvalFast15(w[i], &Ylm[BSDFMATRIX_WIDTH*i]);
            else if (BSDFOrder == 9)
                SHEvalFast9(w[i], &Ylm[BSDFMATRIX_WIDTH*i]);
            else
                assert(false);
        }

        const int nbthreads = std::thread::hardware_concurrency();
        std::vector<BSDFProjection<BSDFClass>*> threads;
        for (int k = 0; k < nbthreads; ++k)
        {
            BSDFProjection<BSDFClass>* th = new BSDFProjection<BSDFClass>(BSDFMATRIX_WIDTH, BSDFMATRIX_HEIGHT, w, Ylm, k, nbthreads, nSamples, bsdf);
            threads.push_back(th);
        }

        for (BSDFProjection<BSDFClass>* th : threads)
        {
            th->join();
            // Gathering Result
            for (int y = 0; y < BSDFMATRIX_HEIGHT; ++y)
            {
                for (int x = 0; x < BSDFMATRIX_WIDTH; ++x)
                {
                    BSDFMatrix[y*BSDFMATRIX_WIDTH + x] += th->Bvec[y][x];
                }
            }
            delete th;
        }

        static_assert(std::is_trivial<float3>::value,
            "Can only serialize POD types with this function");

        std::ofstream ofs(outputFilename + ".dat", std::ios::out | std::ofstream::binary);
        (serialize<float3>(ofs, BSDFMatrix));

        delete[] Ylm;
        delete[] w;
    }
}