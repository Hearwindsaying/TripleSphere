#include "BRDF.h"

#include <SHUtils.h>

using namespace ts;

int main()
{
    RoughMetal roughmetal{ 0.3, make_float3(0.143119f, 0.374957f, 1.442479f), make_float3(3.983160f, 2.385721f, 1.603215f) };
    projectBRDFtoSH<RoughMetal>("roughgold", 15, 15, roughmetal);

    return 0;
}