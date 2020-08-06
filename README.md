## Spherical Light Integration over Spherical Caps via Spherical Harmonics ##

Reference implementation of the paper *Spherical Light Integration over Spherical Caps via Spherical Harmonics*.

### Build ###

The project relies on cmake to construct and a modern C++ compiler with C++11 support.
Code for building:

	git clone https://github.com/Hearwindsaying/TripleSphere.git
	mkdir build
	cd build
	cmake ..
	make

#### Structure ####

The project consists of several headers, testing and selected precomputation data.
Core functions in `include/TripleSphere.h` directly compute sphere light SH coefficients on the fly, without depending on a concrete renderer.
Incoporate these functions into a renderer (either PRT interactive or offline physically based) should be easy, since there are no external dependencies for the implementation.

This code relies on templates to enable a simpler integration into other codebases. We provide examples of our wrappers in `tests/Test.hpp`. To use our code, you will have to define the following class wrappers:

   + `include/Serialize.h`,`include/SHUtils.h`,`include/VecType.h` provide utilities for reading/saving matrix, fundamental types, SH/ZH evaluations and helper functions for projecting BRDF to SH.
   + `include/TripleSphere.h` provides core implementation on computing light SH coefficients.
   + `tests/TestMC.cpp` tests our recurrence formulae and compares with Monte-Carlo method.
   + `tests/BRDF.h`,`tests/ProjectBRDF.cpp` give a sample code for computing BRDF matrix. It is straightforward to extend to other types of BRDF by adding classes with Eval_f(wo,wi) implementation.
   + `data/` has a few sample for BRDF matrix, weighting coefficients and lobe directions for ZHF. Note that for high orders (such as order 16), the matrix conditioning could have an impact on the results. Take care of numerical issues.

#### Usage ####

Preparation:

	  const float3 lobeDirs[2 * 15 + 1]{#include "15order_bluenoise_lobe_dirs.csv"};
	  const float alphaCoeff[256][31]{#include "15order_bluenoise_a_matrix.csv"};
	  const float diffuseCoeff[(15 + 1)*(15 + 1)]{#include "15order_diffuse.csv"};

Assuming it is time for shading BRDF with sphere lights, first figure out the parametrization:
	
	  // Convert light direction to TBN.
	  float3 L_TBN = normalize(frame.toLocal(normalize(sphereLight.center - p));
	  // Compute Light parameters, aperture angle for spherical light cap
      float sinSigma = sqrt(sphereLight.radius * sphereLight.radius / sqr_length(sphereLight.center - p));
      float cosSigma = sqrt(1 - sinSigma * sinSigma);

We could now compute the light SH coefficients:
	  
	  // Allocate space for storage
	  float3 ylmCoeff[(SHIOrder + 1)*(SHIOrder + 1)];
      for (int i = 0; i < (SHIOrder + 1)*(SHIOrder + 1); ++i)
          ylmCoeff[i] = make_float3(0.0f);

	  Sphere_computeCoeff<SHIOrder>(L_TBN, sinSigma, cosSigma, ylmCoeff, make_float3(sphereLight.intensity), lobeDirs, alphaCoeff);

After getting BRDF coefficients (through BRDF matrix or analytical formulae), the magic happens here via a simple dot product:

	  L += (diffuseCoeff[i] * ylmCoeff[i]);
	  

#### Acknowledgement ####

We would like to thanks the authors of [IntegralSH](https://github.com/belcour/IntegralSH) for the implementation of `Integrating Clipped Spherical Harmonics Expansions`.

#### Contact ####

Feel free to contact hearwindsaying@gmail.com if you have any questions.