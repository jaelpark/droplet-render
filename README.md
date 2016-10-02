# Droplet-render
Volumetric cloud modeling and rendering for Blender.

## Features
- <b>Node based approach to implicit surface and fog volume modeling.</b> Initial properties for clouds, fog and velocity fields are modeled with Blender. Surfaces and particle systems are then seamlessly exported, after which they are converted, modified and stored sparsely with OpenVDB grids. Nodes for typical surface and volume operations, such as displacement, composition and advection are provided.

- <b>Physically based cloud rendering.</b> A render engine specifically designed for volumetric rendering. Hybrid code for simultaneous distance field and fog volume sampling. Full global illumination with sky lighting. Mie phase sampling with spectral approximation. Every step from modeling to rendering is highly multithreaded and vectorized. Efficient parallel acceleration is achieved using Intel's TBB library.

- <b>Full blender integration with seamless data export/import.</b>
No temporary files, render targets or additional export steps required. Works with official Blender builds starting from 2.77.

Being relatively new the project is highly WIP. To achieve the desired high realism current focus is to improve the tools for cloud modeling. This also involves lots of time consuming testing with different parameters and potential approaches. Documentation and examples will be added later, when the project is considered ready for public testing.

Improvements to the rendering side will follow the modeling progress. Current implementation is a brute-force Monte-Carlo path tracer with basic hierarchial and MIS optimizations. The renderer is currently CPU-only. However, the code is purposely compact and independend of OpenVDB or any other libraries, and as such readily portable to GPU.

Note that Droplet does not apply any tonemapping or gamma correction. For correct results these will always have to be done manually afterwards. It is advisable to do all the testing with blender started from a terminal. Uncached displacement and advection steps may depending on the resolution and complexity take a significant portion of time, and before rendering starts, the engine does not provide any additional feedback besides stdout.

## Installation
Only Linux operating systems on x86_64 are currently supported. It should be possible to also build and work on Windows, although testing has been less active recently.
1. Install or build the required dependencies:
    - openvdb
    - tbb
    - python3.5 + numpy
    - cmake
2. Create a build directory, for example:
    ```sh
    $ mkdir build-release
    $ cd build-release
    ```
3. Droplet uses Cmake to configure itself. By default a release version with SSE3 support will be built. SSE4.2 can be enabled by passing an additional defition:
    ```sh
    $ cmake -DUSE_SSE4=ON ..
    ```
4. Build the project by typing
    ```sh
    $ make -j
    ```
5. Copy or symlink the produced libdroplet.so to Python library directory, usually /usr/lib/python3.5/site-packages. Copy the contents of 'addon' to Blender addons location by first creating a directory for it, such as 'droplet_render'. Enable the addon from Blender's settings menu.
