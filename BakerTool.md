# Introduction #

The baker tool is a command-line tool that allows you to capture detail attributes from high resolution meshes and bake them into textures, so that they can be applied to lower resolution surfaces.

You can browse the source code of the baker tool here:

http://code.google.com/p/nvidia-mesh-tools/source/browse/#svn/trunk/src/nvmesh/tools/baker

The NVIDIA baker tool supports some unique features not available in other similar tools:

  * Support for different kind of base surfaces, including various subdivision surface approximations, such as Bezier ACC, Gregory ACC, and Pm patches.
  * Capture details using dual parameterizations instead of raycasting. This imposes some constraints on the input assets, but produces much higher quality results.
  * Support for vector displacement maps.


# Command-line Options #

## General Options ##

  * **--width** - Width of the output textures. _Defaults to 1024_.
  * **--height** - Height of the output textures. _Defaults to 1024_.
  * **--threads** - Number of threads to be used to compute ambient occlusion. _Default to 0, which chooses a number equal to the number of processors_.
  * **--output** - File name prefix for the generated textures. _Defaults to "output"_.

## Base Pass Options ##

  * **--lomesh** - File name of the base mesh.
  * **--losub** - Subdivision levels applied to the base mesh. _Defaults to 0_.
  * **--lomode** - Mode of the base mesh:
    * 0 = Bezier ACC surface. _(default)_
    * 1 = Gregory ACC surface.
    * 2 = Triangular mesh.
  * **--disp3D** - Generate vector displacements.
  * **--tangentSpace** - Generate normals and vector displacements in tangent space.

## Detailed Pass Options ##

  * **--himesh** - File name of the detailed mesh.
  * **--hisub** - Subdivision levels applied to the detailed mesh. _Defaults to 0_.
  * **--hidisp** - Displacement map for detailed mesh.
  * **--hinmap** - Normal map for detailed mesh.
  * **--hisupersampling** - Enable supersampling. _Disabled by default_.
  * **--hipos** - Output position map.
  * **--occlusion** - Output occlusion map.
  * **--bentnormal** - Output bent normal map.
  * **--rays** - Number of rays to compute ambient occlusion. _Defaults to 64_.
  * **--maxdistratio** - Maximum occlusion distance with respect to the mesh size. _Defaults to 0.25_.
  * **--ilmocclusion** - Use binary occlusion instead of attenuated visibility.
  * **--hicachedir** - Directory for saving and reading the cached data.
  * **--rebuildcache** - Load original data and recompute cache.

# Implementation #

The baker tool generally loads two meshes, a high resolution _detailed_ mesh and a low resolution _base_ mesh. These two meshes must have the same parameterization. That is, if you rasterize the two meshes in texture space, the output of the two meshes should overlap almost exactly.

The baker tool, relies on this dual parameterization in order to transfer attributes from one mesh to the other.

This is done in two passes. First, in the detailed mesh pass the detailed mesh is loaded and rasterized into a geometry image. Then, in the base mesh pass the base mesh is loaded, rasterized into another geometry image. During the second step, the two geometry images are compared and some attributes are transferred or transformed from the detailed geometry image to the base geometry image. Finally, the attributes stored in both of the geometry images are filtered and stored to disk.


## Geometry Image ##

The [GeometryImage class](http://code.google.com/p/nvidia-mesh-tools/source/browse/trunk/src/nvmesh/tools/baker/GeometryImage.h#14) is just a wrapper over a floating point image. It simply provides some constructors to allocate an image with the appropriate number of channels, and some accessors to easily access these channels.

Two geometry images are used in the baker tool. The detailed mesh pass generates a geometry image that contains position and normal attributes that are both in object space. The base mesh pass generates a geometry image with position and normal attributes too, but in addition to that it also generates displacement maps and tangent space normal maps.

The geometry image simply provides a mechanism to store these attributes, apply some filtering to them and save them to disk.


## Base Surface ##

A base surface is not necessarily a triangular mesh. Triangular meshes are supported, several methods to evaluate the positions, normals and tangent space of the triangular mesh are available. In addition to triangular meshes the baker tool also has support for approximate subdivision surfaces. This is important, because when sampling displacements the distance from the base mesh to the detailed mesh must take the evaluation method into account.

In order to provide support for different surface types, the baker provides a [BaseSurface interface](http://code.google.com/p/nvidia-mesh-tools/source/browse/trunk/src/nvmesh/tools/baker/BaseSurface.h#21) with multiple implementations.


## Rasterization ##

In order to store mesh attributes in texture space the baker uses a [simple rasterizer](http://code.google.com/p/nvidia-mesh-tools/source/browse/#svn/trunk/src/nvmesh/raster). The rasterizer supports rendering of triangles and quads using linear or bilinear interpolation (you can approximate non linear interpolation using tessellation/subdivision). The texture parameterization should be unique, so that polygons do not overlap in texture space. However, if there's any overlap, contributions from both sources can be added together and averaged.

In order to provide high quality anti-aliasing the rasterizer performs analytic estimation of coverage by clipping the triangle to each texel in order to evaluate the area of their intersection. In addition to that, it can optionally perform supersampling although that's usually not needed.

The rasterizer does not support shaders, but instead you can provide a callback function that is invoked on every texel and that allows you to perform arbitrary operations, like sampling surfaces at certain parametric locations or writing the result to textures. The rasterizer does not output anything, it simply invokes the callback.


## Detail Pass ##

During the detail mesh pass the detail mesh is loaded and rasterized. On every texel the detailed triangular mesh is evaluated and the surface attributes are stored in the base geometry image. This is done by the [GeometrySampler](http://code.google.com/p/nvidia-mesh-tools/source/browse/trunk/src/nvmesh/tools/baker/Samplers.h#19), which is a class that provides a callback to the rasterizer. After rasterization, the geometry image is filtered to extrapolate the sampled texels to the empty ones. The [DetailMeshPass class](http://code.google.com/p/nvidia-mesh-tools/source/browse/trunk/src/nvmesh/tools/baker/DetailedMeshPass.h) coordinates all these steps.


## Ambient Occlusion ##

After the detail mesh pass, ambient occlusion can be evaluated, this is done by raytracing rays with a montecarlo distribution for every texel of the geometry image. Raytracing is performed using a [kd-tree acceleration structure](http://code.google.com/p/nvidia-mesh-tools/source/browse/#svn/trunk/src/nvmesh/kdtree). The geometry image is [divided into tiles](http://code.google.com/p/nvidia-mesh-tools/source/browse/trunk/src/nvmesh/tools/baker/TiledTask.h) and multiple threads are used to process these tiles.


## Base Pass ##

The base mesh pass takes a base mesh and the geometry image generated during the detail mesh pass. Then it rasterizes the base surface using the [DisplacementSampler](http://code.google.com/p/nvidia-mesh-tools/source/browse/trunk/src/nvmesh/tools/baker/Samplers.h#44). Despite its name this sampler does more than computing displacements, it also transforms normals from object space to tangent space. The [BaseMeshPass](http://code.google.com/p/nvidia-mesh-tools/source/browse/trunk/src/nvmesh/tools/baker/BaseMeshPass.h) coordinates these steps.


## Filtering ##

Several kind of filters are supported. I generally use several iterations of a [quadratic extrapolation filter](http://code.google.com/p/nvidia-mesh-tools/source/browse/trunk/src/nvimage/HoleFilling.cpp#465) followed by a [push-pull filter](http://code.google.com/p/nvidia-mesh-tools/source/browse/trunk/src/nvimage/HoleFilling.cpp#388) or a [Voronoi filling filter](http://code.google.com/p/nvidia-mesh-tools/source/browse/trunk/src/nvimage/HoleFilling.cpp#153).