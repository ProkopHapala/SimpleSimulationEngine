#!/bin/bash

# ====== Paths

args=""

#name=test_DiffractShader    # diffraction shader demo
#name=test_SphereShader      # raycasted sphere impostor demo
#name=test_ShaderDepth       # depth / shadow map testing
#name=test_DrawOGL3          # general DrawOGL3 tests
#name=test_Instances         # instanced meshes demo
#name=test_Atoms             # atom impostor spheres
#name=test_Sprites           # billboard sprites
#name=test_LandScape        # landscape rendering
#name=test_Vegetation       # vegetation rendering
#name=test_Horizont         # sky/horizon rendering
#name=test_VAOs             # VAO usage demo
#name=test_PatchesOGL3      # patch rendering
#name=test_StencilTextures  # stencil texture demo
#name=test_VolumetricTexture# volumetric textures
#name=test_RenderToTexture  # render-to-texture pipeline
#name=test_AntiAliasing     # anti-aliasing demo
#name=test_GeometryShader   # geometry shader demo
#name=test_Tubes            # tube rendering
#name=test_OrbitalRayMarch  # orbital ray-marching demo
#name=test_Texture          # basic textured rendering
#name=test_SSAO             # SSAO demo
#name=test_ScreenOGL3       # ScreenSDL2OGL3 demo
#name=test_MeshOGL3         # mesh rendering demo
name=meshviewer            # generic mesh viewer using MeshRenderOGL3
args="common_resources/my_mesh.obj"


dir=../../cpp/Build/sketches_SDL/OGL3
ln -s ../../cpp/apps/OrbitalWar/data
ln -s ../../cpp/common_resources

# ====== Multiprocesing
ncpu=`nproc`
ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
echo "compile using ncpu="$ncpu
OMP_NUM_THREADS=$ncpu
export OMP_NUM_THREADS

# ====== Compilation
wd=`pwd`
cd $dir
pwd
rm $name
make -j$ncpu $name   # 2>$wd/compile_err.log
cd $wd
rm $name.x
ln -s $dir/$name ./$name.x

# ====== ASan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

# ====== RUN

echo "test_OGL3.sh args: " $args
./$name.x $args
