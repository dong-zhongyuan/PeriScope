#!/bin/bash
# Blender 4.2 on the CentOS-7 host (glibc 2.17): a copy of the binary patched (patchelf) to use the conda sysroot glibc 2.39
# loader, with libstdc++ from the creed-r env (blender_shim/lib); gconv modules must match that glibc.
SYS=/public/home/mengxl/.conda/pkgs/sysroot_linux-64-2.39-he3f20f0_6/x86_64-conda-linux-gnu/sysroot
export GCONV_PATH=$SYS/lib64/gconv LC_ALL=C LANG=C
exec /public/home/mengxl/dzy/tools/blender-4.2.23-linux-x64/blender_glibc239 "$@"
