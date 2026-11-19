# figures_v2 environment: source before any render. Uses the pd-v3 runtime read-only (R, R library, Arial files);
# every cache and temp file goes under _system so nothing is written outside figures_v2_20261119.
SYS=/public/home/mengxl/dzy/pd-v3/papers/figures_v2_20261119/_system
export FIGV2_SYS=$SYS
export FIGV2_ROOT=$(dirname $SYS)
export PATH=/public/home/mengxl/.conda/envs/creed-r/bin:/usr/bin:$PATH
export R_LIBS_USER=/public/home/mengxl/dzy/pd-v3_runtime/R-library
export FONTCONFIG_FILE=$SYS/fonts.conf
export TMPDIR=$SYS/.tmp
export XDG_CACHE_HOME=$SYS/.cache
export PYTHONDONTWRITEBYTECODE=1
export LC_ALL=en_US.utf8 LANG=en_US.utf8
export FIGV2_RSCRIPT=/public/home/mengxl/.conda/envs/creed-r/bin/Rscript
export FIGV2_PY=/public/home/mengxl/dzy/pd-v3_runtime/bin/python
mkdir -p $TMPDIR $XDG_CACHE_HOME $SYS/.fontcache
