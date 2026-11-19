#!/bin/bash
# usage: run_blender.sh <dir> <stem>   renders <dir>/layout_<stem>.json -> <dir>/<stem>_art.png (+ .rgba for R)
set -e; S=$(cd "$(dirname "$0")" && pwd); cd "$1"
"$S/blender.sh" -b --python "$S/blender_render.py" -- "layout_$2.json" "$PWD/$2_art.png" 2>&1 | grep -E "RENDER_DONE|Error|Traceback|line [0-9]+" | head -8
. "$S/env.sh"; "$FIGV2_PY" "$S/to_rgba.py" "$2_art.png" "$2_art.rgba"
