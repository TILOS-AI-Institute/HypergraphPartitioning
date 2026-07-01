#!/bin/bash
# gpmetis wrapper. Override the binary with the GPMETIS environment variable.
GPMETIS="${GPMETIS:-/home/fetzfs_projects/SpecPart/src/gpmetis}"
"$GPMETIS" "$1" "$2" -ptype=rb -ufactor="$3" -seed="$4" -dbglvl=0 > "$5"