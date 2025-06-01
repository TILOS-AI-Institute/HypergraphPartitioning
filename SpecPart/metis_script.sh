#!/bin/bash
./SpecPart/gpmetis $1 $2 -ptype=rb -ufactor=$3 -seed=$4 -dbglvl=0 > $5
