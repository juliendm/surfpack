#!/bin/sh

SURFPACK=surfpack

rm -f poly1.sps test_data_save.spd test_data_load.spd

echo "Building and saving model"
${SURFPACK} ccd_save.spk
echo "Loading model"
${SURFPACK} ccd_load.spk

diffcount=`diff test_data_save.spd test_data_load.spd | wc -l`
echo "Diff count is: $diffcount"
