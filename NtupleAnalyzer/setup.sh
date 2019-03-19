#!/bin/bash

if [ $MUONHLT_ANALYZER_PATH ]; then
    echo "MUONHLT_ANALYZER_PATH is already defined: use a clean shell!"
    return 1
fi

export MUONHLT_ANALYZER_PATH=$(pwd)

# -- root setup -- #
export ROOT_INCLUDE_PATH=${MUONHLT_ANALYZER_PATH}:${ROOT_INCLUDE_PATH}
export PYTHONPATH=${MUONHLT_ANALYZER_PATH}:${PYTHONPATH}

echo "Enviornment variables"
echo "MUONHLT_ANALYZER_PATH: "$MUONHLT_ANALYZER_PATH
echo "ROOT_INCLUDE_PATH:     "$ROOT_INCLUDE_PATH
echo "PYTHONPATH:            "$PYTHONPATH
echo "Setup for Muon HLT analyzer is done."