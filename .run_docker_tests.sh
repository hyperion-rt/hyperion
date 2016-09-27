#!/bin/bash

GENERATE=0
REFDIR=reference

while getopts "g" opt; do
    case "$opt" in
    g)  GENERATE=1
        ;;
    esac
done

docker info

cat << EOF | docker run -i \
                        -v ${PWD}:/hyperion_src \
                        -a stdin -a stdout -a stderr \
                        astrofrog/hyperion-ci:1.2 \
                        bash || exit $?

set -x

cd /hyperion_src
./configure
make serial
make install
python setup.py install

if [ $GENERATE == 0 ]; then
  python setup.py test --enable-bit-level-tests;
else
  mkdir -p ${REFDIR};
  python setup.py test --enable-bit-level-tests --generate-reference=${REFDIR};
fi

EOF
