#!/bin/bash

#wget https://raw.githubusercontent.com/kward/shunit2/master/shunit2

TEST_DIR=$(cd "$(dirname "$0")" && pwd)/test

. test/test_alignment_minimap2.sh
. test/test_alignment_ultra.sh

testLongReadMinimap2DNA() {
  runLongReadMinimap2DNA $TEST_DIR
}

testLongReadMinimap2RNA() {
  runLongReadMinimap2RNA $TEST_DIR
}

testLongReadULTRA() {
  runLongReadULTRA $TEST_DIR
}

. ./shunit2