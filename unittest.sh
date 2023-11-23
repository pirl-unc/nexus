#!/bin/bash

#wget https://raw.githubusercontent.com/kward/shunit2/master/shunit2

TEST_DIR=$(cd "$(dirname "$0")" && pwd)/test

. test/test_alignment_minimap2.sh

testLongReadMinimap2() {
  runLongReadMinimap2 $TEST_DIR
}

. ./shunit2