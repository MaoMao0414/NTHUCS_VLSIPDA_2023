#!/bin/bash
make clean
make
../bin/hw3 ../testcase/public$1.txt ../output/public$1.floorplan
../verifier/verify ../testcase/public$1.txt ../output/public$1.floorplan