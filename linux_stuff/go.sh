#!/bin/sh

cd build
rm Makefile
rm lima


cmake ../
make
./lima


