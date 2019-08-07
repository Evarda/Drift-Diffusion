#!/bin/bash

make

mkdir Data

mkdir Data/Equilibrium

mkdir Data/Potential

mkdir Data/Diffusion

mkdir Data/Mu

mkdir Data/n

mkdir Data/p

mkdir Data/Non-Equilibrium

./driftDiffusion

mkdir _Data

mv Data _Data/$1_$2