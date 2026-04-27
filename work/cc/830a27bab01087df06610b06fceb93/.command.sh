#!/bin/bash -ue
mkdir -p Datalink/Cartesian_Downsampling/Astro Datalink/Cartesian_Downsampling/Micro Datalink/Cartesian_Downsampling/Oligo
touch Datalink/Cartesian_Downsampling/Astro/astro_cartesian_stub.csv
touch Datalink/Cartesian_Downsampling/Micro/micro_cartesian_stub.csv
touch Datalink/Cartesian_Downsampling/Oligo/oligo_cartesian_stub.csv
echo "stub" > versions.yml
