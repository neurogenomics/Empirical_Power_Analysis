#!/bin/bash -ue
mkdir -p Datalink/Full_Ground_Truth/Astro Datalink/Full_Ground_Truth/Micro Datalink/Full_Ground_Truth/Oligo
touch Datalink/Full_Ground_Truth/Astro/astro_gt_stub.csv
touch Datalink/Full_Ground_Truth/Micro/micro_gt_stub.csv
touch Datalink/Full_Ground_Truth/Oligo/oligo_gt_stub.csv
echo "stub" > versions.yml
