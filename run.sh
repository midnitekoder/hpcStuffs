#!/bin/bash
echo "Tintin"
echo $1
qsub -l gpu=1 -b y -cwd ./a.out $1
