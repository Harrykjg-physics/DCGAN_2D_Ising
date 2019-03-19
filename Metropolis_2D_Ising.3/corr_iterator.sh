#!/bin/bash
for T in $(seq 1.9 0.1 2.6 )
    do 
        python 220A-Ising.py $T
    done
