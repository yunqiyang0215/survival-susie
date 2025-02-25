#!/usr/bin/env dsc

simulate: sim.R
    censor_lvl: 0, 0.2, 0.4, 0.6, 0.8, 0.99
    allele_freq: 0.001, 0.01, 0.1
    n: 10000, 100000
    b1: 0.01, 0.1, 1
    $dat: dat

compute_bf: compute_bf.R
    # module input and variables
    dat: $dat
    # module output
    $result: result

DSC:
    run: simulate * compute_bf
    replicate: 50
    #R_libs: udr,mashr
    exec_path: bin
    output: dsc_result


