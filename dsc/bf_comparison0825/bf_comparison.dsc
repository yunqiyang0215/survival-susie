#!/usr/bin/env dsc

simulate: sim.R
    censor_lvl: 0
    allele_freq: 0.1
    n: 5000, 500000
    b1: 1
    binary: "bin", "nonbin"
    $dat: dat
    $bin: binary

compute_bf: compute_bf.R
    # module input and variables
    dat: $dat
    binary: $bin
    flip: "nonflip", "flip"
    # module output
    $result: result

DSC:
    run: simulate * compute_bf
    replicate: 50
    #R_libs: udr,mashr
    exec_path: bin
    output: dsc_result


