#!/usr/bin/env dsc

simulate: sim.R
    p: 1000
    cor_type: "real", "independent"
    num_effect: 0, 1, 2, 3
    prior_variance: 0.1
    censor_lvl: 0, 0.2, 0.4, 0.6, 0.8, 0.99
    $dat: dat
    $b: b
    $effect_indx: indx_snp
    $is_effect: is_effect

susie: fit_susie.R
    # module input and variables
    dat: $dat
    is_effect: $is_effect
    # module output
    $effect_estimate: effect_estimate
    $pip: pip
    $cs: cs
    $iter: iter

bvsnlp: fit_bvsnlp.R
    # module input and variables
    dat: $dat
    is_effect: $is_effect
    # module output
    $effect_estimate: effect_estimate
    $pip: pip

svb: fit_survivalsvb.R
    # module input and variables
    dat: $dat
    is_effect: $is_effect
    # module output
    $effect_estimate: effect_estimate
    $pip: pip

r2b: fit_r2b.R
    # module input and variables
    dat: $dat
    is_effect: $is_effect
    # module output
    $effect_estimate: effect_estimate
    $pip: pip


susie_rss: fit_susierss.R
    # module input and variables
    dat: $dat
    is_effect: $is_effect
    # module output
    $effect_estimate: effect_estimate
    $pip: pip
    $cs: cs
    $iter: iter

DSC:
    define:
      evaluate: susie, bvsnlp, svb, r2b, susie_rss
    run: simulate * evaluate
    replicate: 10
    #R_libs:
    exec_path: bin
    output: dsc_result

