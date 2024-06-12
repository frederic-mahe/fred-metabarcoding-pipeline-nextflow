#!/usr/bin/env nextflow

params.forward_primer = "CCAGCASCYGCGGTAATTCC"
params.reverse_primer = "ACTTTCGTTCTTGATYRA"  // ARYTAGTTCTTGCTTTCA, should be TYRATCAAGAACGAAAGT
params.reverse_primer_revcomp = "TYRATCAAGAACGAAAGT"
params.fastq_folder = "data"
params.fastq_pattern = "/*_{1,2}.fastq.gz"
params.fastq_encoding = 33
params.threads = 4

process generate_test_data_urls {
    output:
    path "test_data_urls.list"

    shell:
    '''
    #!/bin/bash

    URL="https://github.com/frederic-mahe/BIO9905MERG1_vsearch_swarm_pipeline/raw/main/data"

    (echo "${URL}/MD5SUM"
     for SAMPLE in {B,L}{010..100..10} ; do
         for READ in 1 2 ; do
             echo "${URL}/${SAMPLE}_1_${READ}.fastq.gz"
         done
     done
    ) > test_data_urls.list
    '''
}


process download_list_of_urls {
    input:
    path "urls"

    publishDir params.fastq_folder

    output:
    path "*.fastq.gz"

    shell:
    '''
    #!/bin/bash

    wget --continue --input-file="!{urls}"
    '''
}


process merge_fastq_pairs {
    input:
    tuple val(sampleId), path(fastq_pair)

    output:
    val sampleId
    path "merged_fastq"

    shell:
    '''
    #!/bin/bash

    vsearch \
        --fastq_mergepairs !{fastq_pair[0]} \
        --reverse !{fastq_pair[1]} \
        --threads !{params.threads} \
        --fastq_ascii !{params.fastq_encoding} \
        --fastq_allowmergestagger \
        --quiet \
        --fastqout merged_fastq
    '''
}


process revcomp {
    // reverse-complement a DNA/RNA IUPAC string
    input:
    val "reverse_primer"
    
    output:
    stdout

    shell:
    '''
    #!/bin/bash

    [[ -z !{reverse_primer} ]] && { echo "error: empty string" ; exit 1 ; }
    nucleotides="acgturykmbdhvswACGTURYKMBDHVSW"
    complements="tgcaayrmkvhdbswTGCAAYRMKVHDBSW"
    tr "${nucleotides}" "${complements}" <<< !{reverse_primer} | rev
    '''
}


process trim_primers {
    // search forward primer in both normal and revcomp: now all reads
    // are in the same orientation. Matching leftmost is the default.
    input:
    val sampleId
    path merged_fastq

    output:
    val sampleId
    path trimmed_fastq

    shell:
    '''
    #!/bin/bash

    MIN_F=$(( !{params.forward_primer.length()} * 2 / 3 ))  # match is >= 2/3 of primer length
    MIN_R=$(( !{params.reverse_primer.length()} * 2 / 3 ))
    cutadapt \
        --revcomp \
        --front "!{params.forward_primer};rightmost" \
        --overlap "${MIN_F}" !{merged_fastq} | \
        cutadapt \
            --adapter "!{params.reverse_primer_revcomp}" \
            --overlap "${MIN_R}" \
            --max-n 0 - > trimmed_fastq &
    '''
}


workflow {

    // collect test data
    generate_test_data_urls |
        download_list_of_urls

    // reverse-complement reverse primer
    // params.reverse_primer_revcomp = revcomp(params.reverse_primer)
    
    // merge, trim
    Channel.fromFilePairs(params.fastq_folder + params.fastq_pattern) |
        merge_fastq_pairs |
        trim_primers
}