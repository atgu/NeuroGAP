## Copyright Broad Institute, 2018
##
## This WDL pipeline unmaps a BAM file
##
## Requirements/expectations :
## - A tab delimited file list with two columns: (1) a BAM file ; (2) and BAM name (str) without the extension
##
## Output :
## - An unmapped BAM file
##
## Cromwell version support
## - Successfully tested on v36
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
## For program versions, see docker containers.
##
## LICENSING :
## This script is released under the WDL source code license (BSD-3) (see LICENSE in
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script.

## Notes
# https://support.terra.bio/hc/en-us/articles/360037128072--2-howto-Write-a-simple-multi-step-workflow

# WORKFLOW DEFINITION

workflow BamtoUbam {
  File ref_fasta
  File ref_index
  #File ref_dict
  File input_samples_list

  Array[Array[File]] input_samples = read_tsv(input_samples_list)

  String? gatk_docker_override
  String gatk_docker = select_first([gatk_docker_override, "broadinstitute/gatk:latest"])
  String? gitc_docker_override
  String gitc_docker = select_first([gitc_docker_override, "broadinstitute/genomes-in-the-cloud:2.3.1-1500064817"])
  String? samtools_path_override
  String samtools_path = select_first([samtools_path_override, "samtools"])

  Int? preemptible_attempts

  scatter (sample in input_samples) {

    String inputfile = sample[0]
    String output_name = sample[1]
    # is the input a cram file?
    Boolean is_cram = sub(basename(inputfile), ".*\\.", "") == "cram"

    #### IF IT IS A CRAM FILE ####
    if (is_cram) {
      call CramToBamTask {
        input:
          input_cram = inputfile,
          bam_out_name = output_name,
          #ref_dict = ref_dict,
          ref_fasta = ref_fasta,
          ref_index = ref_index,
          docker = gitc_docker,
          samtools_path = samtools_path
      }

      call BamToUbamTask as crams {
        input:
          input_bam = CramToBamTask.output_bam,
          ubam_out_name = output_name,
          docker = gatk_docker,
          preemptible_attempts = preemptible_attempts
      }
    }

    #### IF IT IS A BAM FILE ####
    if(!is_cram){
      call BamToUbamTask as bams {
        input:
          input_bam = inputfile,
          ubam_out_name = output_name,
          docker = gatk_docker,
          preemptible_attempts = preemptible_attempts
      }
    }
  }
}


# Converts bam to ubam files
task BamToUbamTask {
    # Command parameters
    File input_bam
    String ubam_out_name

    # Runtime parameters
    String docker
    Int? machine_mem_gb
    Int machine_mem = select_first([machine_mem_gb, 18]) # 7
    Int? disk_space_gb
    Int? preemptible_attempts

    Int command_mem_gb = machine_mem - 1
    Int disk_size = ceil( size(input_bam, "GB") * 2) + 100 # 20

    command {
       echo ${input_bam}
       echo ${ubam_out_name}
       gatk RevertSam \
       -I ${input_bam} \
       -O ${ubam_out_name}.ubam \
       -VALIDATION_STRINGENCY LENIENT
    }

    runtime {
        docker: docker
        memory: machine_mem + " GB"
        disks: "local-disk " + select_first([disk_space_gb, disk_size]) + " HDD"
        preemptible: select_first([preemptible_attempts, 3])
    }

    output {
        File output_ubam = "${ubam_out_name}.ubam"
    }
}


# converts CRAM to BAM
task CramToBamTask {
  # Command parameters
  File ref_fasta
  File ref_index
  #File ref_dict
  File input_cram
  String bam_out_name

  # Runtime parameters
  String docker
  Int? machine_mem_gb
  Int? disk_space_gb
  Boolean use_ssd = false
  Int? preemptible_attempts
  String samtools_path

  Float output_bam_size = size(input_cram, "GB") / 0.40
  Float ref_size = size(ref_fasta, "GB") + size(ref_index, "GB")
  Int disk_size = ceil(size(input_cram, "GB") + output_bam_size + ref_size) + 25

  command {
    set -e
    set -o pipefail

    ${samtools_path} view -h -T ${ref_fasta} ${input_cram} |
    ${samtools_path} view -b -o ${bam_out_name}.bam
  }
  runtime {
    docker: docker
    memory: select_first([machine_mem_gb, 15]) + " GB"
    disks: "local-disk " + select_first([disk_space_gb, disk_size]) + if use_ssd then " SSD" else " HDD"
    preemptible: select_first([preemptible_attempts, 3])
 }
  output {
    File output_bam = "${bam_out_name}.bam"
  }
}
