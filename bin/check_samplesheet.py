#!/usr/bin/env python

# This script is based on the example in "assets/samplesheet_test_illumina.csv"
import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/sentieon samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample,rg_id,rg_lb,rg_pl,fastq_1,fastq_2
    NA12878,NA12878.H88WKADXX.1.CGATGT-1,NA12878.U0a,ILLUMINA,SRR11140744_R1.fastq.gz,SRR11140744_R2.fastq.gz
    NA12878,NA12878.H88WKADXX.1.CGATGT-2,NA12878.U0a,ILLUMINA,SRR11140746_R1.fastq.gz,SRR11140746_R2.fastq.gz
    NA13000,NA13000.H88WKADXX.1.CGATGT-1,NA13000.U0a,ILLUMINA,SRR11140748_R1.fastq.gz,SRR11140748_R2.fastq.gz

    For an example see:
    <REPO_ROOT>/assets/samplesheet_test_illumina.csv
    """

    sample_mapping_dict = {}
    rg_ids = []
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 5
        HEADER = ["sample", "rg_id", "rg_lb", "rg_pl", "fastq_1", "fastq_2"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, rg_id, rg_lb, rg_pl, fastq_1, fastq_2 = lspl[: len(HEADER)]
            sample = sample.replace(" ", "_")
            if not sample:
                print_error("'sample' column entry has not been specified!", "Line", line)
            if not rg_id:
                print_error("'rg_id' column entry has not been specified!", "Line", line)
            if not rg_lb:
                print_error("'rg_lb' column entry has not been specified!", "Line", line)
            if not rg_pl:
                print_error("'rg_pl' column entry has not been specified!", "Line", line)

            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )

            ## Auto-detect paired-end/single-end
            sample_info = []  ## [single_end, fastq_1, fastq_2]
            if sample and fastq_1 and fastq_2:  ## Paired-end short reads
                sample_info = ["0", rg_id, rg_lb, rg_pl, fastq_1, fastq_2]
            elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
                sample_info = ["1", rg_id, rg_lb, rg_pl, fastq_1, fastq_2]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)
            
            ## Create sample mapping dictionary = { sample: [ single_end, rg_id, rg_lb, rg_pl, fastq_1, fastq_2 ] }
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample].append(sample_info)
            
            rg_ids.append(rg_id)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        if len(set(rg_ids)) != len(rg_ids):
            dup_ids = [x for n, x in enumerate(rg_ids) if x in rg_ids[:n]]
            print_error("'rg_id' column should be unique across all samples!", "rg_id", "{}".format(', '.join(dup_ids)))

        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample", "single_end", "rg_id", "rg_lb", "rg_pl", "fastq_1", "fastq_2"]) + "\n")
            for sample in sorted(sample_mapping_dict.keys()):

                ## Check that multiple runs of the same sample are of the same datatype
                if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
                    print_error("Multiple runs of a sample must be of the same datatype!", "Sample", "{}".format(sample))

                for idx, val in enumerate(sample_mapping_dict[sample]):
                    fout.write(",".join([sample] + val) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet", "{}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
