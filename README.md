
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)

> THIS IS A PROOF-OF-CONCEPT REPOSITORY THAT IS UNDER ACTIVE DEVELOPMENT. SYNTAX, ORGANISATION AND LAYOUT MAY CHANGE WITHOUT NOTICE!

## Introduction

**nf-sentieon** is a simple, proof-of-concept pipeline to run the [Sentieon](https://www.sentieon.com/) suite of tools.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker containers making installation trivial and results highly reproducible.

Sentieon is licensed software that you can obtain for a [free trial](https://www.sentieon.com/home/free-trial/) via the website. This pipeline has only been tested using AWS Batch and doesn't come with any bundled containers.

## Credits

nf-sentieon was originally written by [Harshil Patel](https://github.com/drpatelh) and [Graham Wright](https://github.com/gwright99), [Seqera Labs](https://seqera.io/) in collaboration with [Don Freed](https://github.com/DonFreed), [Sentieon](https://www.sentieon.com/).

## Citations

The nf-core pipeline template was used to create the skeleton of this pipeline but there are no plans to contribute it to nf-core at this point.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
