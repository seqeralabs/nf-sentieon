#!/bin/bash

set -eu
export SENTIEON_LICENSE=$(mktemp)
echo -e "$sentieon_license_text" > $SENTIEON_LICENSE
