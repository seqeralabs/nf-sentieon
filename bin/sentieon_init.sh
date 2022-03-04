#!/bin/bash
#
# Sentieon initialization script
# This script takes as input the name of a environment 
# variable holding the Sentieon license encoded as Base64 text
# and save it into a temporary file. The file path is set 
# into the variable with name `SENTIEON_LICENSE` as expected 
# by the sentieon binary
#  
set -eu
export SENTIEON_LICENSE=$(mktemp)
echo -e "${!1}" | base64 -d > $SENTIEON_LICENSE
