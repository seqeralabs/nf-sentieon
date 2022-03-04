#!/bin/bash
set -eu
export SENTIEON_LICENSE=$(mktemp)
echo -e "$SENTIEON_LICENSE_BASE64" | base64 -d > $SENTIEON_LICENSE
