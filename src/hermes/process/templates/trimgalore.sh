#!/bin/env bash

set -eu

trimgalore \
  --paired {{ .R1 }} {{ .R2 }} \
  --cores {{ .Threads }} \
  --output_dir output_trimgalore/

