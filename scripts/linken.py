#!/usr/bin/env python
# This is linken container entry script

import os

build_date = os.environ.get("BUILD_DATE", "unknown")
git_hash = os.environ.get("GIT_HASH", "unknown")


## -----------------------------------------------------------------------------
linken_banner = rf"""

  Linken Container (https://github.com/aixnr/linken)
  ==================================================

  This container image was generated by Docker and converted
  to Singularity Image Format using Apptainer.

  There are 3 subprograms available:

    - lunar, a utility script for preparing sequencing workflows
      $ linken lunar
    - hermes, for running the mapping and variant calling
      $ linken hermes -h
    - eevee, utility for analyzing mapped reads
      $ linken eevee
    - ont-vcall, a utility script for ONT sequencing data
      $ linken ont-vcall

  Build date: {build_date}
  Commit hash: {git_hash}
"""


## -----------------------------------------------------------------------------
# main
if __name__ == "__main__":
    print(linken_banner)
