#!/bin/bash --login
set -e

# activate conda environment and let the following process take over
conda activate phyluce
exec "$@"