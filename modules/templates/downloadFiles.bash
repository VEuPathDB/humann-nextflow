#!/usr/bin/env bash

set -euo pipefail
fasterq-dump --split-3 ${id}
