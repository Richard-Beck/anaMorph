#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 1 ]; then
  echo "Usage: $0 path/to/script.R [args ...]" >&2
  exit 1
fi

SCRIPT_PATH=$1
shift

CONTAINER_URI=${CONTAINER_URI:-"docker://dockerhub.moffitt.org/hpc/rocker-rstudio:4.4.2"}
USER_NAME=${USER_NAME:-$(id -un)}
SCRIPT_PATH=$(readlink -f "$SCRIPT_PATH")
SCRIPT_DIR=$(dirname "$SCRIPT_PATH")
WORKDIR=$(pwd)

mkdir -p "$WORKDIR/.tmp" "$WORKDIR/.cache"

export APPTAINERENV_XDG_CACHE_HOME=${APPTAINERENV_XDG_CACHE_HOME:-"$WORKDIR/.cache"}
export APPTAINERENV_TMPDIR=${APPTAINERENV_TMPDIR:-"$WORKDIR/.tmp"}

BIND_PATHS="/home/$USER_NAME,/share,/etc/passwd,/etc/group,$SCRIPT_DIR,$WORKDIR,$APPTAINERENV_XDG_CACHE_HOME,$APPTAINERENV_TMPDIR"
BINDS=${BINDS:-"-B $BIND_PATHS"}

exec apptainer exec $BINDS "$CONTAINER_URI" Rscript "$SCRIPT_PATH" "$@"