#!/bin/bash
#
# Drive the HiPerGator sweeps from this computer. Run from anywhere.
#
#   simulations/hipergator/hpg.sh test          upload, then the three-task smoke test
#   simulations/hipergator/hpg.sh run           upload, submit; keeps existing blocks
#   simulations/hipergator/hpg.sh run --fresh   upload, DELETE all blocks, submit
#   simulations/hipergator/hpg.sh width         upload, submit the width cells; --fresh as above
#   simulations/hipergator/hpg.sh status        queue, block counts, and the pooling log
#   simulations/hipergator/hpg.sh fetch         download the pooled .Rds into results/
#
# --fresh is the difference between finishing an interrupted run and starting a
# new one. Blocks are skipped by cell key, so after changing a statistic the old
# blocks would be reused and the run would look successful while reporting the
# old numbers. Changing what a replication computes means --fresh.
#
# One Duo prompt per invocation: the connection is multiplexed and held open for
# ten minutes, so the upload and the commands that follow share it.

set -euo pipefail

REMOTE=k.ekvall@hpg.rc.ufl.edu
BASE=/blue/k.ekvall/k.ekvall/Simulations/reconf
# This script lives in simulations/hipergator/; LOCAL is the repository root.
LOCAL="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

SSH_OPTS=(-o ControlMaster=auto -o ControlPath=/tmp/hpg-cm-%r@%h:%p
          -o ControlPersist=10m)
ssh_() { ssh "${SSH_OPTS[@]}" "$REMOTE" "$@"; }

push() {
    # The R code and the job scripts land flat in one directory on the cluster.
    rsync -a -e "ssh ${SSH_OPTS[*]}" "$LOCAL/simulations/code/" "$REMOTE:$BASE/"
    rsync -a -e "ssh ${SSH_OPTS[*]}" "$LOCAL"/simulations/hipergator/*.sh "$REMOTE:$BASE/"
    ssh_ "mkdir -p $BASE/Output $BASE/Results/cells $BASE/Results/width_blocks"
}

# submit <job script> <block dir>; a third argument --fresh discards the blocks
submit() {
    push
    if [ "${3:-}" = "--fresh" ]; then
        echo "==> discarding blocks in $2"
        ssh_ "rm -rf $BASE/$2 && mkdir -p $BASE/$2"
    fi
    ssh_ "sbatch $BASE/$1"
}

case "${1:-}" in
test)   push; ssh_ "sbatch $BASE/test_one.sh" ;;
run)    submit sims.sh  Results/cells        "${2:-}" ;;
width)  submit width.sh Results/width_blocks "${2:-}" ;;
status) ssh_ "squeue -u \$USER -o '%.12i %.10P %.12j %.8T %.10M %.5D %R'
              printf 'blocks on disk: %s coverage, %s width\n' \
                  \$(ls $BASE/Results/cells 2>/dev/null | wc -l) \
                  \$(ls $BASE/Results/width_blocks 2>/dev/null | wc -l)
              cat \$(ls -t $BASE/Output/merge*.out 2>/dev/null | head -1) 2>/dev/null \
                  | grep -E '^(all|INCOMPLETE|wrote)' || echo '(no pooling log yet)'" ;;
fetch)  mkdir -p "$LOCAL/results"
        rsync -a --itemize-changes -e "ssh ${SSH_OPTS[*]}" \
              "$REMOTE:$BASE/Results/*.Rds" "$LOCAL/results/" ;;
*)      sed -n '3,18p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
