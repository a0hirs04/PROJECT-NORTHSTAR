#!/bin/bash
# Live monitor for the RC3 Arm C full-window rerun.
# Usage:
#   bash watch_rc3_armc_fullwindow.sh
#   bash watch_rc3_armc_fullwindow.sh --once

set -euo pipefail

export RC3_SWEEP_ROOT="/work/a0hirs04/PROJECT-NORTHSTAR/build/rc3_vismodegib_armc_fullwindow"
exec /bin/bash /home/a0hirs04/PROJECT-NORTHSTAR/watch_rc3_sweep.sh "$@"
