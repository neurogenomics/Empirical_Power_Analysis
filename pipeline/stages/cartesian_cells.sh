# shellcheck shell=bash
# Stage cartesian_cells: see _cartesian_post.sh.
source "$(dirname "${BASH_SOURCE[0]}")/_cartesian_post.sh"
stage_main() { cartesian_post cells; }
