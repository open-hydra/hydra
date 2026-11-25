# checklib.sh: sourced by a case's check.sh, called as `check.sh <run-dir> <exe>` (run.sh does it).
# Sets run, exe and cases (test/MI2/cases), and gives:
#   ok <label> / ko <label>          one PASS / FAIL line
#   expect <label> <command...>      PASS when the command succeeds
#   variant <name> <base-dir>        $run/<name>, a fresh copy of the base case's input.ini, INPUT/ and MESH/
#   runin <dir>                      runs the executable in <dir> (log <dir>/run.log); returns its exit code
#   ended <log>                      the run ended: GAME OVER, no [ERROR]
#   finish                           the summary line; exit 1 when a claim failed or none passed
run=${1:?usage: check.sh <run-dir> <exe>}; exe=${2:?usage: check.sh <run-dir> <exe>}
cases=$(cd "$(dirname "${BASH_SOURCE[0]}")/cases" && pwd)
npass=0; nfail=0
mkdir -p "$run"

ok() { npass=$((npass + 1)); echo "PASS  $*"; }
ko() { nfail=$((nfail + 1)); echo "FAIL  $*"; }
expect() { local label=$1; shift; if "$@"; then ok "$label"; else ko "$label"; fi; }

variant() {
  local d=$run/$1
  rm -rf "$d"; mkdir -p "$d"
  cp -rpL "$2/INPUT" "$2/input.ini" "$d/"
  if [ -d "$2/MESH" ]; then cp -rpL "$2/MESH" "$d/"; fi
  echo "$d"
}

runin() {
  (cd "$1" && rm -rf OUTPUT && mkdir OUTPUT && ls INPUT/*phase.txt > filelist.txt && ulimit -s unlimited \
     && KMP_STACKSIZE=100M timeout "${RUN_TIMEOUT:-3600}" "$exe" > run.log 2>&1)
}

ended() { grep -q 'GAME OVER' "$1" && ! grep -q '\[ERROR' "$1"; }

finish() {
  echo "== $(basename "$(dirname "$(readlink -f "${BASH_SOURCE[1]}")")"): $npass passed, $nfail failed"
  [ $nfail -eq 0 ] && [ $npass -gt 0 ]
  exit
}
