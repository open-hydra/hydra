#!/bin/bash

# utils/print_usage.sh - Displays usage and help for figaro

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
source "$SCRIPT_DIR/utils/common.sh"

print_usage() {
  echo
  echo "Pronto prontissimo son come il fulmine"
  echo "Sono il factotum della citta'"
  echo

  echo -e "${YELLOW}Usage:${NC} figaro <command> [options] [args]"
  echo

  echo -e "${YELLOW}Available commands:${NC}"
  echo

  echo -e "${RED}git submodules tasks${NC}"
  echo -e "  ${CYAN}submodule update [options] <mod>${NC}       Update git submodule. If <mod> not specified, update all submodules"
  echo -e "     ${CYAN}-r, --remote${NC}                            Update to latest remote version"
  echo -e "  ${CYAN}submodule delete [options] <mod>${NC}       Delete git submodule. If <mod> not specified, delete all submodules"
  echo -e "     ${CYAN}-f, --force${NC}                             Delete submodule even if it has uncommitted changes"
  echo

  echo -e "${RED}Build tasks${NC}"
  echo -e "  ${CYAN}build <solver>${NC}                         Build solver (both submodules and hydras)"
  echo -e "     ${CYAN}options${NC}                                 Specific for each solver, see the proper install.sh"
  echo

  echo -e "${RED}Compile tasks${NC}"
  echo -e "  ${CYAN}compile [options]  <solver>${NC}            Full compilation of solver (both submodules and hydras)"
  echo -e "     ${CYAN}-e, --edit${NC}                              Edit the CMakePreset file before compiling"
  echo

  echo -e "${RED}Simulation tasks${NC}"
  echo -e "  ${CYAN}clean${NC}                                  Clean simulation outputs"
  echo -e "  ${CYAN}ATLAS <exe-list>${NC}                       Run list of executables and store results in dummy folder"
  echo -e "  ${CYAN}<solver> [options] solve${NC}               Launch solver"
  echo -e "     ${CYAN}-b, --background${NC}                        Run in background"
  echo -e "     ${CYAN}-p, --parallel <n>${NC}                      Run with <n> threads"
  echo -e "  ${CYAN}kill${NC}                                   Kill currently running simulation"
  echo

  echo -e "${YELLOW}Examples:${NC}"
  echo

  echo "figaro build hydra-MF --compilers=gnu --use-openmp"
  echo "      Build hydra solver 'MF' with GNU compilers and OpenMP enabled"
  echo

  echo "figaro compile hydra-MF"
  echo "      Compile an existing build of hydra-MF using CMake presets"
  echo

  echo "figaro build hydra-MI2 --compilers=intel --use-openmp"
  echo "      Build hydra solver 'MI2' (MOSE + ICE + IGLOO). OpenMP only -- MI2 refuses"
  echo "      to configure with --use-mpi, because IGLOO needs the full gas field on"
  echo "      every rank while MOSE partitions over blocks"
  echo

  echo "figaro build MOSE --use-openmp"
  echo "      Build MOSE solver with OpenMP enabled"
  echo

  echo "figaro submodule update MOSE"
  echo "      Update MOSE git submodule to the commit specified in the main repository"
  echo

  echo "figaro MOSE -p 4 solve"
  echo "      Launch MOSE solver using 4 threads"
  echo

  echo -e "${YELLOW}Notes:${NC}"
  echo "  - Each hydra solver is built in its own directory: build/<suffix>"
  echo "  - CMakePresets.json is generated automatically during build"
  echo "  - hydra executables: MF (MOSE+FUSS), AF (ARES+FUSS), MI2 (MOSE+ICE+IGLOO)"
  
  exit 0
}
