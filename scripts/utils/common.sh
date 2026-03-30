#!/bin/bash

# Color definitions
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;38;5;186m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
  echo -e "${CYAN}[INFO]${NC} $1"
}

log_success() {
  echo -e "${GREEN}[OK]${NC} $1"
}

log_warn() {
  echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
  echo -e "${RED}[ERROR]${NC} $1" >&2
}

# Prompt function
ask_user() {
  echo -ne "${BLUE}$1${NC} [y/N]: "
  read -r reply
  [[ $reply =~ ^[Yy]$ ]]
}

# Assert a required command exists
require_command() {
  command -v "$1" >/dev/null 2>&1 || {
    log_error "Required command '$1' not found. Please install it."
    exit 1
  }
}
