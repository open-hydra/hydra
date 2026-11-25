# figaro_autocomplete.sh -- Bash completion for figaro

_figaro_completions()
{
  local cur prev opts subsolvers solvers hydras atlas_opt
  COMPREPLY=()
  cur="${COMP_WORDS[COMP_CWORD]}"
  prev="${COMP_WORDS[COMP_CWORD-1]}"

  opts="doctor deps clean compile build kill ATLAS MOSE FUSS GROOT ARES hydra-MF hydra-AF hydra-MI2"
  subsolvers="MOSE FUSS GROOT ARES"
  # `build` also takes the hydra executables themselves (figaro build hydra-MI2 ...)
  hydras="hydra-MF hydra-AF hydra-MI2"
  progs="ATLAS GRIB MOSE FUSS GROOT ARES hydra-MF hydra-AF hydra-MI2"
  atlas_opt="GPB ICB BCB import"
  deps_opts="--update ATLAS MOSE FUSS GROOT ARES ICE IGLOO ORION OSLO FLINT"

  case ${prev} in
    doctor)
      COMPREPLY=()
      return 0
      ;;
    deps)
      COMPREPLY=( $(compgen -W "${deps_opts}" -- ${cur}) )
      return 0
      ;;
    compile)
      COMPREPLY=( $(compgen -W "${progs}" -- ${cur}) )
      return 0
      ;;
    build)
      COMPREPLY=( $(compgen -W "${subsolvers} ${hydras}" -- ${cur}) )
      return 0
      ;;
    ATLAS)
      COMPREPLY=( $(compgen -W "${atlas_opt}" -- ${cur}) )
      return 0
      ;;
    hydra*|MOSE|FUSS|GROOT|ARES)
      COMPREPLY=( $(compgen -W "solve" -- ${cur}) )
      return 0
      ;;
    *)
      COMPREPLY=( $(compgen -W "${opts}" -- ${cur}) )
      return 0
      ;;
  esac
}

complete -F _figaro_completions figaro
