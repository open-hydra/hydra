# figaro_autocomplete.sh -- Bash completion for figaro

_figaro_completions()
{
  local cur prev opts subsolvers solvers atlas_opt
  COMPREPLY=()
  cur="${COMP_WORDS[COMP_CWORD]}"
  prev="${COMP_WORDS[COMP_CWORD-1]}"

  opts="doctor submodule clean compile build kill ATLAS MOSE IGLOO ICE FUSS MOSKA Q2D hydra-MM hydra-MI2 hydra-MF"
  subsolvers="MOSE IGLOO ICE FUSS MOSKA Q2D"
  progs="ATLAS GRIB MOSE IGLOO ICE FUSS MOSKA Q2D hydra-MM hydra-MI2 hydra-MF hydra-MQ2D"
  atlas_opt="GPB ICB BCB import"
  git_opts="update delete"

  case ${prev} in
    doctor)
      COMPREPLY=()
      return 0
      ;;
    submodule)
      COMPREPLY=( $(compgen -W "${git_opts}" -- ${cur}) )
      return 0
      ;;
    compile)
      COMPREPLY=( $(compgen -W "${progs}" -- ${cur}) )
      return 0
      ;;
    build)
      COMPREPLY=( $(compgen -W "${subsolvers}" -- ${cur}) )
      return 0
      ;;
    ATLAS)
      COMPREPLY=( $(compgen -W "${atlas_opt}" -- ${cur}) )
      return 0
      ;;
    MOSE|IGLOO|ICE|FUSS|MOSKA|Q2D|hydra*)
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
