# figaro_autocomplete.sh -- Bash completion for figaro

_figaro_completions()
{
  local cur prev opts subsolvers solvers atlas_opt
  COMPREPLY=()
  cur="${COMP_WORDS[COMP_CWORD]}"
  prev="${COMP_WORDS[COMP_CWORD-1]}"

  opts="doctor submodule clean compile build kill ATLAS MOSE FUSS GROOT hydra-MF"
  subsolvers="MOSE FUSS GROOT"
  progs="ATLAS GRIB MOSE FUSS GROOT hydra-MF"
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
    hydra*|MOSE|FUSS|GROOT)
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
