# beckstead.sh: sourced by twins/check.sh and burnout/check.sh.
# beckstead <case dir> <end key line> <end token>: the burning keys of [GPB-Phase2] of the case's input.ini, its phase
# line and its models row set to combustion Beckstead (the twins' constants), ending with the given key (d-stop or
# d-burnout) and token.
BECK='combustion = Beckstead\nK-burn     = 1.39e-6\nn-burn     = 1.8\nX-eff      = 0.5\nT-ign      = 290\nbeta-part  = 0.3\nq-comb     = 1.0e6'
ROW='A combustion Beckstead K-burn=1.39e-06 n-burn=1.8 X-eff=0.5 beta-part=0.3 T-ign=290 q-comb=1000000'
beckstead() {
  awk -v blk="$BECK\n$2" '/^\[/ {s = $0; print; if (s == "[GPB-Phase2]") printf "%s\n", blk; next}
       s == "[GPB-Phase2]" && $1 ~ /^(combustion|T-burn|K-burn|n-burn|X-eff|T-ign|beta-part|q-comb|d-stop|d-burnout)$/ {next}
       {print}' "$1/input.ini" > "$1/ini" && mv "$1/ini" "$1/input.ini"
  sed -i 's/^A 1 combustion=.*/A 1 combustion=Beckstead/' "$1/INPUT/part-phase.txt"
  sed -i "s/^A combustion .*/$ROW $3/" "$1/INPUT/part-models.txt"
}
