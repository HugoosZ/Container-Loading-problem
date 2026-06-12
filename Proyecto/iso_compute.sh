#!/bin/bash
# Comparacion iso-computo: baseline gamma=1 a 90s vs best-of-3 (gamma>1, 3x30s ya en CSV)
# Total equivalente: 90s = 3 semillas x 30s. Corre solo el baseline faltante. ~1h.
CSV=iso_compute.csv
for SET in 1 7 11 15; do
  F="BR${SET}.txt"; [ -f "$F" ] || { echo "FALTA $F"; continue; }
  MF=100; [ "$SET" -ge 8 ] && MF=98
  echo "== BR${SET} baseline a 90s =="
  ./bsg_rcl "$F" --problems 1-10 --time 90 --gamma 1.0 --seed 1 --min_fr $MF --csv $CSV 2>/dev/null
done
echo "Listo: $CSV (compara contra best-of-3 a 30s de resultados_finales.csv)"
