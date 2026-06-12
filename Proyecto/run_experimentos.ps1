$T = 30
$CSV = "resultados_finales.csv"
$sets = @(1, 7, 11, 15)

foreach ($SET in $sets) {
    $F = "BR${SET}.txt"
    if (-not (Test-Path $F)) {
        Write-Host "FALTA $F"
        continue
    }
    
    $MF = 100
    if ($SET -ge 8) { 
        $MF = 98 
    }
    Write-Host "== BR${SET} (min_fr=$MF) =="
    
    # Baseline (gamma = 1.0)
    ./bsg_rcl.exe $F --problems 1-10 --time $T --gamma 1.0 --seed 1 --min_fr $MF --csv $CSV 2>$null
    
    # Experimentación con gamma 1.5 y 2.0 y diferentes semillas
    $gammas = @("1.5", "2.0")
    $seeds = @(1, 2, 3)
    
    foreach ($G in $gammas) {
        foreach ($S in $seeds) {
            ./bsg_rcl.exe $F --problems 1-10 --time $T --gamma $G --seed $S --min_fr $MF --csv $CSV 2>$null
        }
    }
}

Write-Host "Listo. Resultados guardados en $CSV"
