$CSV = "iso_compute.csv"
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
    Write-Host "== BR${SET} baseline a 90s =="
    
    ./bsg_rcl.exe $F --problems 1-10 --time 90 --gamma 1.0 --seed 1 --min_fr $MF --csv $CSV 2>$null
}

Write-Host "Listo: $CSV (compara contra best-of-3 a 30s de resultados_finales.csv)"
