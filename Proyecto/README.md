# BSG-RCL: Beam Search Greedy con Lista Restringida de Candidatos para el SCLP

**Proyecto CIT3352 — Algoritmos Exactos y Metaheurísticas · Entrega 2**

Resumen del diseño, la implementación y el plan experimental de la solución propuesta
para el *Single Container Loading Problem* (SCLP).

---

## 1. Problema

El SCLP consiste en cargar un contenedor de dimensiones L×W×H con un conjunto de cajas
rectangulares, maximizando el porcentaje del volumen del contenedor ocupado. Las cajas
no pueden solaparse, deben quedar completamente dentro del contenedor con caras
paralelas a las de éste, y cada caja tiene restricciones de orientación (qué dimensiones
pueden quedar en el eje vertical). El problema es NP-hard.

## 2. Algoritmo base y declaración de fuente

- **Algoritmo base:** BSG-CLP, de I. Araya y M.-C. Riff, *"A beam search approach to the container loading problem"*, Computers & Operations Research 43 (2014) 100–107.

- **Origen del código:** los autores del paper publican únicamente un ejecutable binario (no código fuente), por lo que **la implementación completa de este proyecto es propia**, escrita en C++ desde el pseudocódigo del paper (Algoritmos 1, 2 y 3) y la descripción de sus seis elementos clave (K1–K6). El diseño del algoritmo base se declara como proveniente del paper; el código y el operador nuevo son contribución del equipo.

Componentes implementados fielmente al paper:

| Elemento | Descripción | Sección del paper |
|---|---|---|
| K1 | Espacio residual con *cover representation* (cuboides solapados + filtro de no maximales) | 3.1.2 |
| K2 | Bloques simples y generales (`GeneralBlockGeneration`, parámetros `max_bl`, `min_fr`) | 3.2.1 |
| K3 | Selección de espacio libre por distancia Manhattan mínima (empate por mayor volumen) | 3.2.5 |
| K4 | Función de evaluación f(b,r) = V(b) − Vloss(b,r), con knapsack (KPA) por programación dinámica, cacheado una vez por estado | 3.2.6 |
| K5 | Colocación del bloque en la esquina ancla del cuboide | 3.2.5 |
| K6 | Beam search + remoción de estados similares + mecanismo de doble esfuerzo (w ← ⌈√2·w⌉) | 3.2.2–3.2.4 |

Simplificación declarada: la combinación por pares de bloques generales limita el pool a los 1.200 bloques de mayor volumen y a 3 rondas (en el original, de Fanslau y Bortfeldt 2010, la combinación itera hasta `max_bl` sin esa cota de pool).

## 3. Contribución propia: el operador Expand-RCL

En BSG-CLP, el paso `expand` selecciona **determinísticamente** los `w` mejores bloques según f(b,r). Nuestra propuesta lo generaliza a una **construcción estocástica sesgada tipo GRASP**:

1. Para el espacio libre seleccionado, se construye una **lista restringida de candidatos (RCL)** con los `K = ⌈γ·w⌉` mejores bloques factibles según f(b,r).
2. Se **muestrean `w` bloques sin reemplazo** de la RCL, con **pesos lineales por rango** (el mejor bloque tiene el mayor peso): aleatorio pero sesgado hacia la calidad, en el espíritu de los *biased random keys* del BRKGA que motivó la propuesta inicial del equipo.
3. **Con γ = 1 el operador se reduce exactamente al BSG-CLP original**, lo que da el baseline determinista en el mismo código y garantiza comparabilidad perfecta.

Justificación bibliográfica: Araya y Riff señalan explícitamente en sus conclusiones que *agregar aleatorización a la selección de bloques/espacios* es una línea de trabajo futuro prometedora para aumentar la diversidad entre estados del beam. Este operador
realiza esa línea.

Decisiones de diseño del operador:

- **RCL por cardinalidad (top-K) y no por umbral de calidad α:** invariante a la escala de f(b,r), que varía órdenes de magnitud entre etapas de la construcción.
- **Pesos lineales por rango:** robustos frente a diferencias extremas de volumen entre bloques (los pesos proporcionales a f degeneran en cuasi-determinismo).
- **La aleatorización NO entra al greedy de evaluación:** el greedy se mantiene determinista para que la evaluación de estados del beam no tenga varianza.

## 4. Hipótesis

Bajo el mismo presupuesto de tiempo, γ > 1 mejora la utilización promedio de volumen respecto del baseline determinista (γ = 1), con mayor efecto en instancias **fuertemente heterogéneas** (BR8–BR15), donde la diversidad de bloques hace más probable que el camino greedy puro descarte alternativas buenas.

## 5. Parámetros

| Parámetro | Valores | Notas |
|---|---|---|
| γ (factor RCL) | {1, 2, 3} | γ=1 es el baseline; γ=2 apuesta principal |
| Semillas | {1, 2, 3} para γ>1 | γ=1 es determinista (1 corrida) |
| Tiempo límite | 30 s/instancia | comparable a la columna "30 s" de la Tabla 2 del paper |
| `max_bl` | 10.000 | igual al paper |
| `min_fr` | 100 (BR0–BR7) / 98 (BR8–BR15) | igual al paper |

## 6. Diseño experimental

- **Instancias:** BR1, BR7 (débilmente heterogéneas), BR11, BR15 (fuertemente heterogéneas); primeras 10 instancias de cada set → 40 instancias.
- **Presupuesto total:** 40 × 30 s × 7 corridas (1 de γ=1 + 3 semillas × 2 valores de γ) ≈ **2,3 horas de cómputo**.
- **Métricas:** % de utilización por instancia; promedio por set; comparación γ=1 vs γ=2 con test de Wilcoxon de rangos con signo (mismo test del paper), usando la media de las 3 semillas por instancia. Reportar además mejor/peor semilla (varianza).

- **Obtención de instancias:** los sets BR1–BR7 corresponden a los archivos `thpack1`–`thpack7` de la OR-Library de Beasley (Bischoff y Ratcliff); BR8–BR15(`thpack8`–`thpack15`) fueron generados por Davies y Bischoff. **Verificar el formato de archivo descargado contra el parser** (formato esperado: nº de problemas; por problema: id y semilla; L W H; nº de tipos; por tipo: índice, d1, flag1, d2, flag2, d3, flag3, cantidad — el flag indica si esa dimensión puede ir vertical).

## 7. Estructura del código (`bsg_rcl.cpp`, autocontenido)

| Sección | Contenido |
|---|---|
| 1 | Parser de instancias formato BR + orientaciones factibles |
| 2 | Geometría: cuboides, solapamiento, contención |
| 3 | Generación de bloques simples y generales |
| 4 | Estado, selección de espacios, colocación + actualización de la cover representation |
| 5 | KPA (subset-sum acotado con bitsets) y f(b,r) |
| 6 | Greedy de evaluación (con poda f(b,r) ≤ V(b) sobre bloques ordenados por volumen) |
| 7 | **Expand-RCL (operador propuesto)** |
| 8 | Beam search con remoción de estados similares |
| 9 | Bucle de doble esfuerzo, salida temprana al óptimo, runner y CSV |
| 10 | Selftest (geometría + end-to-end + determinismo de γ=1) |

## 8. Compilación y uso

```bash
make                       # o: g++ -O2 -std=c++17 -o bsg_rcl bsg_rcl.cpp
./bsg_rcl --selftest       # pruebas unitarias
./bsg_rcl thpack1.txt --problems 1-10 --time 30 --gamma 1.0 --csv res.csv   # baseline
./bsg_rcl thpack1.txt --problems 1-10 --time 30 --gamma 2.0 --seed 1 --csv res.csv
./bsg_rcl thpack11.txt --problems 1-10 --time 30 --gamma 2.0 --seed 1 --min_fr 98 --csv res.csv
bash run_experimentos.sh   # barrido completo del diseño experimental
```

<!-- ## 9. Validación realizada

- Selftest de geometría: residuales correctos al colocar un bloque en esquina
  (3 cuboides maximales esperados). ✔
- End-to-end óptimo: 8 cubos de 5³ en contenedor 10³ → 100% de relleno. ✔
- Determinismo de γ=1: misma solución con semillas distintas. ✔
- Instancia sintética competitiva (volumen de cajas > contenedor): baseline 98,19%;
  con γ=3 una semilla alcanzó **98,37%**, superando al determinista — evidencia
  preliminar a favor de la hipótesis. ✔
- Bloques generales (`min_fr` 98): generación y carga de bloques multi-tipo. ✔ -->

## 9. Pendientes antes de la entrega

1. Descargar los sets BR y validar el parser contra el formato real de los archivos.
2. Hito de sanidad: correr γ=1 en BR1 (10 instancias, 30 s) y verificar relleno en el
   rango 90–95% (la Tabla 2 del paper reporta ~95,4% para las 100 instancias en C++
   optimizado; valores muy por debajo indicarían un bug).
3. Ejecutar `run_experimentos.sh` (≈2,3 h) y tabular el CSV.
4. Redactar el informe: descripción del algoritmo (secciones 2–3 de este documento),
   experimentos y análisis (sección 6 + resultados), conclusiones y trabajo futuro
   (p. ej., extender la aleatorización a la selección de espacios libres con los q
   cuboides de menor distancia Manhattan; hacer adaptativos `min_fr` y γ).
