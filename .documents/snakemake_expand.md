# `expand()` in Snakemake: come funziona e come "looppa"

Riferimento pratico alla regola `featurecounts` in
`pipelines/rna-seq/Snakefile`, che usa `expand()` per costruire dinamicamente
la lista di BAM da passare a `featureCounts`.

## Cos'è un "sample" nella pipeline

`{sample}` non è una singola read o un batch di read consecutive: è un
**campione biologico/sperimentale indipendente**, con il proprio set completo
di FASTQ (R1 + R2). In `pipelines/rna-seq/config.yml`:

```yaml
samples:
  - SRR16039734
```

Se si aggiungessero altri sample (repliche biologiche, o condizioni diverse
come normal/tumor), ognuno verrebbe allineato separatamente da STAR (un BAM
per sample) e poi contato insieme dagli altri in un'unica matrice — è proprio
questo confronto **tra** sample che serve all'analisi differenziale a valle
(vedi `scripts/deseq2_analysis.R`).

## Cosa fa `expand()`

`expand(pattern, **wildcards)` non esegue nulla e non tocca il filesystem: è
puro **string templating**. Concettualmente equivale a:

```python
def expand(pattern, **wildcards):
    result = []
    for combinazione in itertools.product(*wildcards.values()):
        mapping = dict(zip(wildcards.keys(), combinazione))
        result.append(pattern.format(**mapping))
    return result
```

Fa il **prodotto cartesiano** di tutte le liste passate come keyword argument
e per ogni combinazione genera una stringa sostituendo i placeholder
`{nome}` con `str.format()`.

### Caso con una sola wildcard

```python
expand("results/rna-seq/{{run}}/aligned/{sample}_Aligned.sortedByCoord.out.bam",
       sample=SAMPLES)
```

Con `SAMPLES = ["s1", "s2"]`, il loop equivalente è:

```python
result = []
for sample in SAMPLES:
    result.append(f"results/rna-seq/{{run}}/aligned/{sample}_Aligned.sortedByCoord.out.bam")
```

→ lista piatta di 2 stringhe.

### Caso con più wildcard (prodotto cartesiano, loop annidato)

```python
expand("results/rna-seq/{run}/before/fastqc/{sample}_{read}_fastqc.html",
       run=RUN, sample=SAMPLES, read=[1, 2])
```

equivale a:

```python
result = []
for run in RUN:
    for sample in SAMPLES:
        for read in [1, 2]:
            result.append(f"results/rna-seq/{run}/before/fastqc/{sample}_{read}_fastqc.html")
```

Con 1 run, 2 sample, 2 read → 4 path generati.

## Il trucco della doppia graffa `{{run}}`

Dentro `expand()`, le graffe singole (`{sample}`, `{read}`, ...) vengono
sostituite subito con i valori passati come argomenti. Le graffe doppie
(`{{run}}`) sono l'escaping standard di `str.format()`: significano "graffa
letterale", non un placeholder da risolvere ora. Nel risultato restano come
`{run}` singolo, intatto.

Questo serve perché `run` non è un argomento di questa chiamata a `expand` —
è la **wildcard della regola**, che Snakemake risolve *dopo*, in base al
target effettivamente richiesto (es. da `rule all`, che passa
`run=RUN`). Se si scrivesse `{run}` con una sola graffa, `expand()`
cercherebbe un argomento `run=...` da lui stesso e, non trovandolo,
solleverebbe un errore.

**In una frase**: `expand` = genera N righe (una per sample) di un template;
`{{run}}` = "non toccare questa parte, la risolve Snakemake, non `expand`".

## Perché serve `expand` per `featurecounts` e non per le altre regole

Le regole come `samtools_index` girano **una volta per sample**: hanno
`{sample}` sia in input sia in output, e Snakemake risolve la wildcard da
solo per ogni istanza della regola.

`featurecounts` invece gira **una volta sola per run**, prendendo *tutti* i
BAM dei sample come input insieme:

```python
rule featurecounts:
    input:
        bams = expand("results/rna-seq/{{run}}/aligned/{sample}_Aligned.sortedByCoord.out.bam",
                       sample=SAMPLES),
        gtf = GENOME_GTF
    output:
        counts = "results/rna-seq/{run}/counts/gene_counts.txt"
    ...
    shell:
        "featureCounts -a {input.gtf} -o {output.counts} -T {threads} -p {input.bams}"
```

Il suo output (`gene_counts.txt`) non contiene `{sample}`, quindi Snakemake
non ha modo di dedurre automaticamente la lista di BAM: va costruita a mano,
ed è esattamente il compito di `expand`.

## Dalla lista Python al comando shell

1. **In Python**, `input.bams` è una lista di stringhe:
   ```python
   ["results/rna-seq/run1/aligned/s1_Aligned.sortedByCoord.out.bam",
    "results/rna-seq/run1/aligned/s2_Aligned.sortedByCoord.out.bam"]
   ```
2. **Snakemake appiattisce la lista** unendo gli elementi con uno spazio
   quando la interpola dentro `shell: "... {input.bams}"` — è il
   comportamento di default per le liste dentro `{...}`.
3. **A quel punto è solo testo bash**: la shell vede N argomenti posizionali
   separati da spazio, come se fossero stati scritti a mano.
4. **`featureCounts`** (tool esterno, non sa nulla di Snakemake) riceve
   quegli argomenti come lista di BAM da processare *insieme*, e produce
   un'unica tabella con una colonna di conteggi per ciascun BAM in input —
   la matrice `gene × sample` che poi legge `scripts/deseq2_analysis.R`.

## Riepilogo dei livelli di aggregazione

| Livello | Chi lo fa | Cosa succede |
|---|---|---|
| Generazione lista path | `expand()` | Prodotto cartesiano + string templating, nessuna aggregazione dati |
| Gate di dipendenza | Snakemake | La regola aspetta che *tutti* i BAM in `input.bams` esistano prima di girare |
| Costruzione comando | Snakemake | La lista Python diventa argomenti spazio-separati nella riga di shell |
| Aggregazione dati reale | `featureCounts` | Prende N BAM come argomenti multipli e produce un'unica tabella `gene × sample` |

## Dove vederlo in pratica

- Regola Snakemake: `featurecounts` in `pipelines/rna-seq/Snakefile`
- Output: `results/rna-seq/{run}/counts/gene_counts.txt`
- Consumato da: `scripts/deseq2_analysis.R`
