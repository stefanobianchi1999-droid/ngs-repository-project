# Lettura dei grafici di `plot-bamstats`

Questi grafici vengono generati dalla regola `plot_bamstats` della pipeline
RNA-seq (`pipelines/rna-seq/Snakefile`), a partire dall'output di
`samtools stats` sul BAM allineato con STAR. Servono come controllo qualità
(QC) post-allineamento.

## Grafici

- **Insert size** — distribuzione della **lunghezza del frammento originale**
  di DNA/cDNA da cui provengono le read. In sequenziamento paired-end si
  legge solo ciascuna estremità del frammento (R1 da un lato, R2 dall'altro);
  l'insert size non è quindi la distanza tra i punti di partenza delle due
  read, ma la lunghezza totale del frammento, calcolata dopo l'allineamento
  come distanza tra la coordinata più esterna di R1 e quella più esterna di
  R2 sul riferimento — include anche la porzione centrale del frammento non
  coperta direttamente dalle read. Un picco stretto e ben definito indica una
  libreria di sequenziamento consistente (frammentazione uniforme); code
  larghe o distribuzioni bimodali possono segnalare problemi nella
  preparazione della libreria (frammentazione irregolare, contaminazione,
  ecc.).

- **GC content** — distribuzione del contenuto GC (%) nelle read,
  confrontata con la distribuzione attesa. Il grafico mostra curve separate
  per **"first fragment"** (R1) e **"last fragment"** (R2) — stessa
  terminologia usata nei grafici di qualità per ciclo (vedi sopra) — così da
  poter individuare bias di GC specifici di un mate rispetto all'altro.
  Deviazioni marcate dalla curva teorica possono indicare contaminazione da
  altre specie o bias di amplificazione PCR.

- **Per-base sequence content** (`acgt-cycles.png`) — frequenza di ciascuna
  base (A/C/G/T) per posizione lungo la read (ciclo di sequenziamento). Un
  po' di bias nelle prime posizioni è normale (es. da primer random), ma le
  proporzioni dovrebbero stabilizzarsi rapidamente.

- **Quality per cycle** (`quals.png`, `quals2.png`, `quals3.png`,
  `quals-hm.png`) — quattro viste della qualità Phred per ciclo di
  sequenziamento: andamento medio/mediano per posizione, boxplot per
  posizione, e una heatmap di densità dei valori di qualità. Un calo
  progressivo verso la fine della read è fisiologico (tipico Illumina); un
  crollo brusco e precoce indica un problema di sequenziamento.

  **Cosa significa "cycle" (asse X)**: nel sequenziamento Illumina, un
  "ciclo" è il singolo passaggio in cui lo strumento legge una base di ogni
  read — al ciclo 1 legge la 1ª base di tutte le read, al ciclo 2 la 2ª, e
  così via fino alla lunghezza massima della read (150 nel nostro caso).
  Quindi l'asse X di questo grafico (e allo stesso modo di "Per-base
  sequence content" sopra) non è una coordinata sul genoma, ma **la
  posizione della base all'interno della read** (1ª, 2ª, ... 150ª base),
  aggregata su tutte le read del campione.

  `quals2.png` e `quals3.png` in particolare separano la qualità per
  **"first fragment"** e **"last fragment"**: nella terminologia di
  `samtools stats`, il "first fragment" è la prima read della coppia
  paired-end (R1, il mate letto per primo), mentre il "last fragment" è la
  seconda read (R2, l'altra estremità del frammento, letta per ultima). La
  distinzione è utile perché sui sequenziatori Illumina R2 tende
  tipicamente ad avere una qualità leggermente inferiore a R1, soprattutto
  verso la fine della read.

  **Cos'è "quality"**: è il **Phred quality score**, un valore che esprime
  quanto il basecaller è sicuro di aver identificato correttamente quella
  base in quel ciclo, calcolato dall'intensità/purezza del segnale
  fluorescente catturato in quel ciclo. Scala logaritmica:
  `Q = -10 × log10(P)`, dove P è la probabilità stimata di errore (Q20 = 1
  errore ogni 100 basi = 99%; Q30 = 1 ogni 1.000 = 99.9%).

  **`quals-hm.png` (heatmap)**: stessa informazione delle altre viste, ma
  come densità — asse X il ciclo, asse Y il valore di qualità Phred, colore
  di ogni cella quante read hanno esattamente quel valore in quel ciclo. A
  differenza di una linea di media/mediana, mostra la forma completa della
  distribuzione per ciclo: permette di distinguere un calo di qualità
  uniforme su tutte le read (banda stretta che scende gradualmente) da un
  sottogruppo problematico di read a bassa qualità (una seconda banda o
  coda separata, che una semplice media nasconderebbe).

- **Indels per cycle** (`indel-cycles.png`) — frequenza di
  inserzioni/delezioni introdotte durante l'allineamento, per posizione
  nella read. Un aumento agli estremi della read è un normale effetto
  "bordo" dell'allineatore; picchi anomali al centro sono sospetti.

  **Come leggere un picco anomalo concentrato in poche posizioni** (es. un
  picco netto tra i cicli 60-80 molto più alto del rumore di fondo): il
  grafico aggrega gli indel per posizione **nella read**, sommando su
  tutte le read del campione. Perché molte read diverse mostrino un indel
  proprio nella stessa posizione serve una causa sistematica, non rumore
  casuale (che sarebbe distribuito uniformemente). Cause tipiche:
  - un locus del genoma molto coperto con problemi di allineamento locale
    (es. regioni ripetute/multi-copia, in *E. coli* tipicamente operoni
    rRNA/tRNA);
  - un piccolo errore/gap nell'assembly del genoma di riferimento in quel
    punto, che l'allineatore "corregge" con un indel sempre nello stesso
    posto;
  - contaminazione da adapter residuo non rimosso dal trimming, che
    l'allineatore forza come indel per poter comunque allineare la read.

  Da valutare sempre in proporzione al totale delle read (un picco che
  sembra enorme sul grafico può corrispondere a una frazione minima, es.
  <0.1%, coerente con un error rate complessivo comunque basso in
  `samtools stats`). Per indagare: filtrare le read con indel in quella
  posizione (via CIGAR) e verificare su quale coordinata del genoma
  cadono — se convergono su un unico locus, conferma un problema
  localizzato piuttosto che un problema diffuso di qualità dei dati.

- **Indel lengths** (`indel-dist.png`) — distribuzione delle lunghezze degli
  indel osservati. Nella maggior parte dei dataset prevalgono indel di 1-2
  basi.

  **Come leggerlo**: asse X = lunghezza dell'indel in basi (1bp, 2bp,
  3bp, ...), asse Y = quante volte è stato osservato un indel di quella
  lunghezza esatta, in genere con due curve separate per inserzioni e
  delezioni. Nel formato CIGAR di un allineamento BAM corrispondono agli
  operatori `I` (insertion) e `D` (deletion), con il numero che li precede
  che è proprio questa lunghezza. La curva dovrebbe scendere rapidamente
  (quasi esponenziale) man mano che la lunghezza aumenta: gli indel corti
  (1-2bp) sono i più comuni sia come piccoli errori tecnici di
  allineamento sia come eventi biologici reali; una **coda lunga anomala**
  (tanti indel da 10, 20, 50bp) è sospetta e di solito indica soft-clipping
  mal gestito, contaminazione da adapter non trimmata, o discrepanze
  strutturali nel genoma di riferimento.

- **Mapped depth vs GC** (`gc-depth.png`) — relazione tra profondità di
  copertura e contenuto GC locale.

  **Come viene costruito**: il genoma di riferimento viene scandito a
  finestre di lunghezza fissa. Per ciascuna finestra si calcolano due
  valori indipendenti: il **contenuto GC locale** (% di basi G+C nella
  sequenza di riferimento, dipende solo dal genoma) e la **profondità di
  copertura locale** (quante read mediamente coprono quella finestra,
  dipende dai dati allineati). L'asse X di `plot-bamstats` **non è il GC%
  grezzo**, ma il **percentile delle finestre ordinate per GC** (0-100);
  l'etichetta GC% mostrata in alto è solo un riferimento indicativo a
  quale GC% approssimativo corrisponde quel percentile. Per ciascun bin
  percentile si riporta mediana, 25°-75° e 10°-90° percentile della
  profondità osservata (asse Y): in pratica, "tra le finestre che stanno a
  questo percentile di GC, quanto sono coperte in media dai miei dati?".

  **Perché può risultare molto rumoroso/frastagliato invece che una
  campana smussata**: dipende dal genoma di riferimento. Genomi piccoli e
  con GC molto omogeneo (es. *E. coli*, che si aggira quasi tutto intorno
  al 50-51% GC, senza le forti variazioni di un genoma eucariotico) hanno
  pochissime finestre con GC estremo. Quando un bin percentile è
  supportato da poche finestre, la stima di profondità diventa
  statisticamente instabile → picchi o crolli bruschi ai margini della
  distribuzione, invece di una curva liscia. Bande 10°-90° percentile
  molto larghe indicano inoltre alta variabilità di copertura di per sé
  (in RNA-seq riflette il livello di espressione genica: geni molto
  espressi = copertura alta, geni poco espressi = copertura bassa), non
  necessariamente un problema tecnico. Il segnale più affidabile da
  guardare è se la **mediana** resta ragionevolmente stabile nella fascia
  centrale della distribuzione (es. 30°-80° percentile): oscillazioni forti
  solo ai bordi, con pochi dati a supporto, sono attese e non
  necessariamente un campanello d'allarme.

  **Perché esiste questo bias**: nasce nella library prep, non
  nell'allineamento — durante l'amplificazione PCR. Regioni molto ricche
  in GC (>70%) formano strutture secondarie stabili, più difficili da
  denaturare/amplificare → tendono a essere sotto-rappresentate. Regioni
  molto povere in GC (<30%, ricche in A/T) hanno temperature di melting
  basse e amplificano in modo instabile → anche loro spesso
  sotto-rappresentate. Il GC intermedio (40-60%) amplifica in modo più
  efficiente → relativamente sovra-rappresentato. Il risultato atteso è
  quindi una curva a campana (non una linea piatta): bassa profondità agli
  estremi di GC, picco nella zona centrale.

  **Perché conta per l'RNA-seq**: se lo scopo finale è la quantificazione
  dell'espressione (conteggio read per gene, es. con `featureCounts`), un
  gene con GC estremo può risultare sotto- o sovra-stimato per puro bias
  tecnico di libreria, non per differenza biologica reale — uno dei motivi
  per cui strumenti come DESeq2 applicano normalizzazione, e in analisi più
  rigorose si corregge esplicitamente il GC bias prima della
  quantificazione.

  **Cosa guardare**: curva a campana morbida centrata attorno al GC medio
  del genoma di riferimento (~50% per *E. coli*) è normale; una curva
  molto stretta, crolli bruschi a GC intermedio, o un'asimmetria marcata
  (crollo solo da un lato) indicano possibili problemi di libreria o bias
  specifico del kit di prep usato.

- **Read lengths** (`read-length.png`) — distribuzione delle lunghezze delle
  read dopo il trimming con fastp.

  **Cosa significa l'asse Y "Unbinned count / Bin width"**: il valore
  mostrato non è il conteggio grezzo di read per bin, ma il conteggio
  diviso per la larghezza del bin — una densità (read per bp), non un
  istogramma a conteggio puro. Serve perché `samtools stats` può usare bin
  di larghezza non uniforme (più larghi dove i dati sono sparsi, più
  stretti dove sono densi, per contenere la dimensione del file di
  output): un conteggio grezzo confonderebbe un bin largo con tante read
  "diluite" con un bin stretto altrettanto pieno. Dividendo per la
  larghezza, il valore diventa confrontabile tra bin diversi — "unbinned"
  indica che il dato è stato riportato a un'unità standard (per bp)
  invece di restare legato all'ampiezza arbitraria del bin scelto.

  **Cosa aspettarsi**: con read trimmate a lunghezza pressoché fissa (fastp)
  e senza soft-clipping eccessivo in allineamento, ci si aspetta un picco
  singolo, stretto e alto intorno alla lunghezza nominale delle read (150bp
  nel nostro caso, coerente con l'`avg read length: 150` della tabella
  numerica). Più picchi o una coda ampia verso lunghezze minori
  indicherebbero un trimming disomogeneo su un sottoinsieme di read, o
  soft-clipping esteso durante l'allineamento.

  **Perché il grafico ha due assi Y**: la legenda mostra due curve, con
  scale indipendenti.
  - **Asse sinistro ("Read Count", arancione = "unbinned")** — conteggio
    grezzo, non raggruppato, delle read che hanno esattamente quella
    lunghezza intera in bp. Scala molto alta (milioni) perché con dati
    trimmati a lunghezza fissa quasi tutte le read cadono su un unico
    valore (es. 150bp), producendo un picco altissimo e strettissimo.
  - **Asse destro ("count/bin width", blu = "normalized bins")** — la
    stessa distribuzione ma raggruppata in bin più larghi e normalizzata
    per larghezza del bin (vedi sopra), su una scala molto più piccola e
    indipendente.

  Se condividessero lo stesso asse, il picco enorme e stretto sulla
  lunghezza nominale schiaccerebbe a zero qualunque variazione nella coda
  della distribuzione. Con due scale separate il grafico mostra
  contemporaneamente quanto è dominante il picco a lunghezza fissa (asse
  sinistro) e come si distribuisce la piccola minoranza di read a
  lunghezza diversa, es. per trimming (asse destro, "zoomato" per essere
  visibile).

## Tabella numerica (Reads / Bases)

Riepilogo aggregato mostrato sotto ai grafici:

- **total / mapped / zero MQ** — read totali, quante mappate (%), quante con
  MAPQ zero (allineamenti ambigui/multi-mapping).
- **duplicated** — read duplicate marcate (0 se non è stato eseguito un
  passaggio di duplicate marking a monte).
- **error rate** — frazione di basi mismatch rispetto al riferimento negli
  allineamenti.

### Esempio: run `run_2026-07-15`, campione `SRR16039734`

| Metrica | Valore |
|---|---|
| Read totali | 7.807.630 |
| Mappate | 100.0% |
| Duplicate | 0 (0.0%) |
| MAPQ zero | 18.378 (0.2%) |
| Lunghezza media read | 150 bp |
| Basi totali | 1.170.055.200 |
| Basi mappate | 99.9% |
| Error rate | 0.00% |

Numeri di questo tipo (mapped ~100%, error rate ~0%) indicano un allineamento
pulito, coerente con un genoma di riferimento piccolo e ben caratterizzato
come *E. coli* e con dati pre-trimmati con fastp.

## Dove vengono generati

- Regola Snakemake: `plot_bamstats` in `pipelines/rna-seq/Snakefile`
- Dipendenze: `samtools` (per `samtools stats`), `gnuplot` e `perl-uri`
  (richiesti da `plot-bamstats`) — dichiarate in
  `pipelines/rna-seq/envs/rna-seq.yml`
- Output: `results/rna-seq/{run}/aligned/{sample}_bamstats_plots/index.html`
