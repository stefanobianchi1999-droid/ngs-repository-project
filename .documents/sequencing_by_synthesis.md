# Come funziona il sequenziamento Illumina (sequencing-by-synthesis)

Note di studio: come una libreria di frammenti diventa le read che poi
vediamo nei FASTQ e nei grafici QC di `plot-bamstats`
([samtools_bamstats_plots.md](samtools_bamstats_plots.md)).

## 1. Library prep: adapter e indici

Il DNA/cDNA estratto viene frammentato. A **entrambe le estremità** di ogni
frammento viene legata (ligation) una sequenza adapter nota e universale,
uguale per tutte le librerie Illumina, che include anche un **indice di
campione** (barcode):

```
[adapter P5] -- [indice campione] -- [frammento biologico] -- [indice] -- [adapter P7]
```

L'indice serve solo a **multiplexing/demultiplexing**: permette di mescolare
più campioni nella stessa run e poi, dopo il sequenziamento, smistare le read
nel FASTQ del campione giusto leggendo cicli dedicati (index read,
separati dai cicli di R1/R2). Non ha alcun ruolo nella scelta di quale base
si incorpora durante la lettura del frammento biologico.

## 2. Caricamento sulla flow cell

La superficie della flow cell è ricoperta da oligonucleotidi corti ancorati
chimicamente (P5 e P7), complementari agli adapter. Quando la libreria viene
versata sulla flow cell, ogni frammento si appaia (ibridazione) a un oligo
complementare sulla superficie e resta ancorato in un punto preciso.

## 3. Bridge amplification → cluster

Il frammento ancorato si piega verso un oligo vicino complementare all'altro
adapter, la polimerasi lo copia, e il processo si ripete molte volte nello
stesso punto. Il risultato è un **cluster**: non una singola molecola, ma
**migliaia di copie identiche** dello stesso frammento concentrate in un
unico punto della flow cell.

## 4. Sequencing primer

Prima dei cicli di lettura vera e propria, un primer di sequenziamento si
appaia a una porzione dell'adapter su ogni copia del cluster. È questo
primer — non il frammento originale — a fornire l'estremità 3' libera da cui
la polimerasi parte ad aggiungere basi, un ciclo alla volta.

## 5. Un ciclo di sequenziamento

Ad ogni **ciclo**:

1. I 4 nucleotidi (A/C/G/T) vengono versati insieme su tutta la flow cell.
   Ognuno è marcato con un fluoroforo di colore diverso e porta un gruppo
   chimico bloccante ("terminatore reversibile") sulla posizione 3'.
2. Su ogni singola copia di ogni cluster si lega **solo la base
   complementare** al filamento stampo (appaiamento Watson-Crick: A-T,
   C-G) — le altre 3 non formano un legame stabile e vengono lavate via.
3. Anche la base corretta, una volta incorporata, non permette
   l'aggiunta di una seconda base nello stesso ciclo: il terminatore 3'
   blocca fisicamente il sito di aggancio finché non viene rimosso
   chimicamente.
4. Una fotocamera fotografa l'intera flow cell: il colore emesso da
   ciascun cluster in quello scatto indica quale base è stata appena
   aggiunta. Il segnale è la somma della fluorescenza di migliaia di
   copie identiche che si accendono sincronizzate dello stesso colore —
   una singola molecola non sarebbe rilevabile.
5. Il terminatore viene rimosso chimicamente da tutte le basi appena
   aggiunte e si passa al ciclo successivo.

Ogni ripetizione di questi passaggi è un **cycle**, e corrisponde a **una
posizione all'interno della read** (1ª base, 2ª base, ... fino alla
lunghezza massima, es. 150). Non è una coordinata sul genoma: è per questo
che l'asse X dei grafici "Quality per cycle" e "Per-base sequence content"
in `plot-bamstats` va letto come "posizione nella read", aggregata su tutte
le read del campione.

## Perché la qualità cala nei cicli finali

Ad ogni ciclo si accumula un po' di rumore chimico/ottico: alcune copie nel
cluster perdono sincronia con le altre (restano "indietro" o "avanti" di un
ciclo), i fluorofori sbiadiscono, il segnale si indebolisce. È un effetto
fisiologico della chimica Illumina — un calo progressivo verso la fine della
read è normale, un crollo brusco e precoce indica un problema di
sequenziamento.

## Schema visivo

Diagramma pubblicato: vedi conversazione — schema con flow cell, i 4
passaggi di un ciclo, e la corrispondenza ciclo → colonna della read.
