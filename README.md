# Esperimentazioni di Fisica Nucleare e Subnucleare

## Effetto Compton

In questa repository sono contenuti i codici utilizzati per le analisi dati e per la relazione finale.

## Librerie:

* myLib.h
* myLib.cc

## Programmi:

* graph.cpp: riproduce in un istogramma i dati acquisiti.
  Utile per capire i range di fit per i 2 picchi.
  Necessita un file *.txt contenente i dati in argv[1]
* fit.cpp: come graph.cpp ma esegue i fit chiedendo i set da input.
  Utile per vedere se i range di fit sono adeguati.
* risol.cpp: riproduce un grafico risoluzione vs V_bias.
