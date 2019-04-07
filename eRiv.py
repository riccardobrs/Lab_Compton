#!/usr/bin/env python
# -*- coding: utf-8 -*- 

# commenti precedenti necessari per il corretto encoding del file di input
# non è necessaria la compilazione. Eseguire da terminale: python eRiv.py energie.txt
# il file in argv[1] deve essere così strutturato:  |Energia(keV)|Channel|Channel_err|

import os
import sys # per accedere ad argv[1]
import ROOT as rt

gr = rt.TGraphErrors() # se avessi inserito soltanto "import ROOT" dovrei scrivere ROOT.TGraphErrors()
gr.SetMarkerSize(0.4)
gr.SetMarkerStyle(20)
gr.SetMarkerColor(4)
gr.GetXaxis().SetTitle('E (keV)')
gr.GetYaxis().SetTitle('Channel')

filein = open(sys.argv[1], 'r') # ozpione 'r' = modalità 'read'
row = filein.readlines() # leggo il contenuto di una intera riga

i = 0

for line in row: # è fondamentale l'indentazione (o equivalentemente 4 spazi)
    var = line.split() # con il metodo 'split' creo una lista contenente i valori relativi ad ogni colonna
    gr.SetPoint(i, float(var[0]), float(var[1]))
    gr.SetPointError(i, 0., float(var[2]))
    if i == 0:
        xmax = float(var[0])
    else:
        xmax = max(xmax, float(var[0]))
    i += 1
    
filein.close()

c = rt.TCanvas('Energia', '', 100, 200, 700, 500)
f = rt.TF1('fitfunc', '[0]+[1]*x', 0., xmax+0.1*xmax)
f.SetLineColor(2)
#f.SetLineWidth(2)

c.cd()
gr.Draw('AP')
gr.Fit('fitfunc', 'RS')
c.SaveAs('energia.png')