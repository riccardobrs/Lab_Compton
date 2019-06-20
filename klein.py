# python klein.py

import ROOT as rt
from math import *

Ei = '(511000*1.60217649e-19)'
den = '(1+('+Ei+'/((9.10938215e-31)*299792458*299792458))*(1-cos(x*3.1415927/180)))'
Ef = '('+Ei+'/'+den+')'
r2 = '((2.817940297322e-15)*(2.817940297322e-15)*0.5)'
Ef_i = '('+Ef+'/'+Ei+')'
Ei_f = '('+Ei+'/'+Ef+')'
klein = r2+'*'+Ef_i+'*'+Ef_i+'*('+Ei_f+'+'+Ef_i+'-'+'sin(x*3.1415927/180)*sin(x*3.1415927/180))'

t = rt.TF1('t', klein,0,180)
c = rt.TCanvas('c', '', 100, 200, 700, 500)
g=rt.TGraphErrors()
g.SetPoint(0,30,2.31676e-24)
g.SetPoint(1,45,4.99561e-26)
g.SetPoint(2,60,1.08756e-26)
g.SetPoint(3,75,6.54125e-27)
g.SetPoint(4,130,3.17227e-26)
g.SetPointError(0,6,1.07508e-25)
g.SetPointError(1,6,4.64029e-27)
g.SetPointError(2,6,6.37536e-28)
g.SetPointError(3,6,9.80522e-28)
g.SetPointError(4,6,6.84152e-27)
#g.SetTitle("Sezione d'urto differenziale")
g.GetXaxis().SetTitle("#theta [#circ]")

g.GetYaxis().SetTitle("#frac{d#sigma}{d#Omega} [m^{2}]")

c.cd()
g.Draw('AP')
t.Draw('SAME')
c.SaveAs('klein.png')
c.SaveAs('klein.root')
