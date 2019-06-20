import ROOT as r
from math import *

def gaussian(x, par):
    return par[0] * exp( -( (x[0] - par[1])**2) / ( 2*(par[2]**2)))

g=r.TGraph()
g.SetPoint(0,0,0)
g.SetPoint(1,5000,0)

a1 = '4469.34' 
a1_err = '4.86941'
m1 = '3399.62' 
m1_err = '0.104307'
s1 = '114.789' 
s1_err = '0.0786632'
P0 = '(1/(sqrt(2*3.1415927)*'+s1+'))*exp(-(x-'+m1+')*(x-'+m1+')/(2*'+s1+'*'+s1+'))'

a2 = '85.887' 
a2_err = '4.02752'
m2 = '3012.35' 
m2_err = '4.95253'
s2 = '140.173' 
s2_err = '7.45731'
P30 = '(1/(sqrt(2*3.1415927)*'+s2+'))*exp(-(x-'+m2+')*(x-'+m2+')/(2*'+s2+'*'+s2+'))'

a3 = '112.223' 
a3_err = '2.97132'
m3 = '2594.8' 
m3_err = '3.85113'
s3 = '136.363'
s3_err = '4.35201'
P45 = '(1/(sqrt(2*3.1415927)*'+s3+'))*exp(-(x-'+m3+')*(x-'+m3+')/(2*'+s3+'*'+s3+'))'

a4 = '86.5305' 
a4_err = '7.82837'
m4 = '2225.18' 
m4_err = '11.1355'
s4 = '109.668' 
s4_err = '13.4881'
P60 = '(1/(sqrt(2*3.1415927)*'+s4+'))*exp(-(x-'+m4+')*(x-'+m4+')/(2*'+s4+'*'+s4+'))'

a5 = '320.282' 
a5_err = '37.0501'
m5 = '2037.26' 
m5_err = '7.24593'
s5 = '167.536' 
s5_err = '17.1476'
P75 = '(1/(sqrt(2*3.1415927)*'+s5+'))*exp(-(x-'+m5+')*(x-'+m5+')/(2*'+s5+'*'+s5+'))'

a6 = '1410.09' 
a6_err = '1052.85'
m6 = '1261.8' 
m6_err = '7.08438'
s6 = '263.515'
s6_err = '56.2324'
P130 ='(1/(sqrt(2*3.1415927)*'+s6+'))*exp(-(x-'+m6+')*(x-'+m6+')/(2*'+s6+'*'+s6+'))'

c = r.TCanvas("c", "", 100, 200, 700, 500)
t0 = r.TF1('0^{#circ}', P0, 0,5000)
t30 = r.TF1('30^{#circ}', P30, 0,5000)
t45 = r.TF1('45^{#circ}', P45, 0,5000)
t60 = r.TF1('60^{#circ}', P60, 0,5000)
t75 = r.TF1('75^{#circ}', P75, 0,5000)
t130 = r.TF1('130^{#circ}', P130, 0,5000)
c.cd()
g.Draw('AP')
t0.Draw('SAME')
t30.Draw('SAME')
t45.Draw('SAME')
t60.Draw('SAME')
t75.Draw('SAME')
t130.Draw('SAME')

c.SaveAs('picchi.root')



