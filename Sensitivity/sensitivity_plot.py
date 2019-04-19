# Programma per disegnare un plot di esclusione (Valerio Gentile)

print "\nExclusion plot for WIMP"
print "\nInserire i seguenti parametri:\n"

# Importo dei moduli
import math
import ROOT
from ROOT import TCanvas, TGraph
from array import array

# Scrive un file
out_file = open("nit1000.txt","w")


# Inizializzazione degli indici
i=0
j=0
k=0
l=0
n=0

# Input iniziali
Nel = 8
#Nel = input("Numero di elementi del target: ")                                    # Num elementi del target
#print
el = [s for s in range(Nel)]                                                      # Array per indicizzazione degli elementi
#A = [str(input("Numero di massa del target " + str(u+1) + ": ")) for u in el]     # Array per il numero di massa di ciascun elemento
#print 
#f = [str(input("Fraction mass del target " + str(u+1) + ": ")) for u in el]       # Array per la frazione in massa di ogni elemento
#A2 = [0]*len(A)
#f2 = [0]*len(A)
#while l < len(A):                                                                 # Trasformo gli array di stringhe in array numerici 
#    A2[l] = float(A[l])
#    f2[l] = float(f[l])
#    l +=1

A2 = [108,80,127,12,14,16,1,32]
f2 = [0.44,0.32,0.019,0.101,0.027,0.074,0.016,0.003]
    
Kg = input("Quanti chili: ")                                                      # Chili del target 
day = input("Quanti giorni: ")                                                    # Giorni di esposizione


# Definizione delle variabili
pigreco = math.pi
densdm = 0.4                                      # (GeV)(c^-2)(cm^-3)      (densita di DM)
numav = 6.02*pow(10,26)                           # (Kg^-1)                 (num Avogadro)
v0 = 230.0                                        # (Km/s)                  (velocita media di una wimp)  
vE = 244.0
crsect_norm = pow(10,-36)                         # (pb)                    (per normalizzare R0)
vesc = 600.0                                      # (Km/s)                  (velocita di fuga della galassia)
k0suk1 = 0.9965                                   #                         (costante di normalizzazione)
c1 = 0.751                                        #                         (parametri fattore di forma Helm SI)    
c2 = 0.561                                        #                         (      "         "         "       )
s = 0.9                                           # (fm)                    (      "         "         "       )
rn = [1.14*pow(x,1./3.) for x in A2]                # (fm)                    (      "         "         "       )                         
rate = 0.0                                        # (Kg^-1)(d^-1)           (Rate di eventi integrato sull'energia per ciascun elemento)
s_in_d = 86400                                    # (d)                     (Secondi in giorni)
Mt = [0.932*x for x in A2]                        # (GeV/c^2)               (massa target) 
sum_rate = 0.0                                    # (Kg^-1)(d^-1)           (Rate totale di eventi su tutto il target)

c = 299792.458

# creazione degli array
dE=1
E = range(1,1500,1)                             # Array per le energie di rinculo
#E = [x*pow(10,-6) for x in E]                    # Porta le energie in KeV
M1 = range(10,100)                                # Array per la massa delle WIMP
M2 = range(10,100)
M3 = range(100,1000,10)
crsect = [0]*3*len(M1)                            # Array per la sezione d'urto
coupling = pow(10,-41)

a = math.erf(vesc/v0)
b = (2*vesc)/(pow(math.pi,0.5)*v0)
d = math.exp(-pow(vesc/v0,2))

RperF = 0
point = 0
sigma = 0

k0k1 = pow(a-b*d,-1)
#print k0k1
#print rn

g = ROOT.TGraph()
#print E                                            

# Ciclo per il calcolo della sezione d'urto al variare di Md  (1-10000 GeVc^-2)
while k < 3*len(M1):
    i=0
    sum_rate=0
    if(k<len(M1)):
        Md = M1[k]*0.1
    if(k<2*len(M1) and k>=len(M1)):
        Md = M2[k-len(M1)]
    if(k>=2*len(M1)):
        Md = M3[k-2*len(M1)]        
    E0 = 0.5*Md*pow(v0/c,2)*pow(10,6)
    #E0 = 0.5*Md*pow(v0,2)*pow(10,4)/9
    #print v0/c, "\n"
    #print Md

# Ciclo per ogni elemento che compone il target
    while i<len(el):                     
        j = 0 
        #print i
        rate = 0.0
        Mtar = Mt[i]
        Atar = A2[i]
        ftar = f2[i]
        rnt = rn[i]

        
        MN = (Md*0.932*Atar)/(Md+Atar*0.932)
        Mp = (Md*0.938)/(Md+0.938)
        xsec = pow(Atar*(MN/Mp),2)#*coupling
        #print xsec,"\n"
        r = (4*Md*Mtar)/(pow(Md+Mtar,2))
        R0 = (2./pow(pigreco,0.5))*(numav/Atar)*(densdm/Md)*v0*pow(10,5)*s_in_d*xsec
        #R0 = (0.5/pow(pigreco,0.5))*(numav/Atar)*(densdm/Md)*(v0*pow(10,5)*s_in_d)
      
        #print E0, "\t",r,"\t",R0,"\t", R0/(E0*r),"\n"

        #print E0, "\t",Md,"\t",v0,"\n"

       
         
        # Ciclo per l'integrazione sullo spettro di energia  (10-100 KeV)
        while j < len(E):

             Er = E[j]
             # print Er, "\n"
             #qr_n = [6.92*pow(10,-3)*pow(Atar,0.5)*pow(x,0.5)*rnt for x in E]          # q adimensionale
             #qs = [6.92*pow(10,-3)*pow(Atar,0.5)*pow(x,0.5)*s for x in E]              # q adimensionale
             #F = [(3*(math.sin(x)-x*math.cos(x))/pow(x,3))*math.exp(-pow(y,2)/2) for x,y  in zip(qr_n,qs)]
             #F = [(3*(math.sin(x)-x*math.cos(x))/pow(x,3)) for x in qr_n]

             qr_n =6.92*pow(10,-3)*pow(Atar,0.5)*pow(Er,0.5)*rnt          # q adimensionale
             #Ftar = F[j]
             Ftar = (3*(math.sin(qr_n)-qr_n*math.cos(qr_n))/pow(qr_n,3))
           

             vmin = v0*pow(Er/(E0*r),0.5)
             #print qr_n[j], " ",qr_n2," 1\n"
             #print Ftar," ",Ftar2," 2\n"
             # print Ftar, "\n"

             #exp1 = c1*math.exp(-(c2*Er)/(E0*r))#-())
             exp1 = (R0/(E0*r))*(pow(pigreco,0.5)/4.)*(v0/vE)*(math.erf((vmin+vE)/v0)-math.erf((vmin-vE)/v0))
             #print exp1," ",E0*r, "\n"
             #exp2 = math.exp(-pow(vesc/v0,2))

             #Rate2 = (k0suk1)*(R0/(E0*r))*(exp1- exp2)
             Rate2 = exp1
            # print Md," ",Er," ",i," ", R0," ", exp1, " ", exp2," ",qr_n[j]," ",rnt,"\n"
             
             #RperF = Rate2*Ftar*Ftar*ftar                                   # Moltiplica il rate ed il fattore di forma fissati Md e Er

             
             if(i==0 and Er >= 245.1 and Er <=2332.8):
                 RperF = Rate2*Ftar*Ftar*ftar*dE                             # Moltiplica il rate ed il fattore di forma fissati Md e Er
                 if(RperF>0):
                     rate += RperF                    
                # print i," ", Er," ", rate, " ", Rate2, "\n"
                 
             if(i==1 and Er >= 179 and Er <=1683.5):
                 RperF = Rate2*Ftar*Ftar*ftar*dE                                # Moltiplica il rate ed il fattore di forma fissati Md e Er
                 if(RperF>0):
                     rate += RperF 
             if(i ==2 and Er >= 273.4 and Er <=2750.0):
                 RperF = Rate2*Ftar*Ftar*ftar*dE                                 # Moltiplica il rate ed il fattore di forma fissati Md e Er
                 if(RperF>0):
                     rate += RperF 
             if(i ==3 and Er >= 34.9 and Er <=444.1):
                 RperF = Rate2*Ftar*Ftar*ftar*dE                                  # Moltiplica il rate ed il fattore di forma fissati Md e Er
                 if(RperF>0):
                     rate += RperF
                     #if(Md==30):
                         #print Ftar*Ftar,"\t",Er,"\t",E0,"\t",R0,"\t",(math.erf((vmin+vE)/v0)-math.erf((vmin-vE)/v0)),"\t",(R0/(E0*r))*(pow(pigreco,0.5)/4.)*(v0/vE),"\t",Rate2,"\n"
             if(i ==4 and Er >= 41.1 and Er <=501.1):
                 RperF = Rate2*Ftar*Ftar*ftar*dE                                  # Moltiplica il rate ed il fattore di forma fissati Md e Er
                 if(RperF>0):
                     rate += RperF 
             if(i ==5 and Er >= 45.3 and Er <=518.1):
                 RperF = Rate2*Ftar*Ftar*ftar*dE                                  # Moltiplica il rate ed il fattore di forma fissati Md e Er
                 if(RperF>0):
                     rate += RperF
             """        
             if(i==6 and Er >= 7.6 and Er <=121.9):
                 RperF = Rate2*Ftar*Ftar*ftar*dE                                  # Moltiplica il rate ed il fattore di forma fissati Md e Er
                 if(RperF>0):
                     rate += RperF 
             """
             if(i ==7 and Er >= 89 and Er <=918.8):
                 RperF = Rate2*Ftar*Ftar*ftar*dE                                  # Moltiplica il rate ed il fattore di forma fissati Md e Er
                 if(RperF>0):
                     rate += RperF             
             
             #rate += RperF
             

             j+=1  
        
        #print Rate, "\t"
        sum_rate += rate                                                       # Somma del rate integrato per ciascun elemento
        #print sum_rate, "\t", Md, "\n"
        i +=1
    #sigma = 3/(sum_rate*day*Kg-math.exp(-(pow(vesc/v0,2))))                     # Formula inversa per il calcolo della sezione d'urto
    if(sum_rate>0):
        sigma = 2.44/(sum_rate*Kg*day)
        crsect[k] = sigma
        g.SetPoint(point,Md,sigma)
        point +=1
    k +=1
    #print Md," ", sigma," ",sum_rate,"\n"
    out_file.write(str(Md)+" "+str(sigma)+"\n")

#creazione del grafico
nBins = len(M1)
#print nBins

#g = ROOT.TGraph(nBins, array('d', M), array('d', crsect))
c = TCanvas()
g.SetMarkerStyle(0)
g.SetMarkerSize(0)
ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()
g.SetTitle("Exclusion Plot")
g.GetXaxis().SetTitle("WIMP Mass [GeV/c^{2}]")
g.GetYaxis().SetTitle("WIMP-nucleon cross section[cm^{2}]")
g.GetXaxis().CenterTitle()
g.GetYaxis().CenterTitle()
g.GetXaxis().SetTitleOffset(1.4);
g.GetYaxis().SetTitleOffset(1.4);
g.Draw("APL")
g.SetLineColor(4)
g.SetLineWidth(3)


      
out_file.close()
raw_input("Press enter to exit")




