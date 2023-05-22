# -*- coding: utf-8 -*-
import numpy as np
import dm4bem
#PARTIE 1

Text = 5+273 #température extérieure
Tc = 20+273 #température imposée à l'intérieur de la maison

emext = 30*(10**(-3)) #épaisseur mur extérieur en m
emint = 10*(10**(-3)) #épaisseur mur intérieur en m
ep = 10*(10**(-3)) #épaisseur porte en m
ef = 3*(10**(-3)) #épaisseur fenetre double vitrage en m
Spf = 2 #surface fenetre et porte
Sm1 = 30-2*Spf #surface mur 1
Sm = 30-Spf #surface mur 2 et 3

#conductivité bois
lambois = 0.077 
hconvext = 10 #coef de convection à l'extérieur
hconvint = 2 #coef de onvection à l'intéieur

#coductivité thermique
Gconvext1 = hconvext*Sm1 #sur le mur 1 à l'extérieur
Gconvint1 = hconvint*Sm1 #sur le mur 1 à l'intérieur
Gconvext2 = hconvext*Sm #sur le mur 2 à l'extérieur
Gconvint2 = hconvint*Sm #sur le mur 2 à l'intérieur

Gcondm1 = lambois*Sm1*2/emext #conduction de la moitié de mur 1
Gcondm2 = lambois*Sm*2/emext #conduction de la moitié de mur 2

Gint = lambois*Spf*2/ep + 1/((2/(Sm*hconvint))+lambois*Sm*2/emint)#conductance du mur entre les 2 pièces

#advection
Va = 10*10*3
ACH = 1
Va_dot = ACH/3600 * Va
rhoair = 1.3
Cair = 1

Gad = rhoair*Cair*Va_dot #conductivité advection


Gc = 10**4 #gain proportionnel controlleur

Up = 2.5 #coef thermique de la porte
Uf = 1.5 #coef thermique de la fenetre
G1 = Up*Spf+Uf*Spf #conductance porte+fenetre mur 1
G2 = Uf*Spf#conductance fenetre mur 2

#capacité
mv = 370 #kg/m3
massem1 = mv*Sm1*emext #masse du mur extérieur 1
massem2 = mv*Sm*emext #masse du mur extérieur 1
Csp = 1250 #capacité spécifique
C1 = massem1*Csp
C2 = massem2*Csp

#PARTIE 2

B = np.zeros(13)
B[0]=Text
B[4]=Text
B[9]=Text
B[10]=Text
B[11]=Tc
B[12]=Tc
print(B)

A = np.zeros((13,8))
A[0,0]=1
A[1,0]=-1
A[1,1]=1
A[2,1]=-1
A[2,2]=1
A[3,2]=-1
A[3,3]=1
A[4,3]=1
A[5,3]=-1
A[5,4]=1
A[6,4]=1
A[6,5]=-1
A[7,5]=1
A[7,6]=-1
A[8,6]=1
A[8,7]=-1
A[9,7]=1
A[10,4]=1
A[11,3]=1
A[12,4]=1
print(A)

G=np.diag([Gconvext1,Gcondm1,Gcondm1,Gconvint1,G1,Gint,Gconvint2,Gcondm2,Gcondm2,Gconvext2,G2,Gc,Gc])
print(G)

C = np.zeros((8,8))
C[1,1]=C1
C[6,6]=C2
print(C)

f = np.zeros(8)
Q1= 300 #émissions chaleur gens et appareils dans la pièce 1 (1 personne=80W)
Q2= 300 #émissions chaleur gens et appareils dans la pièce 2
E = 100 #éclairement du soleil en W/m2
phi1 =E*Sm1*0.8 #corps gris
phi2 =E*Sm*0.8 #corps gris

f[0]=phi1
f[3]=Q1
f[4]=Q2
f[7]=phi2
#ajouter les radiations si possible après
print(f)

#Résolution en régime permanent donc theta_dot = 0. On peut alors résoudre :
AT = np.transpose(A)
INV = np.linalg.inv(np.dot(np.dot(AT,G),A))
ATGB = np.dot(np.dot(AT,G),B)
theta0 = np.dot(INV,ATGB+f)
print(theta0-273) #en degré

U = np.zeros(10)
U[0]=Text
U[1]=Text
U[2]=Text
U[3]=Text
U[4]=Tc
U[5]=Tc
U[6]=phi1
U[7]=phi2
U[8]=Q1
U[9]=Q2
Y = np.zeros(8)
Y[3]=1
Y[4]=1
[As, Bs, Cs, Ds] = dm4bem.tc2ss(A, G, B, C, f, Y)
theta = -np.dot(np.dot(np.linalg.inv(As),Bs),U)
print(theta-273) #en degré
#résultat : theta = theta0 (2 manière de le calculer)

#PARTIE 3 : modèle en dynamique

#trouver le pas de temps
eig = np.linalg.eig(As)[0]
dtmax = min(-2/eig) 
print(dtmax) #en s
