#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  2 16:16:46 2022

@author: akel
"""
import numpy as np
r=2
arc=3.13
As=np.pi*r**2
A_arc=arc*r**2    #area do arco formado)
y=r*np.sin(arc/2) #Altura do triangulo central)
x=r*np.cos(arc/2) #base do triangulo central
St=(x*y/2)*2            #Área
Ai=A_arc-2*St     #interseção
print('angulo:', 2.501*180/np.pi)
print('mag:',1-x/r)
print('OBS:',Ai/As*100,'%')
