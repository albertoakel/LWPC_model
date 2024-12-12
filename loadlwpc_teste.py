#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:44:42 2022

@author: akel
Leitura log LWPM
"""
import numpy as np


#localizador PRFL_HGTS:

filename='lwpm.log'
arq=open(filename,"rt")

for i in range(20):  ##linhas de cabeçalhos
    arq.readline()
line1=arq.readline()
temp=np.array(line1.split(),dtype=float)
v=temp[1:10]
print('1-',temp[1])
n=1
while True:
    n+=1
    line=arq.readline()
    if not line:
        break
    line=line[:-1]
    #line=arq.readline()
    if line==str(' '):
        print('Ué')
        break
    temp=np.array(line.split(),dtype=float)
    if len(temp)==9:
        v=np.append(v,temp)
        print(n,temp[0])

    
n=int(len(v)/9)
V=np.reshape(v,(n,9))    
    
    
    
    # (arq.readline() !=' \n'):
    # n+=1
    # line=arq.readline()
    # temp=np.array(line.split(),dtype=float)
    # #print(n,'->',temp[0])
    # if len(temp)==9:
    #     n+=1
    #     print(n,'->',temp[0])
    #     v=np.append(v,temp)
   


#V=np.reshape(v,(n,9))    
    #print(line)
    #print('n',n)
    #print(line)

#line2=arq.readline()
#line3=arq.readline()


#    print(arq.readline())
#while (name )
# for i in range(348):
#     line=arq.readline()
#     temp=np.array(line.split(),dtype=float)
#     if len(temp)==9:
#         n+=1
#         v=np.append(v,temp)

# V=np.reshape(v,(n,9))

# line3=arq.readline()
# print(line3)



# while (i <= 5):
#     print(i)
#     i += 1
