#!/usr/bin/env python

import numpy as np
from scipy import linalg as la

def construct_mat(H,gp,fo) :

  si_z = len(H)
  H = 1.0*H
  
  mata = np.zeros((si_z,si_z));
  
  mata[0,0] =  1/(H[0]*gp[0]);
  mata[0,1] = -1/(H[0]*gp[0]);
  

  if si_z > 2 :
    for nz in range(1,si_z-1):
      mata[nz,nz] = 1/H[nz]*(1/gp[nz]+1/gp[nz-1]);
      mata[nz,nz-1] = -1/(H[nz]*gp[nz-1]);
      mata[nz,nz+1] = -1/(H[nz]*gp[nz]);
      
  mata[si_z-1,si_z-1]  =  1/(H[si_z-1]*gp[si_z-2]);
  mata[si_z-1,si_z-2]  = -1/(H[si_z-1]*gp[si_z-2]);
      
  mata = fo**2*mata;

  return mata

#####################################################
def cal_rad(H,gp,fo) :
  
  mata = construct_mat(H,gp,fo)

  iRd2 = la.eig(mata,right=False)

  iRd2 = iRd2.real
  iRd2 = np.sort(iRd2)
  
  Rd = 1./np.sqrt(iRd2[:])/1000  

  return Rd

###########################################################
def cal_transfo(H,gp,fo) :

  si_z = len(H)
  mata = construct_mat(H,gp,fo)
  
  iRd2, eigl,eigr= la.eig(mata,left=True)
    
  iRd2 = iRd2.real
  idx = np.argsort(iRd2)
  eigl = eigl[:,idx]
  eigr = eigr[:,idx]

  #  iRd2 = iRd2[idx]
  #  Rd = 1./np.sqrt(iRd2[:])/1000  
  
  # Normalize eigenvectors
  Ht = np.sum(H)
  H = np.reshape(H,(si_z,1))
  
  scap = np.sum(H*eigr*eigr,0)
  eigr = eigr*np.sqrt(Ht/scap)*np.sign(eigr[0,:])
  
  # scalar product
  check = np.sum(H.T*eigr[:,0]*eigr[:,0])/Ht
  
  scap2 =  np.sum(eigl*eigr,0)
  eigl = eigl/scap2
  
  m_l2m = eigl.T
  m_m2l = eigr


  return m_l2m, m_m2l

###########################################################
