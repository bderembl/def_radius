#!/usr/bin/env python

import numpy as np
from scipy import linalg as la
from scipy.sparse import diags

# def construct_mat(H,gp,fo) :

#   si_z = len(H)
#   H = 1.0*H
  
#   mata = np.zeros((si_z,si_z));
  
#   mata[0,0] =  1/(H[0]*gp[0]);
#   mata[0,1] = -1/(H[0]*gp[0]);
  

#   if si_z > 2 :
#     for nz in range(1,si_z-1):
#       mata[nz,nz] = 1/H[nz]*(1/gp[nz]+1/gp[nz-1]);
#       mata[nz,nz-1] = -1/(H[nz]*gp[nz-1]);
#       mata[nz,nz+1] = -1/(H[nz]*gp[nz]);
      
#   mata[si_z-1,si_z-1]  =  1/(H[si_z-1]*gp[si_z-2]);
#   mata[si_z-1,si_z-2]  = -1/(H[si_z-1]*gp[si_z-2]);
      
#   mata = fo**2*mata;

#   return mata

#####################################################
def gamma_stretch(dh,N2,fo,**kwargs) :
  
  gp = N2*0.5*(dh[1:] + dh[:-1])
  construct_mat(dh,gp,fo, kwargs)

  return mata


#####################################################
def construct_mat(dh,gp,fo,**kwargs) :
  
  sparse = kwargs.get('sparse', False)

  diag_p1 = -1/(dh[:-1]*gp)
  diag_m1 = -1/(dh[1:]*gp)
  
  diag0 = -np.append(diag_p1,0.)
  diag0[1:] = diag0[1:] - diag_m1
  
  diagonals = [diag0,diag_m1,diag_p1]
  mata = diags(diagonals, [0, -1, 1])

  mata = fo**2*mata

  if not sparse:
    mata = mata.toarray()

  return mata

#####################################################
def cal_rad(dh,gp,fo) :
  
  mata = construct_mat(dh,gp,fo)

  iRd2 = la.eig(mata,right=False)

  iRd2 = iRd2.real
  iRd2 = np.sort(iRd2)
  with np.errstate(divide='ignore', invalid='ignore'):
    Rd = 1./np.sqrt(iRd2[:])

  return Rd

###########################################################
def cal_transfo(dh,gp,fo) :

  # compute matrices for the mode/layer conversion
  # m_m2l[:,0] is the barotropic mode: should be 1..1
  # m_m2l[:,i] is the ith baroclinic mode

  # to convert from physical to modal
  #  u_mod = np.dot(m_l2m[:,:],u_lev)

  # to go back to the physical space
  # u_lev = np.dot(m_m2l[:,:],u_mod)



  
  si_z = len(dh)
  mata = construct_mat(dh,gp,fo)
  
  iRd2, eigl,eigr= la.eig(mata,left=True)
    
  iRd2 = iRd2.real
  idx = np.argsort(iRd2)
  eigl = eigl[:,idx]
  eigr = eigr[:,idx]

  #  iRd2 = iRd2[idx]
  #  Rd = 1./np.sqrt(iRd2[:])/1000  
  
  # Normalize eigenvectors
  Ht = np.sum(dh)
  dh = np.reshape(dh,(si_z,1))
  
  scap = np.sum(dh*eigr*eigr,0)
  eigr = eigr*np.sqrt(Ht/scap)*np.sign(eigr[0,:])
  
  # scalar product
  check = np.sum(dh.T*eigr[:,0]*eigr[:,0])/Ht
  
  scap2 =  np.sum(eigl*eigr,0)
  eigl = eigl/scap2
  
  m_l2m = eigl.T
  m_m2l = eigr


  return m_l2m, m_m2l

###########################################################
