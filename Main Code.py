# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 18:48:00 2020

@author: apoor
"""
###to clear variable in spyder console#####
from IPython import get_ipython
get_ipython().magic('reset -sf')
############################################



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from decimal import Decimal
import math
import xlsxwriter
import csv
import matplotlib.ticker as mticker
############################### Lifetime functions  #####################################################

class Lifetime:
    
    def __init__(self,sig_p, sig_n, Eg_minus_E_defects, Impurity, T, n0, N_defects,J0,Jl,Rs, sweep_dn=False):  #do all calulations
        #Eg=1.124  #### in eV
        N0_eeh=3.3e17  #### in cm-3
        N0_ehh=7e17  #### in cm-3

        kb=1.3807e-23     ###### in (J/K)
        k = 8.61758e-5 #[eV/K]
        q=1.6021765e-19   ######### in C

        self.Eg=1.1785-((9.025e-5)*T)-((3.05e-7)*(T**2))
        #self.Eg=1.170-((0.0004734*(T**2))/(T+636))
        self.E_defects = self.Eg - Eg_minus_E_defects
        self.ni=2.9135*(10**15)*(T**1.6)*np.exp((-self.Eg)/(2*k*T))

        self.V = np.arange(0.2, 0.95, 0.0002) 
        self.Vlist = self.V.tolist()   ################# NEEDS CORRECTION

        self.J01 = J0*(T/300)**3*(np.exp(-self.Eg/(k*T)))/(np.exp(-1.124/(k*300)))
        #self.J01 = J0*(np.exp((-q*self.Eg)/(k*T)))/(np.exp((-q*1.124)/(k*300)))

        Nd=n0
        Na=1
        p0 = self.ni**2/Nd
        dn1 = (- Nd + np.sqrt((Nd**2)+4*((self.ni**2)*np.exp(self.V/(k*T)))))/2
        n1=(n0 + dn1)
        p1=(p0 + dn1)        
        
        ###### Schenk   #############
        
      
        Ry = 1.655e-2
        a_ex = 3.7185e-7
        Tn = T*k/Ry
        
        #model parameters
        alpha_e = 0.5187 
        alpha_h = 0.4813 
        b_e = 8
        b_h = 1 
        c_e = 1.3346 
        c_h = 1.2365 
        d_e = 0.893 
        d_h = 1.153 
        g_e = 12 
        g_h = 4 
        h_e = 3.91 
        h_h = 4.2 
        j_e = 2.8585 
        j_h = 2.9307 
        k_e = 0.012 
        k_h = 0.19 
        p_e = 7/30
        p_h = 7/30 
        q_e = 0.75 
        q_h = 0.25
        
        #normalization to Bohr radius
        p_n = p1*a_ex**3 
        n_n = n1*a_ex**3 
        c = p_n + n_n 
        n_p = alpha_e*n_n + alpha_h*p_n 
        Nd_n = Nd*a_ex**3 
        Na_n = Na*a_ex**3 
        Nc = Nd_n + Na_n
        
        #exchange-correlation 
        xcn = ( (4*np.pi)**3*c**2*((48*n_n/(np.pi*g_e))**(1/3)+c_e*np.log(1+d_e*n_p**p_e)) + 8*np.pi*alpha_e/g_e*n_n*Tn**2 + np.sqrt(8*np.pi*c)*Tn**(5/2)  ) / ( (4*np.pi)**3*c**2 + Tn**3 + b_e*np.sqrt(c)*Tn**2 + 40*c**(3/2)*Tn   ) * Ry * 1000                                    
        xcp = ( (4*np.pi)**3*c**2*((48*p_n/(np.pi*g_h))**(1/3)+c_h*np.log(1+d_h*n_p**p_h)) + 8*np.pi*alpha_h/g_h*p_n*Tn**2 + np.sqrt(8*np.pi*c)*Tn**(5/2)  ) / ( (4*np.pi)**3*c**2 + Tn**3 + b_h*np.sqrt(c)*Tn**2 + 40*c**(3/2)*Tn   ) * Ry * 1000                                      
    
        #ionic
        iD = Nc *(1+c**2/Tn**3) / ( np.sqrt(Tn*c/(2*np.pi))*(1+ h_e*np.log(1+np.sqrt(c)/Tn)) + j_e*c**2/Tn**3*n_p**(3/4)*(1+k_e*n_p**q_e)  )  *Ry*1000
        iA = Nc *(1+c**2/Tn**3) / ( np.sqrt(Tn*c/(2*np.pi))*(1+ h_h*np.log(1+np.sqrt(c)/Tn)) + j_h*c**2/Tn**3*n_p**(3/4)*(1+k_h*n_p**q_h)  )  *Ry*1000
    
        self.dEg = (xcn + xcp + iD + iA)*1e-3


          

    ##################################################################################


        
        #############################################################
        self.ni_eff=self.ni*np.exp((self.dEg*q)/(2*kb*T))
        self.dn = (- n0 + np.sqrt(n0**2+4*(self.ni_eff**2*np.exp(self.V/(k*T)))) )/2
        
        if sweep_dn:
            self.dn = np.logspace(10,17,1775)
        
        self.p0=(self.ni_eff**2)/n0 
        self.geeh=1+13*(1-np.tanh((n0/N0_eeh)**0.66))
        self.gehh=1+7.5*(1-np.tanh((self.p0/N0_ehh)**0.63))
        
        self.Vth_e=0.06761*(T**3)-98.513*(T**2)+72696*T+5709700
        self.Vth_h=-25.507*(T**2) + 45810*T + 5E+06
        
        self.tau_n0=1/(self.Vth_e*N_defects*sig_n)
        self.tau_p0=1/(self.Vth_h*N_defects*sig_p)
        
#        self.tau_n0=1/(self.Vth_e*N_defects*sig_n)
#        self.tau_p0=1/(self.Vth_h*N_defects*sig_p)
#      
#
#        self.Vth=0.06761*(T**3)-98.513*(T**2)+72696*T+5709700
#        self.tau_n0=1/(self.Vth*N_defects*sig_nT)
#        self.tau_p0=1/(self.Vth*N_defects*sig_pT)

        self.Nc=2.86*(10**19)*((T/300)**1.58)
        self.Nv=3.1*(10**19)*((T/300)**1.85)
        self.n1=self.Nc*np.exp((-q*Eg_minus_E_defects)/(kb*T))
        self.p1=self.Nv*np.exp((-q*self.E_defects)/(kb*T))


 ############### B(Alterman 2005)  ####################

        bmax=1
        rmax=0.2
        rmin=0
        smax=1.5e18
        smin=1e07
        wmax=4e18
        wmin=1e9
        b2=0.54
        b4=1.25
        r1=320
        r2=2.5
        s1=550
        s2=3
        w1=365
        w2=3.54


        b1=smax+(smin-smax)/(1+(T/s1)**s2)
        b3=wmax+(wmin-wmax)/(1+(T/w1)**w2)
        bmin=rmax+(rmin-rmax)/(1+(T/r1)**r2)
   
        self.Blow=10**(-9.65614-((8.05258e-2)*T)+(6.02695e-4*(T**2))-((2.29844E-6)*(T**3))+(4.31934E-9*(T**4))-(3.16154E-12*(T**5)))
        
        self.n=n0+self.dn
        self.p=self.p0+self.dn
        self.B=self.Blow*(bmin+(bmax-bmin)/(1+((self.n+self.p)/(2*b1))**b2+((self.n+self.p)/(2*b3))**b4))
        self.tau_rad=(self.dn*1000)/((self.n*self.p-(self.ni_eff**2))*self.B)
        self.tau_aug=(self.dn*1000)/((self.n*self.p-(self.ni_eff**2))*(2.5e-31*self.geeh*n0+8.5e-32*self.gehh*self.p0+3e-29*self.dn**0.92))
        self.tau_srh=(self.dn*1000)/((self.n*self.p-(self.ni_eff**2))/(self.tau_p0*(n0+self.n1+self.dn)+self.tau_n0*(self.p0+self.p1+self.dn)))
             
        self.tau_eff_no_surf=1/((1/self.tau_srh)+(1/self.tau_aug)+(1/self.tau_rad))
       
    ####################
        self.U_aug= ((self.n*self.p-(self.ni_eff**2))*(2.5e-31*self.geeh*n0+8.5e-32*self.gehh*self.p0+3e-29*self.dn**0.92))
        self.U_rad= ((self.n*self.p-(self.ni_eff**2))*self.B)
        self.U_srh=((self.n*self.p-(self.ni_eff**2))/(self.tau_p0*(n0+self.n1+self.dn)+self.tau_n0*(self.p0+self.p1+self.dn)))
        
    ####################
        self.Js = self.J01*(np.exp(self.V/(n02*k*T))-1)*1e3 #Surface recombination current
        self.Us = self.J01*(np.exp(self.V/(n02*k*T))-1)/q/width #rate
        
        self.tau_surf=(self.dn*1000)/self.Us
        self.tau_eff=1/((1/self.tau_srh)+(1/self.tau_aug)+(1/self.tau_rad)+(1/self.tau_surf))

            
        #Auger, radiative calculated according to Richter parameterization
        #SRH calculated according to SRH eqn above
        self.Jaug = q*width*self.U_aug*1e3 #Auger
        self.Jrad = q*width*self.U_rad*1e3 #Radiative
        self.Jsrh = q*width*self.U_srh*1e3 #SRH
        
        
        #self.Jint = q*width*self.U_int*1e3 #Intrinsic (Radiative + Auger)
        
        
        self.Jr = Jl-self.Jrad-self.Jaug-self.Jsrh-self.Js ####  Current density of everything
        ind = np.where(np.diff(np.sign(self.Jr)))[0] #Voc injection level index, where Jr crosses 0
        self.Voc = self.V[ind]*1000
        self.P = self.Jr*self.V
        ind2 = np.where(self.P==np.max(self.P)) #MPP injection level index
        self.Vmp=self.V[ind2]*1000
        self.dn_mpp = self.dn[ind2]#MPP injection level
        self.dn_oc=self.dn[ind]
        self.FF = (max(self.P)/(self.Voc/1000)/max(self.Jr))*100
        self.Eff = max(self.P)
        self.Jsc = max(self.Jr)
        self.Rch = self.Voc/(self.Jsc*1e-3);
        self.rs = Rs/self.Rch
        self.FFr = self.FF*(1-1.1*self.rs)+self.rs**2/5.4
        self.Effr = self.Jsc*self.Voc*self.FFr
        self.Eff_check = self.Jsc*self.Voc*self.FF;
                

        ###   Intrincisc + SRH-limited values
        self.Jr_no_surf=Jl-self.Jrad-self.Jaug-self.Jsrh
        ind_no_surf = np.where(np.diff(np.sign(self.Jr_no_surf)))[0] #Voc injection level index, where Jr crosses 0, Andre add -Js
        self.Voc_no_surf  = (self.V[ind_no_surf])*1000
        self.P_no_surf = (self.Jr_no_surf)*self.V # André add - Js
        ind2_no_surf = np.where(self.P_no_surf==np.max(self.P_no_surf)) #MPP injection level index
        self.Vmp_no_surf=(self.V[ind2_no_surf])*1000
        self.Eff_no_surf = (max(self.P_no_surf))
        self.FF_no_surf = (max(self.P_no_surf)/(self.Voc_no_surf/1000)/max(self.Jr_no_surf))*100
        self.dn_oc_no_surf=self.dn[ind_no_surf]
        self.dn_mpp_no_surf = self.dn[ind2_no_surf]
        
        
        ##### Intrinsic - aug+rad

        self.Jr_intrinsic=Jl-self.Jrad-self.Jaug
        ind_intrinsic = np.where(np.diff(np.sign(self.Jr_intrinsic)))[0] #Voc injection level index, where Jr crosses 0, Andre add -Js
        self.Voc_intrinsic  = (self.V[ind_intrinsic])*1000
        self.P_intrinsic = (self.Jr_intrinsic)*self.V # André add - Js
        ind2_intrinsic = np.where(self.P_intrinsic==np.max(self.P_intrinsic)) #MPP injection level index
        self.Vmp_intrinsic=(self.V[ind2_intrinsic])*1000
        self.Eff_intrinsic = (max(self.P_intrinsic))
        self.FF_intrinsic = (max(self.P_intrinsic)/(self.Voc_intrinsic/1000)/max(self.Jr_intrinsic))*100
        self.dn_intrinsic=self.dn[ind_intrinsic]
        self.dn_mpp_intrinsic = self.dn[ind2_intrinsic]
        
        
        
    

def plot_lifetime(tau_rad_list, tau_aug_list, tau_srh_list, tau_eff_no_surf_list, tau_eff_list,tau_surf_list, dn_list, plotlabel, savelabel, color,dash=False):
    f1 = plt.figure(1)
    #plt.plot(dn_list,tau_rad_list,'-',label = 'T ='"{:.0f}\u00b0C".format(param-273.15))
    if dash==False:
        plt.plot(dn_list,tau_rad_list,'-',label = plotlabel, color = color,linewidth=3)
    else:
        plt.plot(dn_list,tau_rad_list,'--',label = plotlabel, color = color,linewidth=3)
    plt.yscale('log')
    plt.xscale('log')
  #  plt.title("Radiative Lifetime vs Injection level \n/n("+savelabel+")")
    plt.ylabel('Radiative lifetime (ms)',fontsize = 19)
    plt.xlabel(r'Excess carrier density ($\mathregular{cm^{-3}}$)',fontsize = 18)
    
#    ax =plt.gca()
#    
#    ticks=[]
#    for i in ax.get_yticks().tolist():
#        if int(i)>0:
#            ticks.append('{:,}'.format(int(i)))
#        else:
#            ticks.append('{:,}'.format(float(i)))
#    ax.set_yticklabels(ticks)

    #ax.tick_params(direction = 'in', which = 'both')
    plt.legend(fontsize=16,labelspacing=0.3)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.ylim(1e-1,1e3)
#
    plt.xlim(1e14,1e17)
    f1.savefig("Radiative ("+savelabel+").png", bbox_inches='tight', bbox_to_anchor = (0,0), dpi = 300)

    
    f2 = plt.figure(2)
    if dash==False:
        plt.plot(dn_list,tau_aug_list,'-',label = plotlabel, color=color,linewidth=3)
    else:
        plt.plot(dn_list,tau_aug_list,'--',label = plotlabel, color=color,linewidth=3)
    plt.yscale('log')
    plt.xscale('log')
   # plt.title("Auger Lifetime vs Injection level ("+savelabel+")")
    plt.ylabel('Auger lifetime (ms)', fontsize = 19)
    plt.xlabel(r'Excess carrier density ($\mathregular{cm^{-3}}$)',fontsize = 18)
    
#    ax =plt.gca()
#    ticks=[]
#    for i in ax.get_yticks().tolist():
#        if int(i)>0:
#            ticks.append('{:,}'.format(int(i)))
#        else:
#            ticks.append('{:,}'.format(float(i)))
#    ax.set_yticklabels(ticks)

    #ax.tick_params(direction = 'in', which = 'both')

    plt.legend(fontsize=16,labelspacing=0.3)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.ylim(1e-1,1e3)
#
    plt.xlim(1e14,1e17)
    filename2='Auger ('+savelabel+').png'
    f2.savefig(filename2, bbox_inches="tight", bbox_to_anchor = (0,0), loc = 'best', dpi = 300)
    
    f3 = plt.figure(3)
    if dash==False:
        plt.plot(dn_list,tau_srh_list,'-',label = plotlabel, color=color,linewidth=3)
    else:
        plt.plot(dn_list,tau_srh_list,'--',label = plotlabel, color=color,linewidth=3)
    plt.yscale('log')
    plt.xscale('log')
  #  plt.title("SRH Lifetime vs Injection level ("+savelabel+")")
    plt.ylabel('Bulk SRH lifetime (ms)',fontsize = 19)
    plt.xlabel(r'Excess carrier density ($\mathregular{cm^{-3}}$)',fontsize = 18)

#    ax =plt.gca()
#    ticks=[]
#    for i in ax.get_yticks().tolist():
#        if int(i)>0:
#            ticks.append('{:,}'.format(int(i)))
#        else:
#            ticks.append('{:,}'.format(float(i)))
#    ax.set_yticklabels(ticks)

    #ax.tick_params(direction = 'in', which = 'both')
    plt.legend(fontsize=16,labelspacing=0.14,loc='lower right')
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    #plt.ylim(1e-1,1e3)
    plt.xlim(1e14,1e17)
    f3.savefig("SRH ("+savelabel+").png", bbox_inches='tight', bbox_to_anchor = (0,0), dpi = 300)
    
    f4 = plt.figure(4)
    if dash==False:
        plt.plot(dn_list,tau_eff_no_surf_list,'-',label = plotlabel, color=color,linewidth=3)
    else:
        plt.plot(dn_list,tau_eff_no_surf_list,'--',label = plotlabel, color=color,linewidth=3)
    plt.yscale('log')
    plt.xscale('log')
  #  plt.title("Effective lifetime (No surface) vs Injection level ("+savelabel+")\n")
    plt.ylabel('Effective lifetime (ms)',fontsize = 19)

#    ax =plt.gca()
#    ticks=[]
#    for i in ax.get_yticks().tolist():
#        if int(i)>0:
#            ticks.append('{:,}'.format(int(i)))
#        else:
#            ticks.append('{:,}'.format(float(i)))
#    ax.set_yticklabels(ticks)

    #ax.tick_params(direction = 'in', which = 'both')

    plt.xlabel(r'Excess carrier density ($\mathregular{cm^{-3}}$)',fontsize = 18)
    plt.legend(fontsize=16,labelspacing=0.3,loc='best')
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.ylim(1e-1,1e3)
    plt.xlim(1e14,1e17)
    f4.savefig("Effective (no surface) ("+savelabel+").png", bbox_inches='tight', bbox_to_anchor=(1, 0.5), dpi = 300)
    
    f5 = plt.figure(5)
    if dash==False:
        plt.plot(dn_list,tau_eff_list,'-',label = plotlabel, color=color,linewidth=3)
    else:
        plt.plot(dn_list,tau_eff_list,'--',label = plotlabel, color=color,linewidth=3)
    plt.yscale('log')
    plt.xscale('log')
  #  plt.title("Effective lifetime  vs Injection level ("+savelabel+")\n")
    plt.ylabel('Effective lifetime (ms)',fontsize = 19)
    plt.xlabel(r'Excess carrier density ($\mathregular{cm^{-3}}$)',fontsize = 18)

#    ax =plt.gca()
#    ticks=[]
#    for i in ax.get_yticks().tolist():
#        if int(i)>0:
#            ticks.append('{:,}'.format(int(i)))
#        else:
#            ticks.append('{:,}'.format(float(i)))
#    ax.set_yticklabels(ticks)
    #ax.tick_params(direction = 'in', which = 'both')
    plt.legend(loc='center left',bbox_to_anchor=(1,0.5), fontsize=14)
    plt.legend(fontsize=16,labelspacing=0.3)
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.ylim(1e-1,1e3)
#
    plt.xlim(1e14,1e17)
    f5.savefig("Effective ("+savelabel+").png", bbox_inches='tight', bbox_to_anchor=(1, 0.5), dpi = 300)
    
    f6 = plt.figure(6)
    if dash==False:
        plt.plot(dn_list,tau_surf_list,'-',label = plotlabel, color=color,linewidth=3)
    else:
        plt.plot(dn_list,tau_surf_list,'--',label = plotlabel, color=color,linewidth=3)
    plt.yscale('log')
    plt.xscale('log')
   # plt.title("Surf Lifetime vs Injection level ("+savelabel+")")
    plt.ylabel('Surface lifetime (ms)',fontsize = 19)
    plt.xlabel(r'Excess carrier density ($\mathregular{cm^{-3}}$)',fontsize = 18)

#    ax =plt.gca()
#    ticks=[]
#    for i in ax.get_yticks().tolist():
#        if int(i)>0:
#            ticks.append('{:,}'.format(int(i)))
#        else:
#            ticks.append('{:,}'.format(float(i)))
#    ax.set_yticklabels(ticks)
    #ax.tick_params(direction = 'in', which = 'both')
    plt.legend(fontsize=16,labelspacing=0.3, loc='upper right')
    plt.xticks(fontsize = 15)
    plt.yticks(fontsize = 15)
    plt.ylim(1e-1,1e3)
#
    plt.xlim(1e14,1e17)
    f6.savefig("Surf ("+savelabel+").png", bbox_inches='tight', bbox_to_anchor = (0,0), dpi = 300)
    

#####################################################################################################################################################################
def plot_contour(Rho_list,Y_data,Eff_no_surf,FF_no_surf,Voc_no_surf,Vmp_no_surf,Y_label, _log=True):
    f6 = plt.figure(6)    
    contour = plt.contourf(Rho_list,Y_data,Eff_no_surf, cmap = 'plasma_r', levels = np.linspace(np.min(Eff_no_surf), np.max(Eff_no_surf),100))
    plt.xscale('log')
    if _log:
        plt.yscale('log')
    plt.xticks(fontsize = 13)
    plt.yticks(fontsize = 13)
   # plt.title("Efficiency Vs Wafer resistivity \n (N-type Wafer) ",fontsize=14)
    plt.ylabel(Y_label,fontsize=14)
    plt.xlabel(r'Wafer resistivity ($\Omega$ cm)',fontsize=14)
    cbar = plt.colorbar(contour,format='%.0f')
    cbar.set_label('Efficiency (%)',fontsize=13, labelpad=+6)
    f6.savefig('Efficiency Vs Wafer resistivity - N-type.png', bbox_inches="tight")

    f7 = plt.figure(7)    
    contour2 = plt.contourf(Rho_list,Y_data,FF_no_surf, cmap = 'plasma_r',levels = np.linspace(np.min(FF_no_surf), np.max(FF_no_surf),100))
    plt.xscale('log')
    if _log:
        plt.yscale('log')
    plt.xticks(fontsize = 13)
    plt.yticks(fontsize = 13)
    
    #plt.title("FF Vs Wafer resistivity \n (N-type Wafer) ",fontsize=14)
    plt.ylabel(Y_label,fontsize=14)
    plt.ylabel(r'Bulk lifetime $\tau$$\mathregular{_{bulk}}$ (ms)',fontsize=14)
    plt.xlabel(r'Wafer resistivity ($\Omega$ cm)',fontsize=14)
    cbar2 = plt.colorbar(contour2,format='%.0f')
    cbar2.set_label('FF (%)',fontsize=13,labelpad=+6)
    f7.savefig('FF Vs Wafer resistivity - N-type.png', bbox_inches="tight")

    f8 = plt.figure(8)    
    contour3 = plt.contourf(Rho_list,Y_data,Voc_no_surf, cmap = 'plasma_r', levels = np.linspace(np.min(Voc_no_surf), np.max(Voc_no_surf),100))
    plt.xscale('log')
    if _log:
        plt.yscale('log')
    plt.xticks(fontsize = 13)
    plt.yticks(fontsize = 13)
  #  plt.title("Voc Vs Wafer resistivity \n (N-type Wafer) ",fontsize=14)
    plt.ylabel(Y_label,fontsize=14)
    plt.xlabel(r'Wafer resistivity ($\Omega$ cm)',fontsize=14)
    cbar3 = plt.colorbar(contour3,format='%.0f')
    cbar3.set_label('$\mathregular{V_{oc}}$ (mV)',fontsize=13, labelpad=+6)
    f8.savefig('Voc Vs Wafer resistivity - N-type.png', bbox_inches="tight")

    f9 = plt.figure(9)    
    contour4 = plt.contourf(Rho_list,Y_data,Vmp_no_surf, cmap = 'plasma_r', levels = np.linspace(np.min(Vmp_no_surf), np.max(Vmp_no_surf),100))
    plt.xscale('log')
    if _log:
        plt.yscale('log')
    plt.xticks(fontsize = 13)
    plt.yticks(fontsize = 13)
 #   plt.title("Vmp Vs Wafer resistivity \n (N-type Wafer) ",fontsize=14)
    plt.ylabel(Y_label,fontsize=14)
    plt.xlabel(r'Wafer resistivity ($\Omega$ cm)',fontsize=14)
    cbar4 = plt.colorbar(contour4,format='%.0f')
    cbar4.set_label('$\mathregular{V_{mp}}$ (mV)',fontsize=13, labelpad=+6)
    f9.savefig('Vmp Vs Wafer resistivity - N-type.png', bbox_inches="tight")
    
#######################################################################################################################################################
################## Lambertian light trapping   ###################

#computes Lambertian light absorption from wafer thickness
df_jl = pd.read_csv(r'lambertian light inputs.csv')
import scipy.integrate as it

W_list = np.linspace(20,400,39) #list of wafer thicknesses
alpha = df_jl['alpha cSi']
n = df_jl['n']
Jl_list=[]

for w in W_list:
    w = w*1e-6
    lambertian = (1-np.exp(-4*alpha*w)+np.exp(-4*alpha*w)*(1-1/n**2)*(1-np.exp(-4*alpha*w))/(1-np.exp(-4*alpha*w)+1/n**2*np.exp(-4*alpha*w)))
    lambertian_w_r = (1 - df_jl['primary R'])*lambertian #reflection is typical of SHJ cell
    jl = it.simps(lambertian*df_jl['AM1.5g(l2)']*df_jl['WL'])*2*0.000040328 #integrate lambertian absorpbtion with AM1.5 spectrum
    #This seems suspect, but it gives the same result as the EPFL loss equations spreadsheet
    Jl_list.append(jl)
    
#######  Calculate Jr for a given J0, and a given thickness #########
    
W_list = W_list #list of wafer thicknesses from above in cm
Jl_dict = dict(zip(W_list, Jl_list)) #Jl from above





############################################################################################



##################### MAIN - PLOTS ########################
     ################    Inputs #############

n02=1
Rs=0.9
W=180          ######### Thickness of the cell - 180 um    ############
kb=1.3807e-23     ###### in (J/K)
k = 8.61758e-5 #[eV/K]

Jl = Jl_dict[W]
width = W*1e-4 #(in cm)

J0=10e-15



#####data from paper######
fname='impurity.txt'
impurities=[]
with open(fname, newline = '') as games:                                                                                          
    filedata = csv.reader(games, delimiter=',')
    for data in filedata:
        impurities.append(data)



colormap = plt.cm.gist_ncar
colors = [colormap(i) for i in np.linspace(0, 0.95, 8)] #play with the linsapce to vary the color gradients
#plt.gca().set_color_cycle(colors)

    ###############################################################################
                               ############## CHOOSE  ###################
if 1:
    ##########  Sweep temperature ##########
    defect_num = 12  ### see defect list text file or defect list from paaper file
    N_defects=5e9   ##### impurity defect density in cm-3
    n0=1.549e15  ##### doping density in cm-3 FOR 3 OHM -CM
    
    #n0=2.203e11  ## 20k ohm-cm
    
    init_t = 30 # in Celcius
    fin_t = 80 #i n Celcius
    div=10  # interval of the temp needed for output
    
    
    
    ini_T=init_t+273.15  ######## lowest temp in the range in Kelvins K
    fin_T=fin_t+273.15######## highest temp in the range in Kelvins K
    T_div=(fin_T-ini_T)/(div)+1
    T_list= np.linspace(ini_T,fin_T,T_div)  #### list of temp in the range in Kelvins K
    
    
    J01_list=[]
    tau_surf_list=[]
    
    workbook = xlsxwriter.Workbook('data.xlsx')
    
    for i in range(len(T_list)):
        T = T_list[i]
        sig_n = eval(impurities[defect_num][3])
        sig_p = eval(impurities[defect_num][4])
        Eg_minus_E_defects = float(impurities[defect_num][2])
        defect = impurities[defect_num][1]
        #calc =Lifetime(sig_p[defect_num-1], sig_n[defect_num-1], Eg_minus_E_defects[defect_num-1],defect[defect_num-1], T_list[i], n0, N_defects, J0, Jl, Rs, sweep_dn=False)
        calc =Lifetime(sig_p, sig_n, Eg_minus_E_defects,defect, T_list[i], n0, N_defects, J0, Jl, Rs, sweep_dn=False)
        print('Calculating for '+'T ='"{:.0f}\u00b0C".format(T_list[i]-273.15))
        ind1=min(range(len(calc.dn)), key=lambda i: abs(calc.dn[i]-1e13))
        ind2=min(range(len(calc.dn)), key=lambda i: abs(calc.dn[i]-1e15))
        ind3=min(range(len(calc.dn)), key=lambda i: abs(calc.dn[i]-1e17))
        dn_list = [calc.dn[ind1], calc.dn[ind2],calc.dn[ind3]]
        tau_rad_List=[calc.tau_rad[ind1], calc.tau_rad[ind2],calc.tau_rad[ind3]]
        tau_srh_List=[calc.tau_srh[ind1], calc.tau_srh[ind2],calc.tau_srh[ind3]]
        tau_surf_List=[calc.tau_surf[ind1], calc.tau_surf[ind2],calc.tau_surf[ind3]]
        tau_eff_no_surf_List=[calc.tau_eff_no_surf[ind1], calc.tau_eff_no_surf[ind2],calc.tau_eff_no_surf[ind3]]
        tau_eff_List=[calc.tau_eff[ind1], calc.tau_eff[ind2],calc.tau_eff[ind3]]
        p0_List=[calc.p0[ind1], calc.p0[ind2],calc.p0[ind3]]
        n_List=[calc.n[ind1], calc.n[ind2],calc.n[ind3]]
        p_List=[calc.p[ind1], calc.p[ind2],calc.p[ind3]]
        tau_n0_List=[calc.tau_n0]
        tau_p0_List=[calc.tau_p0]
        n1_List=[calc.n1]
        p1_List=[calc.p1]
        J01_List=[calc.J01]
        ni_eff_List=[calc.ni_eff[ind1], calc.ni_eff[ind2],calc.ni_eff[ind3]]
       
#        sig_nT_list=[calc.sig_nT[ind1], calc.sig_nT[ind2],calc.sig_nT[ind3]]
#        sig_pT_list=[calc.sig_pT[ind1], calc.sig_pT[ind2],calc.sig_pT[ind3]]


        worksheet = workbook.add_worksheet(str(T_list[i]))
        worksheet.write(0,0,'dn')
        worksheet.write(0,1,'tau_rad')
        worksheet.write(0,2,'tau_srh')
        worksheet.write(0,3,'tau_surf')
        worksheet.write(0,4,'tau_eff_no_surf')
        worksheet.write(0,5,'tau_eff')
        worksheet.write(0,6,'tau_n0')
        worksheet.write(0,7,'tau_p0')
        worksheet.write(0,8,'n0')
        worksheet.write(0,9,'p0')
        worksheet.write(0,10,'n1')
        worksheet.write(0,11,'p1')
        worksheet.write(0,12,'n')
        worksheet.write(0,13,'p')
        worksheet.write(0,14,'J01')
        worksheet.write(0,15,'ni_eff')
        




        worksheet.write_column(1,0,dn_list)
        worksheet.write_column(1,1,tau_rad_List)
        worksheet.write_column(1,2,tau_srh_List)
        worksheet.write_column(1,3,tau_surf_List)
        worksheet.write_column(1,4,tau_eff_no_surf_List)
        worksheet.write_column(1,5,tau_eff_List)
        worksheet.write_column(1,6,tau_n0_List)
        worksheet.write_column(1,7,tau_p0_List)
#        worksheet.write_column(1,8,n0_List)
        worksheet.write_column(1,9,p0_List)
        worksheet.write_column(1,10,n1_List)
        worksheet.write_column(1,11,p1_List)
        worksheet.write_column(1,12,n_List)
        worksheet.write_column(1,13,p_List)
        worksheet.write_column(1,14,J01_List)
        worksheet.write_column(1,15,ni_eff_List)
        



        plot_lifetime(calc.tau_rad, calc.tau_aug, calc.tau_srh, calc.tau_eff_no_surf,calc.tau_eff,calc.tau_surf, calc.dn, 'T ='"{:.0f}\u00b0C".format(T_list[i]-273.15), defect, colors[i])
        dn=calc.dn
        B = calc.B
        Blow=calc.Blow
        ni_eff=calc.ni_eff
        dEg=calc.dEg
        J01_list.append(calc.J01)
        tau_surf_list.append(calc.tau_surf)
        
    workbook.close()
###############################################################################
                               ############## CHOOSE  ###################
########################## sweep all defects     ####################
if 0:
    
    N_defects=1e9  ##### impurity defect density in cm-3
    #n0=1.549e15
    n0=2.203e11##### doping density in cm-3
    t = 25     ### temp in C    
    defects_list=[1,2,6,7,14,10,11]
    
    
    T =t+273.15 
    i=0
    for defect_num in defects_list:
        calc = Lifetime(eval(impurities[defect_num][3]), eval(impurities[defect_num][4]), float(impurities[defect_num][2]),impurities[defect_num][1], T, n0, N_defects, J0, Jl, Rs, sweep_dn=False)
        print('Calculating for '+impurities[defect_num][1])
        plot_lifetime(calc.tau_rad, calc.tau_aug, calc.tau_srh, calc.tau_eff_no_surf, calc.tau_eff,calc.tau_surf, calc.dn, impurities[defect_num][1], 'T ='"{:.0f}\u00b0C".format(T-273.15),colors[i])
        i+=1
##################################################################################

                               ############## CHOOSE  ###################
##########################  sweep N_defect - number of defects or defect concentration     ####################
if 0:
    t = 25
    defect_num = 10### see defect list text file
      ##### doping density in cm-3
    n0=1.549e15
    N_defect_list = np.logspace(8,10,3)

    T = t+273.15
    colors = ['limegreen','orangered','gold']
    f = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    g = lambda x,pos : "${}$".format(f._formatSciNotation('%1.10e' % x))
    fmt = mticker.FuncFormatter(g)
    
    Eff_no_surf_result=[]    
    Eff_no_surf_result2=[]
    
    Eff_result=[]    
    Eff_result2=[]
 
    Eff_intrinsic_result=[]    
    Eff_intrinsic_result2=[]
    
    for i in range(len(N_defect_list)):
        calc = Lifetime(eval(impurities[defect_num][3]), eval(impurities[defect_num][4]), float(impurities[defect_num][2]),impurities[defect_num][1], T, n0, N_defect_list[i], J0, Jl, Rs, sweep_dn=False)
        print('Calculating for number of defects = '+"{:.0E}".format(Decimal(N_defect_list[i])))
        plot_lifetime(calc.tau_rad, calc.tau_aug, calc.tau_srh, calc.tau_eff_no_surf,calc.tau_eff, calc.tau_surf,calc.dn, r"{} ".format(fmt(N_defect_list[i])), impurities[defect_num][1]+' at '+'T ='"{:.0f}\u00b0C".format(T-273.15),colors[i])
        Eff_no_surf_result.append(calc.Eff_no_surf)
        Eff_result.append(calc.Eff)
        Eff_intrinsic_result.append(calc.Eff_intrinsic)
        
    n0=2.203e11
    
    for i in range(len(N_defect_list)):
        calc = Lifetime(eval(impurities[defect_num][3]), eval(impurities[defect_num][4]), float(impurities[defect_num][2]),impurities[defect_num][1], T, n0, N_defect_list[i], J0, Jl, Rs, sweep_dn=False)
        print('Calculating for number of defects = '+"{:.0E}".format(Decimal(N_defect_list[i])))
        plot_lifetime(calc.tau_rad, calc.tau_aug, calc.tau_srh, calc.tau_eff_no_surf,calc.tau_eff, calc.tau_surf,calc.dn, r"{}".format(fmt(N_defect_list[i])), impurities[defect_num][1]+' at '+'T ='"{:.0f}\u00b0C".format(T-273.15),colors[i],dash=True)
        Eff_no_surf_result2.append(calc.Eff_no_surf)
        Eff_result2.append(calc.Eff)
        Eff_intrinsic_result2.append(calc.Eff_intrinsic)
################################################################################


                               ############## CHOOSE  ###################
# ######################### sweep n0 - doping concentration   ####################
if 0:
    t = 25  #### temp in C
    defect_num = 3  ### see defect list text file
    N_defects=1e9  ##### impurity defect density in cm-3
    
    
    n0_list = np.logspace(11,16,6)
    
    T = t+273.15
    
    for i in range(len(n0_list)):
        calc = Lifetime(sig_p[defect_num-1], sig_n[defect_num-1], Eg_minus_E_defects[defect_num-1],defect[defect_num-1], T, n0_list[i], N_defects, J0, Jl, Rs, sweep_dn=True)
        print('Calculating for doping concentration = '+"{:.0E}".format(Decimal(n0_list[i])))
        plot_lifetime(calc.tau_rad, calc.tau_aug, calc.tau_srh, calc.tau_eff_no_surf, calc.dn, "{:.0E}".format(Decimal(n0_list[i])), defect[defect_num-1].decode()+' at '+'T ='"{:.0f}\u00b0C".format(T-273.15),colors[i])

#######################################################
 
    ###################### PLot tau_bulk vs doping and vs resistivity ################
if 0:
    
    t = 25  #### temp in C
    defect_num = 10  ### see defect list text file
    N_defects=1e9  ##### impurity defect density in cm-3
    
    doping_input = pd.read_csv(r'ntype.csv')
    n0_list = doping_input['Dopant concentration']
    Rho_list =doping_input['Resistivity']
    Rho_dict = dict(zip(n0_list, Rho_list))

#    n0_list = np.logspace(11,16,20)
    
    T = t+273.15
    
    bulk_lifetime=[]
    lifetime=[]
    SRH=[]
    RAD=[]
    AUG=[]
    surf=[]
    Voc=[]
    Vmp=[]
    Eff=[]
    FF=[]
    Voc_no_surf=[]
    Vmp_no_surf=[]
    Eff_no_surf=[]
    FF_no_surf=[]
    
    for i in range(len(n0_list)):
 #       Rho = Rho_dict[Nd(i)]
        calc = Lifetime(eval(impurities[defect_num][3]), eval(impurities[defect_num][4]), float(impurities[defect_num][2]),impurities[defect_num][1], T, n0_list[i], N_defects, J0, Jl, Rs, sweep_dn=False)
        #calc = Lifetime(sig_p[defect_num-1], sig_n[defect_num-1], Eg_minus_E_defects[defect_num-1],defect[defect_num-1], T, n0_list[i], N_defects, J0, Jl, Rs, sweep_dn=False)
        dn = list(calc.dn)
        inj= min(dn, key=lambda x:abs(x-1e15))
        ind = dn.index(inj)
        bulk_lifetime.append(calc.tau_eff_no_surf[ind])
        lifetime.append(calc.tau_eff[ind])
        SRH.append(calc.tau_srh[ind])
        RAD.append(calc.tau_rad[ind])
        AUG.append(calc.tau_aug[ind])
        surf.append(calc.tau_surf[ind])
        Voc.append(calc.Voc)
        Vmp.append(calc.Vmp)
        FF.append(calc.FF)
        Eff.append(calc.Eff)
        Voc_no_surf.append(calc.Voc_no_surf)
        Vmp_no_surf.append(calc.Vmp_no_surf)
        Eff_no_surf.append(calc.Eff_no_surf)
        FF_no_surf.append(calc.FF_no_surf)

    f6 = plt.figure(6)    
    plt.plot(n0_list, bulk_lifetime)
    plt.xscale('log')
    plt.xlabel('Doping concentration ($\mathregular{cm^{-3}}$)',fontsize=14)
    plt.ylabel('Bulk lifetime (ms)',fontsize=14.5)
    plt.xticks(fontsize = 13)
    plt.yticks(fontsize = 13)
    ax =plt.gca()
    #ax.set_xticklabels(['{:,}'.format(int(x)) for x in ax.get_xticks().tolist()])
    ticks=[]
    for i in ax.get_xticks().tolist():
        if int(i)>0:
            ticks.append('{:,}'.format(int(i)))
        else:
            ticks.append('{:,}'.format(float(i)))
    ax.set_xticklabels(ticks)
    ax.tick_params(direction = 'in', which = 'both')
    f6.savefig('Lifetime vs doping conc',bbox_inches='tight')


    f7 = plt.figure(7)
    plt.plot(Rho_list, SRH,label='Bulk SRH',linestyle='-.',color='g',linewidth=3)
    plt.plot(Rho_list, AUG,label='Auger',linestyle=':',color='m',linewidth=3)
    plt.plot(Rho_list, RAD,label='Radiative',linestyle='--',color='b',linewidth=3 )
    plt.plot(Rho_list, bulk_lifetime,label='Effective',linestyle='-',linewidth=4, color='r')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='lower right', fontsize=16)
    ax=plt.gca()
    ticks=[]
    for i in ax.get_xticks().tolist():
        if int(i)>0:
            ticks.append('{:,}'.format(int(i)))
        else:
            ticks.append('{:,}'.format(float(i)))
    ax.set_xticklabels(ticks)
    plt.xlabel(r'Wafer resistivity ($\Omega$cm)',fontsize=20)
    plt.ylabel('Lifetime (ms)',fontsize=23)
    
    #plt.ylabel(r'{\fontsize{29pt}{3em}\selectfont{}{Lifetime (ms)\n}{\fontsize{18pt}{3em}\selectfont{}(At 1E15 $\mathregular{cm^{-3}}$ MCD)}')
    plt.yticks(fontsize = 18)
    plt.xticks(fontsize = 18)
    plt.xlim(5e-2,5e3)
    plt.ylim(1e-2,1e3)
    ax =plt.gca()
    ticks=[]
    for i in ax.get_xticks().tolist():
        if int(i)>0:
            ticks.append('{:,}'.format(int(i)))
        else:
            ticks.append('{:,}'.format(float(i)))
    ax.set_xticklabels(ticks)
    ticks=[]
    for i in ax.get_yticks().tolist():
        if int(i)>0:
            ticks.append('{:,}'.format(int(i)))
        else:
            ticks.append('{:,}'.format(float(i)))
    ax.set_yticklabels(ticks)
    #ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])
    ax.tick_params(direction = 'in', which = 'both')
    f7.savefig('Lifetime (no surf)vs resistivity',bbox_inches='tight')
    
    
    f8 = plt.figure(8)
    plt.plot(Rho_list, SRH,label='Bulk SRH',linestyle='-.',color='g',linewidth=3)
    plt.plot(Rho_list, AUG,label='Auger',linestyle=':',color='m',linewidth=3)
    plt.plot(Rho_list, RAD,label='Radiative',linestyle='--',color='b',linewidth=3 )
    plt.plot(Rho_list, surf,label='Surface',linestyle='-',color='c',linewidth=3 )
   
    plt.plot(Rho_list, lifetime,label='Effective',linestyle='-',linewidth=4, color='r')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='lower right', fontsize=16)
    plt.xlabel(r'Wafer resistivity ($\Omega$cm)',fontsize=20)
    plt.ylabel('Lifetime (ms)',fontsize=23)
    
    #plt.ylabel(r'{\fontsize{29pt}{3em}\selectfont{}{Lifetime (ms)\n}{\fontsize{18pt}{3em}\selectfont{}(At 1E15 $\mathregular{cm^{-3}}$ MCD)}')
    plt.yticks(fontsize = 18)
    plt.xticks(fontsize = 18)
    plt.xlim(5e-2,5e3)
    plt.ylim(1e-2,1e3)
    ax =plt.gca()
    ticks=[]
    for i in ax.get_xticks().tolist():
        if int(i)>0:
            ticks.append('{:,}'.format(int(i)))
        else:
            ticks.append('{:,}'.format(float(i)))
    ax.set_xticklabels(ticks)
    ticks=[]
    for i in ax.get_yticks().tolist():
        if int(i)>0:
            ticks.append('{:,}'.format(int(i)))
        else:
            ticks.append('{:,}'.format(float(i)))
    ax.set_yticklabels(ticks)
    #ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])
    ax.tick_params(direction = 'in', which = 'both')
    f8.savefig('Lifetime(surf) vs resistivity',bbox_inches='tight')
    
    f9 = plt.figure(9)
#    plt.plot(Rho_list, Voc_no_surf,label='Effective',linestyle='-',linewidth=4, color='b')
    plt.plot(Rho_list, Voc,label='$\mathregular{iV_{OC}}$',linestyle='-',linewidth=4, color='r')
    plt.xscale('log')
    #plt.yscale('log')
    #plt.legend(loc='lower right', fontsize=16)
    plt.xlabel(r'Bulk resistivity ($\Omega$cm)',fontsize=20)
    plt.ylabel('$\mathregular{iV_{OC}}$ (mV)',fontsize=23)
    
    #plt.ylabel(r'{\fontsize{29pt}{3em}\selectfont{}{Lifetime (ms)\n}{\fontsize{18pt}{3em}\selectfont{}(At 1E15 $\mathregular{cm^{-3}}$ MCD)}')
    plt.yticks(fontsize = 18)
    plt.xticks(fontsize = 18)
    plt.xlim(5e-2,5e3)
    plt.ylim(660,760)
    ax =plt.gca()
    ticks=[]
    for i in ax.get_xticks().tolist():
        if int(i)>0:
            ticks.append('{:,}'.format(int(i)))
        else:
            ticks.append('{:,}'.format(float(i)))
    ax.set_xticklabels(ticks)
    ticks=[]
    for i in ax.get_yticks().tolist():
        if int(i)>0:
            ticks.append('{:,}'.format(int(i)))
        else:
            ticks.append('{:,}'.format(float(i)))
    ax.set_yticklabels(ticks)
    #ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])
    ax.tick_params(direction = 'in', which = 'both')
    f9.savefig('Voc vs resistivity',bbox_inches='tight')       
    
    f10 = plt.figure(10)
#    plt.plot(Rho_list, Vmp_no_surf,label='Effective',linestyle='-',linewidth=4, color='b')
    plt.plot(Rho_list, Vmp,label='$\mathregular{iV_{mp}}$',linestyle='-',linewidth=4, color='r')
    plt.xscale('log')
    #plt.yscale('log')
    #plt.legend(loc='lower right', fontsize=16)
    plt.xlabel(r'Bulk resistivity ($\Omega$cm)',fontsize=20)
    plt.ylabel('$\mathregular{iV_{MP}}$ (mV)',fontsize=23)
    
    #plt.ylabel(r'{\fontsize{29pt}{3em}\selectfont{}{Lifetime (ms)\n}{\fontsize{18pt}{3em}\selectfont{}(At 1E15 $\mathregular{cm^{-3}}$ MCD)}')
    plt.yticks(fontsize = 18)
    plt.xticks(fontsize = 18)
    plt.xlim(5e-2,5e3)
    plt.ylim(580,680)
    ax =plt.gca()
    ticks=[]
    for i in ax.get_xticks().tolist():
        if int(i)>0:
            ticks.append('{:,}'.format(int(i)))
        else:
            ticks.append('{:,}'.format(float(i)))
    ax.set_xticklabels(ticks)
    ticks=[]
    for i in ax.get_yticks().tolist():
        if int(i)>0:
            ticks.append('{:,}'.format(int(i)))
        else:
            ticks.append('{:,}'.format(float(i)))
    ax.set_yticklabels(ticks)    #ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])
    ax.tick_params(direction = 'in', which = 'both')
    f10.savefig('Vmp vs resistivity',bbox_inches='tight')       
    
    f11 = plt.figure(11)
#    plt.plot(Rho_list, FF_no_surf,label='Effective',linestyle='-',linewidth=4, color='b')
    plt.plot(Rho_list, FF,label='iFF',linestyle='-',linewidth=4, color='r')
    plt.xscale('log')
    #plt.yscale('log')
    #plt.legend(loc='lower right', fontsize=16)
    plt.xlabel(r'Bulk resistivity ($\Omega$cm)',fontsize=20)
    plt.ylabel('iFF (%)',fontsize=23)
    
    #plt.ylabel(r'{\fontsize{29pt}{3em}\selectfont{}{Lifetime (ms)\n}{\fontsize{18pt}{3em}\selectfont{}(At 1E15 $\mathregular{cm^{-3}}$ MCD)}')
    plt.yticks(fontsize = 18)
    plt.xticks(fontsize = 18)
    plt.xlim(5e-2,5e3)
    plt.ylim(83,88)
    ax =plt.gca()
    ticks=[]
    for i in ax.get_xticks().tolist():
        if int(i)>0:
            ticks.append('{:,}'.format(int(i)))
        else:
            ticks.append('{:,}'.format(float(i)))
    ax.set_xticklabels(ticks)
    ticks=[]
    for i in ax.get_yticks().tolist():
        if int(i)>0:
            ticks.append('{:,}'.format(int(i)))
        else:
            ticks.append('{:,}'.format(float(i)))
    ax.set_yticklabels(ticks)    #ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])
    ax.tick_params(direction = 'in', which = 'both')
    f11.savefig('FF vs resistivity',bbox_inches='tight')       
    
    f12 = plt.figure(12)
#    plt.plot(Rho_list, FF_no_surf,label='Effective',linestyle='-',linewidth=4, color='b')
    plt.plot(Rho_list, Eff,label='iEfficiency',linestyle='-',linewidth=4, color='r')
    plt.xscale('log')
    #plt.yscale('log')
    #plt.legend(loc='lower right', fontsize=16)
    plt.xlabel(r'Bulk resistivity ($\Omega$cm)',fontsize=20)
    plt.ylabel('iEfficiency (%)',fontsize=23)
    
    #plt.ylabel(r'{\fontsize{29pt}{3em}\selectfont{}{Lifetime (ms)\n}{\fontsize{18pt}{3em}\selectfont{}(At 1E15 $\mathregular{cm^{-3}}$ MCD)}')
    plt.yticks(fontsize = 18)
    plt.xticks(fontsize = 18)
    plt.xlim(5e-2,5e3)
    plt.ylim(24,30)
    ax =plt.gca()
    ticks=[]
    for i in ax.get_xticks().tolist():
        if int(i)>0:
            ticks.append('{:,}'.format(int(i)))
        else:
            ticks.append('{:,}'.format(float(i)))
    ax.set_xticklabels(ticks)
    ticks=[]
    for i in ax.get_yticks().tolist():
        if int(i)>0:
            ticks.append('{:,}'.format(int(i)))
        else:
            ticks.append('{:,}'.format(float(i)))
    ax.set_yticklabels(ticks)    #ax.set_yticklabels(['{:,}'.format(int(x)) for x in ax.get_yticks().tolist()])
    ax.tick_params(direction = 'in', which = 'both')
    f12.savefig('Eff vs resistivity',bbox_inches='tight') 
######################33











if 0:
##################################################### CHOOSE - CONTOURS PLOTS  ###############################33
  
    ################ Sweep defects concentration , constant T ###############
    defect_num = 3  ### see defect list text file
    t=25  #### temp in C


    
    doping_input = pd.read_csv(r'ntype.csv')
    n0_list = doping_input['Dopant concentration']
    Rho_list =doping_input['Resistivity']
    
    Rho_dict = dict(zip(n0_list, Rho_list))
    
    T =t+273.15 
    N_defect_list=list(np.logspace(6,11,5))
    
    
    
    Voc_no_surf_counter=[]
    Vmp_no_surf_counter=[]
    Eff_no_surf_counter=[]
    FF_no_surf_counter=[]
    dn_mpp_no_surf_counter=[]
    dn_oc_no_surf_counter=[]
    
    for N_defect in N_defect_list:
    
        Vmp_intermediate_no_surf_counter=[]
        Voc_intermediate_no_surf_counter=[]
        Eff_intermediate_no_surf_counter=[]
        FF_intermediate_no_surf_counter=[]
        Jr_intermediate_no_surf_counter=[]
        
        for n0 in n0_list:
            #B,tau_rad, tau_aug, tau_srh, tau_eff_no_surf, dn,Voc_no_surf, Vmp_no_surf,Eff_no_surf, FF_no_surf,dn_mpp_no_surf,dn_oc_no_surf,tau_n0,tau_p0 = calculation(sig_p[defect_num-1], sig_n[defect_num-1], Eg_minus_E_defects[defect_num-1],defect[defect_num-1], T, n0, N_defect,J0,Jl,Rs )
            calc = Lifetime(sig_p[defect_num-1], sig_n[defect_num-1], Eg_minus_E_defects[defect_num-1],defect[defect_num-1], T, n0, N_defect,J0,Jl,Rs)
            #print('tau_n0={}  tau_p0 ={}'.format(calc.tau_n0*1000,calc.tau_p0*1000))
            Vmp_intermediate_no_surf_counter.extend(calc.Vmp_no_surf)
            Voc_intermediate_no_surf_counter.extend(calc.Voc_no_surf)
            Eff_intermediate_no_surf_counter.append(calc.Eff_no_surf)
            FF_intermediate_no_surf_counter.extend(calc.FF_no_surf)
            
            
        Vmp_no_surf_counter.append(Vmp_intermediate_no_surf_counter)
        Voc_no_surf_counter.append(Voc_intermediate_no_surf_counter)
        Eff_no_surf_counter.append(Eff_intermediate_no_surf_counter)
        FF_no_surf_counter.append(FF_intermediate_no_surf_counter)
            
    Y_data=N_defect_list
    Y_label=r'Concentration of defects (cm-3)'
    plot_contour(Rho_list,Y_data,Eff_no_surf_counter,FF_no_surf_counter,Voc_no_surf_counter,Vmp_no_surf_counter,Y_label)


if 0:
##################################################### CHOOSE - CONTOURS  ###############################33

    ################ Sweep T, constant N defects ########
    defect_num = 3  ### see defect list text file
    N_defects=1e9  ##### impurity defect density in cm-3

    
    doping_input = pd.read_csv(r'ntype.csv')
    n0_list = doping_input['Dopant concentration']
    Rho_list =doping_input['Resistivity']
    
    Rho_dict = dict(zip(n0_list, Rho_list))
    
    T_list=list(np.linspace(0,40,10))
    
    
    
    Voc_no_surf_counter=[]
    Vmp_no_surf_counter=[]
    Eff_no_surf_counter=[]
    FF_no_surf_counter=[]
    dn_mpp_no_surf_counter=[]
    dn_oc_no_surf_counter=[]
    
    for T in T_list:
    
        Vmp_intermediate_no_surf_counter=[]
        Voc_intermediate_no_surf_counter=[]
        Eff_intermediate_no_surf_counter=[]
        FF_intermediate_no_surf_counter=[]
        Jr_intermediate_no_surf_counter=[]
        
        for n0 in n0_list:
            calc = Lifetime(sig_p[defect_num-1], sig_n[defect_num-1], Eg_minus_E_defects[defect_num-1],defect[defect_num-1], T+273.15, n0, N_defects,J0,Jl,Rs)
            #print('tau_n0={}  tau_p0 ={}'.format(calc.tau_n0*1000,calc.tau_p0*1000))
            Vmp_intermediate_no_surf_counter.extend(calc.Vmp_no_surf)
            Voc_intermediate_no_surf_counter.extend(calc.Voc_no_surf)
            Eff_intermediate_no_surf_counter.append(calc.Eff_no_surf)
            FF_intermediate_no_surf_counter.extend(calc.FF_no_surf)
            
        Vmp_no_surf_counter.append(Vmp_intermediate_no_surf_counter)
        Voc_no_surf_counter.append(Voc_intermediate_no_surf_counter)
        Eff_no_surf_counter.append(Eff_intermediate_no_surf_counter)
        FF_no_surf_counter.append(FF_intermediate_no_surf_counter)
            
    Y_data=T_list
    Y_label=r'Temperature (°C)'
    plot_contour(Rho_list,Y_data,Eff_no_surf_counter,FF_no_surf_counter,Voc_no_surf_counter,Vmp_no_surf_counter,Y_label, _log = False)



