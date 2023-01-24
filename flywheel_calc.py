import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import signal

import math
import time

#flywheel_car omdat dit model oorspronkelijk bedoeld was voor het vliegwiel model
def flywheel_car(k,u,r_ig,r_og,g_w,r_a,d_g,m_t,t_in,t_end,dt,nrg):

    #wiel eigenschappen
    r_ws = 0.036
    r_wb = 0.055
    m_ws = 0.0625
    m_wb = 0.165
    
    #wrijving's coefficienten
    c_w = 0.002
    c_b = 0.001
    c_g = 0.04
    c_st = 0.7

    #tandwiel ratio
    t_r = (r_og/r_ig)**-1

    max_dist = u/(2*math.pi*r_a)*t_r*2*math.pi*r_wb*2
    slipped = False

    n = int((t_end-t_in)/dt)
    t = np.zeros(n+1)
    x = np.zeros(n+1)
    e = np.zeros(n+1)

    m_ig = d_g*g_w*math.pi*r_ig**2
    m_og = d_g*g_w*math.pi*r_og**2
    
    I_ws = 0.5*m_ws*r_ws**2
    I_wb = 0.5*m_wb*r_wb**2
    I_ig = 0.5*m_ig*r_ig**2
    I_og = 0.5*m_og*r_og**2
    u_i = u
    x_i = 0
    o_fw_i =0
    o_w_i =0.0
    v_i = 0
    x_total =0
    t_i = t_in
    t[0] = t_i
    x[0] = x_i
    e[0] = 0.5*k*u**2
    
    #constanten
    t_wf = (c_w*m_t*9.81/4.0)*(2.0*r_wb+2.0*r_ws)

    div  = I_ig +(2.0*I_wb+(2.0*I_ws*(r_wb/r_ws)**2)+m_t*r_wb**2+I_og)*t_r**2 
    div = div**-1
    
    cosa =math.cos(20/180*math.pi)
    cosa = cosa**-1
    tana = math.tan(20/180*math.pi)
    dtinv = dt**-1
    for i in range(1,n+1):
        fr_dir = math.atan((t_r*o_w_i)/(r_wb))
        #moment van de veer
        t_sp = -k*u_i*r_a#
        
        #krachten van de tandwielen opelkaar en op de as
        ft = (t_sp/r_ig)*cosa
        fr = ft*tana
        fn = math.sqrt(fr**2 + ft**2)

        #lagerwrijving
        t_bf = c_b*r_a*((9.81*m_ig+abs(k*u_i))+((fn+9.81)*(m_t-2.0*m_ws-2.0*m_wb)) + 9.81*(m_t-2.0*m_ws-2.0*m_wb))
        #tandwiel wrijving
        t_gf = fn*c_g*(r_ig+r_og)
        
        t_total = t_sp-((t_bf+t_wf+t_gf)*fr_dir)
        alpha = t_total*div

        o_fw_i  +=alpha*dt
        o_w_i  += alpha*t_r*dt
        e1 = 0.5*k*u_i**2 
        u_i += o_fw_i*r_a*dt
        e2 = 0.5*k*u_i**2
        x_i += o_w_i*r_wb*dt

        p = (e2-e1)*dtinv
        if(p/o_w_i > (c_st*9.81*m_t*0.25*2.0*r_wb)):
            slipped = True

        if(nrg == True):
        #energie
            e[i]  = 0.5*k*u_i**2
            e[i] += 0.5*2*I_wb*(o_w_i**2)
            e[i] += 0.5*2*I_ws*(o_w_i*(r_wb/r_ws))**2
            e[i] += 0.5*m_t*((o_w_i*r_wb)**2)
            e[i] += 0.5*I_ig*((o_fw_i)**2)
            e[i] += 0.5*I_og*((o_w_i)**2)
        t_i     += dt
        t[i]     = t_i
        x[i]     = x_i
    return t,x,e,max_dist,slipped


#het simpele model met alleen een veer en een wiek
def flywheel_spring(k,u,c,r_a,r_w,m,t_in,t_end,dt):
    n = int((t_end-t_in)/dt)
    t = np.zeros(n+1)
    x = np.zeros(n+1)
    omega = np.zeros(n+1)
    # -k*u(t)*r_a -c*omega= I_w*alfa 
    I = 0.5*(m*r_w**2)
    x_i = u
    o_i = 0
    ti = t_in
    
    t[0] = ti
    x[0] = x_i
    omega[0] = o_i
    for i in range(1,n+1):
        alpha  = ((-k*x_i)*r_a)/I
        o_i += alpha*dt
        x_i += o_i*r_a*dt
        ti+=dt 
        t[i] = ti
        omega[i] = o_i
        x[i] = abs(x_i)
    return t,x,omega

def find_peaks(fp):
    arr = fp
    r1 = min(arr)
    s = np.argmin(arr)
    arr = arr.tolist()
    r2 =max(arr[s:])
    return abs(r1),abs(r2)

def flywheel_opt(k,u,min_ra,max_ra,min_rig,max_rig,min_rog,max_rog,min_mt,max_mt,min_frm,st_ra,st_rig,st_rog,st_mt,dt,nrg, opt):
    #tandwiel eigenschappenx
    w_g = 0.008
    d_g = 1420
    
    #wiel eigenschappen
    m_wb = 0.165#165
    m_ws = 0.0625#62.5
    
    m_ig = 0
    m_og = 0
    x_t =0
    x_dif=0
    eff = 0
    x_s1 =0.000000001
    x_s2 =1000000000
    
    max_d = 8.75
    ddist = 0
    
    t_in =0
    t_end = 20
    ra_res = 0.0
    rig_res = 0.0
    rog_res = 0.0
    mt_res = 0.0
    xtotal = 0.0
    

    
    s1_a = np.zeros(int((max_ra-min_ra)/st_ra*(max_rig-min_rig)/st_rig*(max_rog-min_rog)/st_rog*(max_mt-min_mt)/st_mt))
    s2_a = np.zeros(int((max_ra-min_ra)/st_ra*(max_rig-min_rig)/st_rig*(max_rog-min_rog)/st_rog*(max_mt-min_mt)/st_mt))
    eff_a = np.zeros(int((max_ra-min_ra)/st_ra*(max_rig-min_rig)/st_rig*(max_rog-min_rog)/st_rog*(max_mt-min_mt)/st_mt))
    
    st = time.perf_counter()
    for i in range(0,100):
        flywheel_car(346,0.251,0.03,0.01,0.008,0.005,2700,1,t_in,t_end,dt,nrg)
    end = time.perf_counter()
    est = end-st
    est /=50
    st = time.perf_counter()
    print(str((max_ra-min_ra)/st_ra*(max_rig-min_rig)/st_rig*(max_rog-min_rog)/st_rog*(max_mt-min_mt)/st_mt*est/60) + " min. verwachte rekentijd")
    #itereer over alle waarden
    arrc =0;
    for r_a in np.arange(min_ra,max_ra,st_ra):
        for r_ig in np.arange(min_rig,max_rig,st_rig):
            for r_og in np.arange(min_rog, max_rog, st_rog):
                for m_t in np.arange(min_mt, max_mt, st_mt):
                    m_ig = d_g*w_g*math.pi*r_ig**2
                    m_og = d_g*w_g*math.pi*r_og**2
                    m_s = float(m_t)-(2*m_wb + 2*m_ws +m_ig +m_og+ min_frm)
                    if m_s > 0.0:
                        #vind de piekken en vergelijk ze
                        t,x,e,md,sl=flywheel_car(k,u,r_ig,r_og,w_g,r_a,d_g,m_t,t_in,t_end,dt,nrg)
                        tx_s1 , tx_s2 = find_peaks(x)
                        tx_s2 = (tx_s1-tx_s2)
                        tx_t = tx_s1+tx_s2
                        t_eff = tx_t/(2*md)
                        if(t_eff > 1.0):
                            t_eff = 1/t_eff

                        if(sl == False and opt == False):
                            s1_a = np.append(s1_a,tx_s1)
                            s2_1 = np.append(s2_a,tx_s2)
                            eff_a = np.append(eff_a,t_eff)
                        if(tx_t > x_t and tx_s2/tx_s1 < x_s2/x_s1 and tx_s1 <= max_d and sl == False and opt == True):
                            ra_res = r_a
                            rig_res = r_ig
                            rog_res = r_og
                            mt_res = m_t
                            x_t = tx_t
                            x_s1 = tx_s1
                            x_s2 = tx_s2
                            ddist = md
                            eff = t_eff
                            s1_a = np.append(s1_a,tx_s1)
                            s2_1 = np.append(s2_a,tx_s2)
                            eff_a = np.append(eff_a,t_eff)

                        else:
                            t_in =0
    
    end = time.perf_counter()
    print("tijd: " + str(end-st)+ "\n")
    print("totaal: "+ str(x_t))
    print("s1: "+ str(x_s1))
    print("s2: "+ str(x_s2))
    print("rendement:"+str(eff))
    print("radius as: "+str(ra_res))
    print("radius 1ste tandwiel: "+ str(rig_res))
    print("radius 2de tandwiel: " +str(rog_res))
    print("tandwiel ratio: " + str((rog_res/rig_res)**-1))
    print("totale massa: " +str(mt_res))
    return x_t,x_s1,x_s2,ra_res,rig_res,rog_res,eff,mt_res,s1_a,s2_a,eff_a


xt,s1,s2,ra,rig,rog,mt,ef,s1a,s2a,effa =flywheel_opt(346,0.251,0.004,0.005,0.01,0.04,0.01,0.04,1.0,2.0,0.5,0.0005,0.0005,0.0005,0.05,1,True,True)
x,y,z,u,v = flywheel_car(346,0.251,rig,rog,0.008,ra,1420,mt,0,100,0.0001,True)

fig,ax = plt.subplots(2,1)

ax[0].set_xlabel("tijd in seconden",fontsize = 18)
ax[0].set_ylabel("verplaatsing in meter",fontsize =18)
ax[1].set_xlabel("tijd in seconden",fontsize = 18)
ax[1].set_ylabel("Energie in Joule", fontsize = 18)

ax[0].spines['left'].set_linewidth(2.5)
ax[0].spines['bottom'].set_linewidth(2.5)
ax[1].spines['left'].set_linewidth(2.5)
ax[1].spines['bottom'].set_linewidth(2.5)

ax[0].tick_params(axis='both', which='major', labelsize=18, width=2.5, length=5)
ax[1].tick_params(axis='both', which='major', labelsize=18, width=2.5, length=5)

ax[0].grid()
ax[1].grid()

ax[0].plot(x, y, linewidth = 4)
ax[1].plot(x,z,linewidth =4, color ="red")

print(u)
print(v)
plt.show()


