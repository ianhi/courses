import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation
import numexpr as ne
mpl.style.use(['fivethirtyeight', 'figsize'])


# define constants
lx =128
ly = lx
Nsteps = 2000
ex=[0,1,0,-1,0 ,1,-1,-1, 1]
ey=[0,0,1,0 ,-1,1, 1,-1,-1]
w = [4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36]
G = -5
omega = 1
rho_psi = 1

cs = 1/np.sqrt(3)
cs2  = 1.00 / 3.00
cs22 = 2.00 * cs2
cssq = 2.0 / 9.00
#plotting functions:
def flip_T(a):
    return np.flipud(a.T)
def calc_macro_f(f):
    #Assumes f at a given time step
    #  this fxn shouldn't be used in the loops as it takes times to create new numpy arrays
    # useful after running for getting macro quantities w/o having the rho,u, and v arrays handy
    u = np.zeros([lx,ly])
    v = np.zeros([lx,ly])

    rho = np.sum(f,axis=0)
    for k in range(9):
        u += f[k]*ex[k]
        v += f[k]*ey[k]
    u=u/rho
    v=v/rho
    return rho, u ,v

# simulation functions
def equilibrium(rho,u,v):
    feq = np.zeros([9,lx,ly])
    usq = u*u
    vsq = v*v
    sumsq = (usq+vsq)/cs22
    sumsq2 = sumsq*(1-cs2)/cs2
    u22 = usq/cssq
    v22 = vsq/cssq
    ui = u/cs2
    vi = v/cs2
    uv = ui*vi
    feq[0] = (4.0/9.0)*rho*(1.00 - sumsq)

    feq[1] = (1.0/9.0)*rho*(1.00 - sumsq + u22 + ui)
    feq[2] = (1.0/9.0)*rho*(1.00 - sumsq + v22 + vi)
    feq[3] = (1.0/9.0)*rho*(1.00 - sumsq + u22 - ui)
    feq[4] = (1.0/9.0)*rho*(1.00 - sumsq + v22 - vi)

    feq[5] = (1.0/36.0)*rho*(1.00 + sumsq2 + ui + vi + uv)
    feq[6] = (1.0/36.0)*rho*(1.00 + sumsq2 - ui + vi - uv)
    feq[7] = (1.0/36.0)*rho*(1.00 + sumsq2 - ui - vi + uv)
    feq[8] = (1.0/36.0)*rho*(1.00 + sumsq2 + ui - vi - uv)
    return feq
def calc_force(rho):
    psi = ne.evaluate('rho_psi*(1-exp(-rho/rho_psi))')
    # a negative roll velocity implied i+1

    f_x = psi*(np.roll(psi,(-1,0),(0,1))-np.roll(psi,(1,0),(0,1)))
    f_y = psi*(np.roll(psi,(0,-1),(0,1))-np.roll(psi,(0,1),(0,1)))

    f2_x = psi*(np.roll(psi,(-1,1),(0,1))+np.roll(psi,(-1,-1),(0,1))-np.roll(psi,(1,1),(0,1))-np.roll(psi,(1,-1),(0,1)))


    f2_y = psi*(np.roll(psi,(-1,-1),(0,1))+np.roll(psi,(1,-1),(0,1))
                -np.roll(psi,(1,1),(0,1))-np.roll(psi,(-1,1),(0,1)))

    force_x = -G*(f_x/9+f2_x/36)
    force_y = -G*(f_y/9+f2_y/36)
    return force_x, force_y
def fast(rho, Nsteps):
    f = np.zeros([Nsteps, 9,lx,ly])
    u = np.zeros([lx,ly])
    v = np.zeros([lx,ly])

    for kk in range(0,9):
        f[0,kk,:,:] = w[kk]*rho

    for ts in range(1,Nsteps):
        # if ts%100==0:
        #     print(ts)
#         prog_bar(ts,Nsteps)
        fout = np.copy(f[ts-1])
        for kk in range(0,9):
            f[ts,kk,:,:] = np.roll(fout[kk,:,:],(ex[kk],ey[kk]),(0,1))

        u[:,:] = 0
        v[:,:] = 0

        rho = np.sum(f[ts],axis=0)
        for k in range(9):
            u += f[ts,k]*ex[k]
            v += f[ts, k]*ey[k]
        u=u/rho
        v=v/rho


        force_x, force_y = calc_force(rho)
        u += force_x/(omega*rho)
        v += force_y/(omega*rho)

        #=======================================
        cs2  = 1.00 / 3.00
        cs22 = 2.00 * cs2
        cssq = 2.0 / 9.00
        feq = equilibrium(rho,u,v)
        for k in range(9):
            f[ts,k] = feq[k]
    return f[-2]
def init_rho():
    rho= np.zeros([lx,ly])
    disp_x = lx/2
    disp_y = ly/2
    y,x = np.ogrid[-disp_x:lx-disp_x,-disp_y:ly-disp_y]
    radius = 20
    mask = x*x +y*y < radius**2
    rho[mask] = 2.4
    rho *= 1+.05*(np.random.rand(lx,ly)-.5)
    rho[mask == 0] = .125
    return rho
def determine_radius(rho):
    rho_cent = rho[int(lx/2),int(ly/2)]
    rho_out = rho[int(lx/8),int(ly/8)]
    target = (rho_cent+rho_out)/2
    idx = np.argmin(np.abs(rho[:,int(ly/2)] - target))
    return np.abs(lx/2-idx)
def pressure(rho):
    psi = rho_psi*(1-np.exp(-rho/rho_psi))
    press = (rho*cs2+.5*G*cs2*psi**2)
    pmax = np.max(press)
    pmin = np.min(press)
    return pmin, pmax

rho_init = init_rho()
plt.imshow(flip_T(rho_init))
plt.colorbar()
plt.grid('off')
plt.show()


a = fast(rho_init, Nsteps)
rho = calc_macro_f(a)[0]
print("Radius: {:}".format(determine_radius(rho)))
print("Pressures:")
print(pressure(rho))
plt.imshow(flip_T(rho))
plt.colorbar()
plt.grid('off')
plt.show()
