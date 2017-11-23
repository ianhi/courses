import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as animation
import numexpr as ne
mpl.style.use(['fivethirtyeight', 'figsize'])


# define constants
# lx =128
# ly = lx
# Nsteps = 5000
ex=[0,1,0,-1,0 ,1,-1,-1, 1]
ey=[0,0,1,0 ,-1,1, 1,-1,-1]
w = [4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36]
# G = -5
omega = 1
rho_psi = 1

cs = 1/np.sqrt(3)
cs2  = 1.00 / 3.00
cs22 = 2.00 * cs2
cssq = 2.0 / 9.00
#plotting functions:
def flip_T(a):
    return np.flipud(a.T)
def calc_macro_f(f,lx=128,ly=128):
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
def equilibrium(feq, rho,u,v):
    feq[:,:,:]=0
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
def shan_chen(rho,G):
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

def shan_chen_grav(rho,G):
    psi = ne.evaluate('rho_psi*(1-exp(-rho/rho_psi))')
    # a negative roll velocity implied i+1

    f_x = psi*(np.roll(psi,(-1,0),(0,1))-np.roll(psi,(1,0),(0,1)))
    f_y = psi*(np.roll(psi,(0,-1),(0,1))-np.roll(psi,(0,1),(0,1)))

    f2_x = psi*(np.roll(psi,(-1,1),(0,1))+np.roll(psi,(-1,-1),(0,1))-np.roll(psi,(1,1),(0,1))-np.roll(psi,(1,-1),(0,1)))


    f2_y = psi*(np.roll(psi,(-1,-1),(0,1))+np.roll(psi,(1,-1),(0,1))
                -np.roll(psi,(1,1),(0,1))-np.roll(psi,(-1,1),(0,1)))

    force_x = -G*(f_x/9+f2_x/36)
    force_y = -G*(f_y/9+f2_y/36) + .2*rho
    return force_x, force_y
def double_belt_shan_chen(rho,G):
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
def fast(rho, Nsteps,G, lx,ly,force_method):
    f = np.zeros([2, 9,lx,ly])
    u = np.zeros([lx,ly])
    v = np.zeros([lx,ly])
    feq = np.zeros([9,lx,ly])

    for kk in range(0,9):
        f[0,kk,:,:] = w[kk]*rho

    for ts in range(1,Nsteps):
        fout = np.copy(f[0])
        for kk in range(0,9):
            f[1,kk,:,:] = np.roll(fout[kk,:,:],(ex[kk],ey[kk]),(0,1))

        u[:,:] = 0
        v[:,:] = 0

        rho = np.sum(f[1],axis=0)
        for k in range(9):
            u += f[1,k]*ex[k]
            v += f[1, k]*ey[k]
        u=u/rho
        v=v/rho


        force_x, force_y = force_method(rho,G)
        u += force_x/(omega*rho)
        v += force_y/(omega*rho)

        #=======================================
        cs2  = 1.00 / 3.00
        cs22 = 2.00 * cs2
        cssq = 2.0 / 9.00
        feq = equilibrium(feq, rho,u,v)
        for k in range(9):
            f[0,k] = feq[k]
    return f
def init_rho(radius=20, lx=128, ly = 128):
    rho= np.zeros([lx,ly])
    disp_x = lx/2
    disp_y = ly/2
    y,x = np.ogrid[-disp_x:lx-disp_x,-disp_y:ly-disp_y]
    mask = x*x +y*y < radius**2
    rho[mask] = 2.4
    rho *= 1+.01*(np.random.rand(lx,ly)-.5)*2
    rho[mask == 0] = .125
    return rho
def determine_radius(rho, lx=128,ly=128):
    rho_cent = rho[int(lx/2),int(ly/2)]
    rho_out = rho[int(lx/16),int(ly/16)]
    target = (rho_cent+rho_out)/2
    idx = np.argmin(np.abs(rho[:,int(ly/2)] - target))
    return np.abs(lx/2-idx)
def determine_radius2(rho):
    lx = rho.shape[0]
    ly = rho.shape[1]
    m_tot = lx*ly*np.average(rho)
    rho_min = np.min(rho)
    m_drop = m_tot - lx*ly*rho_min
    rho_drop = np.max(rho)
    return np.sqrt(m_drop/(np.pi*rho_drop))

def pressure(rho,G):
    psi = rho_psi*(1-np.exp(-rho/rho_psi))
    press = (rho*cs2+.5*G*cs2*psi**2)
    pmax = np.max(press)
    pmin = np.min(press)
    return pmin, pmax
def pressure2(rho,G,lx=128,ly=128):
    psi = rho_psi*(1-np.exp(-rho/rho_psi))
    press = (rho*cs2+.5*G*cs2*psi**2)
    pmid = press[int(lx/2),int(ly/2)]
    pout = press[int(lx/16),int(ly/16)]
    return pmid-pout
def run(G = -5, radius_init = 20, Nsteps = 5000, save = True,
                                    savename=None, lx=128,ly=128):
    rho_init = init_rho(radius_init, lx,ly)
    out = fast(rho_init,Nsteps, G, lx, ly,shan_chen)

    rho = calc_macro_f(out[1],lx,ly)[0]
    radius = determine_radius(rho,lx,ly)
    radius2 = determine_radius2(rho)
    press = pressure(rho,G)
    delta_P2 = pressure2(rho,G,lx,ly)
    print("G: {:}".format(G))
    print("Radius: {:}".format(radius))
    print("Radius2: {:}".format(radius2))
    print("Delta P: {:}".format(press[1]-press[0]))
    print("Pmin: {:}".format(press[0]))
    print("Pmax: {:}".format(press[1]))

    if save:
        if savename is None:
            savename = "muphase_{:}_{:}_{:}".format(G,Nsteps,radius_init)
        np.save(savename,rho)
    out_dict = {'G':G,'radius':radius,'radius2':radius2,'pmin':press[0],'pmax':press[1],'f':out[1],
                'deltaP':press[1]-press[0],'rho':rho,'deltaP2':delta_P2}
    return out_dict
class LBM():
    def __init__(self, G= -5, lx=128, ly=128, omega = 1):
        self.lx = lx
        self.ly = ly
        self.G = G
        self.omega = 1
        self.f = np.zeros([2, 9,lx,ly])
        self.u = np.zeros([lx,ly])
        self.v = np.zeros([lx,ly])
        self.feq = np.zeros([9,lx,ly])

        self.rho = 1+.05*((np.random.rand(self.lx,self.ly)-.5)*2)

        for k in range(0,9):
            self.f[0,k,:,:] = w[k]*self.rho
    def force_shan_chen(self):
        rho_psi = 1
        r = self.rho

        psi = ne.evaluate('rho_psi*(1-exp(-r/rho_psi))')
        # a negative roll velocity implied i+1

        f_x = psi*(np.roll(psi,(-1,0),(0,1))-np.roll(psi,(1,0),(0,1)))
        f_y = psi*(np.roll(psi,(0,-1),(0,1))-np.roll(psi,(0,1),(0,1)))

        f2_x = psi*(np.roll(psi,(-1,1),(0,1))+np.roll(psi,(-1,-1),(0,1))-np.roll(psi,(1,1),(0,1))-np.roll(psi,(1,-1),(0,1)))


        f2_y = psi*(np.roll(psi,(-1,-1),(0,1))+np.roll(psi,(1,-1),(0,1))
                    -np.roll(psi,(1,1),(0,1))-np.roll(psi,(-1,1),(0,1)))

        force_x = -self.G*(f_x/9+f2_x/36)
        force_y = -self.G*(f_y/9+f2_y/36)
        return force_x, force_y
    def force_shan_chen_grav(self):
        rho_psi = 1
        r = self.rho

        psi = ne.evaluate('rho_psi*(1-exp(-r/rho_psi))')
        # a negative roll velocity implied i+1

        f_x = psi*(np.roll(psi,(-1,0),(0,1))-np.roll(psi,(1,0),(0,1)))
        f_y = psi*(np.roll(psi,(0,-1),(0,1))-np.roll(psi,(0,1),(0,1)))

        f2_x = psi*(np.roll(psi,(-1,1),(0,1))+np.roll(psi,(-1,-1),(0,1))-np.roll(psi,(1,1),(0,1))-np.roll(psi,(1,-1),(0,1)))


        f2_y = psi*(np.roll(psi,(-1,-1),(0,1))+np.roll(psi,(1,-1),(0,1))
                    -np.roll(psi,(1,1),(0,1))-np.roll(psi,(-1,1),(0,1)))

        force_x = -self.G*(f_x/9+f2_x/36)
        force_y = -self.G*(f_y/9+f2_y/36) +10**-4
        return force_x, force_y

    def _equilibrium(self):
        usq = self.u*self.u
        vsq = self.v*self.v
        sumsq = (usq+vsq)/cs22
        sumsq2 = sumsq*(1-cs2)/cs2
        u22 = usq/cssq
        v22 = vsq/cssq
        ui = self.u/cs2
        vi = self.v/cs2
        uv = ui*vi
        self.feq[0] = (4.0/9.0)*self.rho*(1.00 - sumsq)

        self.feq[1] = (1.0/9.0)*self.rho*(1.00 - sumsq + u22 + ui)
        self.feq[2] = (1.0/9.0)*self.rho*(1.00 - sumsq + v22 + vi)
        self.feq[3] = (1.0/9.0)*self.rho*(1.00 - sumsq + u22 - ui)
        self.feq[4] = (1.0/9.0)*self.rho*(1.00 - sumsq + v22 - vi)

        self.feq[5] = (1.0/36.0)*self.rho*(1.00 + sumsq2 + ui + vi + uv)
        self.feq[6] = (1.0/36.0)*self.rho*(1.00 + sumsq2 - ui + vi - uv)
        self.feq[7] = (1.0/36.0)*self.rho*(1.00 + sumsq2 - ui - vi + uv)
        self.feq[8] = (1.0/36.0)*self.rho*(1.00 + sumsq2 + ui - vi - uv)
    def step(self,Nsteps,force_method=force_shan_chen):
        for ts in range(0,Nsteps):
            fout = np.copy(self.f[0])
            for kk in range(0,9):
                self.f[1,kk,:,:] = np.roll(fout[kk,:,:],(ex[kk],ey[kk]),(0,1))

            self.u[:,:] = 0
            self.v[:,:] = 0

            self.rho = np.sum(self.f[1],axis=0)
            for k in range(9):
                self.u += self.f[1,k]*ex[k]/self.rho
                self.v += self.f[1, k]*ey[k]/self.rho
            # u /= rho
            # v /= rho


            force_x, force_y = force_method()
            self.u += force_x/(self.omega*self.rho)
            self.v += force_y/(self.omega*self.rho)

            #=======================================
            cs2  = 1.00 / 3.00
            cs22 = 2.00 * cs2
            cssq = 2.0 / 9.00
            self._equilibrium()
            for k in range(9):
                self.f[0,k] = self.feq[k]
            # Calc macro quantites for output
            self.u[:,:] = 0
            self.v[:,:] = 0

            self.rho = np.sum(self.f[1],axis=0)
            for k in range(9):
                self.u += self.f[1,k]*ex[k]/self.rho
                self.v += self.f[1, k]*ey[k]/self.rho
saving = True
if __name__ == "__main__":
    sim = LBM()
    sim.step(1)
    plt.imshow(sim.rho,cmap=plt.cm.viridis)
    plt.grid('off')
    cb = plt.colorbar()
    cb.set_label('Rho')
    ax = plt.gca()
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    plt.tight_layout()
    if saving:
        plt.savefig('Figures/single_belt_init.png')
    plt.show()

    sim.step(100)
    plt.imshow(sim.rho,cmap=plt.cm.viridis)
    plt.grid('off')
    cb = plt.colorbar()
    cb.set_label('Rho')
    ax = plt.gca()
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    plt.tight_layout()
    if saving:
        plt.savefig('Figures/single_belt_100.png')
    plt.show()

    sim.step(1000)
    plt.imshow(sim.rho,cmap=plt.cm.viridis)
    plt.grid('off')
    cb = plt.colorbar()
    cb.set_label('Rho')
    ax = plt.gca()
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    plt.tight_layout()
    if saving:
        plt.savefig('Figures/single_belt_1000.png')
    plt.show()
