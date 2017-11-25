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
def add_circ_to_mask(mask,x,y,cent_x,cent_y,radius):
    new = (x-cent_x)**2 + (y-cent_y)**2 < radius**2
    return mask+new
def mask_to_rho(mask, lx, ly):
    rho= np.zeros([lx,ly])
    rho[mask !=0] = 2.4
    rho *= 1+.01*(np.random.rand(lx,ly)-.5)*2
    rho[mask==0] = .125
    return rho
class LBM():
    def __init__(self, G = -5, lx=128, ly=128, omega = 1,G1 = None, G2=0,
                init_type = 'rand',g = 0,g_theta=-np.pi/2, rho=None):
        self.lx = lx
        self.ly = ly
        self.G = G
        #Make 2 belt default to 1 belt if no G1 passed
        if G1 is not None:
            self.G1 = G1
        else:
            self.G1 = G
        self.G2 = G2

        #Gravity
        self.g = g
        self.g_theta = g_theta

        self.omega = 1
        self.f = np.zeros([2, 9,lx,ly])
        self.u = np.zeros([lx,ly])
        self.v = np.zeros([lx,ly])
        self.feq = np.zeros([9,lx,ly])

        # Intiliaze rho
        if rho is not None:
            self.rho=rho
        elif init_type.lower() == 'rand':
            self.rho = 1+.05*((np.random.rand(self.lx,self.ly)-.5)*2)
        elif init_type == 'mult_circ':
            y,x = np.ogrid[0:self.lx,0:self.ly]
            mask = np.zeros([self.lx,self.ly])
            mask = add_circ_to_mask(mask,x,y,25,100,40)
            mask = add_circ_to_mask(mask,x,y,82,52,26)
            mask = add_circ_to_mask(mask,x,y,25,20,30)
            self.rho = mask_to_rho(mask, self.lx, self.ly)
        else:
            Print("No input. Random initial rho")
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
        force_x, force_y = self.force_shan_chen()
        force_x += np.cos(self.g_theta)*self.g
        force_y += np.sin(self.g_theta)*self.g
        return force_x, force_y
    def force_shan_chen_double_belt(self):
        rho_psi = 1
        r = self.rho

        psi = ne.evaluate('rho_psi*(1-exp(-r/rho_psi))')
        # a negative roll velocity implied i+1

        f_x = psi*(np.roll(psi,(-1,0),(0,1))-np.roll(psi,(1,0),(0,1)))
        f_y = psi*(np.roll(psi,(0,-1),(0,1))-np.roll(psi,(0,1),(0,1)))

        f2_x = psi*(np.roll(psi,(-1,1),(0,1))+np.roll(psi,(-1,-1),(0,1))
                    -np.roll(psi,(1,1),(0,1))-np.roll(psi,(1,-1),(0,1)))
        f2_y = psi*(np.roll(psi,(-1,-1),(0,1))+np.roll(psi,(1,-1),(0,1))
                    -np.roll(psi,(1,1),(0,1))-np.roll(psi,(-1,1),(0,1)))
        w1 =4/21/3
        w2 =4/45/3

        w4 =1/60/3
        w5 =2/315/3
        w8 =1/5040/3

        w0=1-(4*w1+4*w2+4*w4+8*w5+4*w8)

        force_x = -self.G1*(f_x/9+f2_x/36)-self.G2*(f_x*w1+f2_x*w2)
        force_y = -self.G1*(f_y/9+f2_y/36)-self.G2*(f_y*w1+f2_y*w2)


        # Second belt of forces
        # Here assume that psi for G1 is same as for G2. Everyone makes this
        # assumption. 2 is from increased c for farther lattice point

        # 2nd belt on axes
        f3_x = 2*psi*(np.roll(psi,(-2,0),(0,1))-np.roll(psi,(2,0),(0,1)))
        f3_y = 2*psi*(np.roll(psi,(0,-2),(0,1))-np.roll(psi,(0,2),(0,1)))

        #2nd belt corners
        f4_x = 2*psi*(np.roll(psi,(-2,2),(0,1))+np.roll(psi,(-2,-2),(0,1))
                    -np.roll(psi,(2,2),(0,1))-np.roll(psi,(2,-2),(0,1)))
        f4_y = 2*psi*(np.roll(psi,(-2,-2),(0,1))+np.roll(psi,(2,-2),(0,1))
                    -np.roll(psi,(2,2),(0,1))-np.roll(psi,(-2,2),(0,1)))
        # 2nd belt remaining points
        f5_x = 2*psi*(np.roll(psi,(-2,1),(0,1))+np.roll(psi,(-2,-1),(0,1))
                    -np.roll(psi,(2,1),(0,1))-np.roll(psi,(2,-1),(0,1))
                    +np.roll(psi,(-1,2),(0,1))+np.roll(psi,(-1,-2),(0,1))
                    -np.roll(psi,(2,-1),(0,1))-np.roll(psi,(1,-2),(0,1)))
        # I don't want to type that out carefully again, so just swap the axis
        # argument, which is the same as s
        f5_y = 2*psi*(np.roll(psi,(-2,1),(1,0))+np.roll(psi,(-2,-1),(1,0))
                    -np.roll(psi,(2,1),(1,0))-np.roll(psi,(2,-1),(1,0))
                    +np.roll(psi,(-1,2),(1,0))+np.roll(psi,(-1,-2),(1,0))
                    -np.roll(psi,(2,-1),(1,0))-np.roll(psi,(1,-2),(1,0)))

        force_x += -self.G2*(f3_x*w4+f4_x*w8+f5_x*w5)
        force_y += -self.G2*(f3_y*w4+f4_y*w8+f5_y*w5)
        return force_x, force_y
    def force_shan_chen_double_belt_grav(self):
        force_x, force_y = self.force_shan_chen_double_belt()
        force_x += np.cos(self.g_theta)*self.g
        force_y += np.sin(self.g_theta)*self.g
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
                #TODO Check that fout is necessary
                #TODO check that the loop is necessary. One numpy roll?
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
