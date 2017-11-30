from visualization import Canvas as viz
from muphase_1belt_wall import LBM_1belt_wall
from vispy import app

# sim = LBM(G1 = -5, G2 = 1.5,g=10**-4, init_type='mult_circ')
sim = LBM_1belt_wall(G1 = -6, G2 = 1.5,g=10**-4, init_type='rand',boundaries='can')
# sim = LBM(G1 = -7, G2 = 2.65,g=10**-4, init_type='given',input_rho='rand_rho.npy')
# sim = LBM(G1 = -6.5, G2 = 0,g=10**-4, init_type='rand')
# sim = LBM(G1 = -6.5, G2 = 2,g=10**-4, init_type='rb')
c = viz(sim,force_method= sim.force_shan_chen_grav,steps_per_draw = 1)
# c = viz(sim,force_method= sim.force_shan_chen_grav,steps_per_draw = 1)

# c = viz(sim,force_method= sim.force_shan_chen,steps_per_draw = 1,
#         render_folder='render_prep/', save_images = True, max_steps = 25*60)

# c = viz(sim,force_method= sim.force_shan_chen_double_belt,steps_per_draw = 1,
#         render_folder='render_prep/', save_images = True, max_steps = 24*60)

c.measure_fps()

app.run()
