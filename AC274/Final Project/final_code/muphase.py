from visualization import Canvas as viz
from python_muphase_no_hist_obj import LBM
from vispy import app

# sim = LBM(G1 = -6.5, G2 = 2,g=10**-4, init_type='mult_circ')
sim = LBM(G1 = -6, G2 = .65,g=10**-4, init_type='rand')
# sim = LBM(G1 = -6.5, G2 = 0,g=10**-4, init_type='rand')
# sim = LBM(G1 = -6.5, G2 = 2,g=10**-4, init_type='rb')
# c = viz(sim,force_method= sim.force_shan_chen,steps_per_draw = 1)
# c = viz(sim,force_method= sim.force_shan_chen_double_belt,steps_per_draw = 1)
c = viz(sim,force_method= sim.force_shan_chen,steps_per_draw = 1,
        render_folder='render_prep/', save_images = True, max_steps = 24*60)
c.measure_fps()

app.run()
