from visualization import Canvas as viz
from python_muphase_no_hist_obj import LBM
from vispy import app

sim = LBM(G1 = -6, G2 = 5,g=10**-4, init_type='mult_circ')
# c = viz(sim,force_method= sim.force_shan_chen_double_belt_grav,steps_per_draw = 1)
c = viz(sim,force_method= sim.force_shan_chen_double_belt,steps_per_draw = 1,
        render_folder='render_prep/', save_images = True)
c.measure_fps()

app.run()
