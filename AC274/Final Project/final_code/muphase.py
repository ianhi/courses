from visualization import Canvas as viz
from python_muphase_no_hist_obj import LBM
from vispy import app

sim = LBM()
c = viz(sim, force_method= sim.force_shan_chen)
c.measure_fps()
app.run()
