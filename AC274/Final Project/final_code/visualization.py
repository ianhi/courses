#!/usr/bin/env python
# -*- coding: utf-8 -*-
# vispy: gallery 2
"""
Example demonstrating showing a, image with a fixed ratio.
"""

import numpy as np

from vispy.util.transforms import ortho
from vispy import gloo
from vispy import app
from vispy import io
import matplotlib.pyplot as plt




VERT_SHADER = """
// Uniforms
uniform mat4 u_model;
uniform mat4 u_view;
uniform mat4 u_projection;
uniform float u_antialias;

// Attributes
attribute vec2 a_position;
attribute vec2 a_texcoord;

// Varyings
varying vec2 v_texcoord;

// Main
void main (void)
{
    v_texcoord = a_texcoord;
    gl_Position = u_projection * u_view * u_model * vec4(a_position,0.0,1.0);
}
"""

FRAG_SHADER = """
uniform sampler2D u_texture;
varying vec2 v_texcoord;

uniform float scale_factor;
uniform float max_magnitude;
uniform sampler1D colormap_array;

void main()
{
    float i_value = texture2D(u_texture, v_texcoord).r;
    float original = i_value;
    i_value *= scale_factor;

    // Calculate the position of i in the colormap
    if (i_value <= 0 ){
        gl_FragColor = texture1D(colormap_array, 0.000001);
    }
    else if (i_value >= max_magnitude){
        gl_FragColor = texture1D(colormap_array, 0.9999999);
    }
    else {
        float color_value = i_value;//(i_value + max_magnitude)/(2*max_magnitude);
        gl_FragColor = texture1D(colormap_array, color_value);
    }
}

"""


class Canvas(app.Canvas):

    def __init__(self,sim,scaling_factor=1.0, max_magnitude=1.0,force_method=None,
                steps_per_draw = 1, render_folder = "./",
                save_images = False,max_steps = 2147483647):
        self.sim = sim
        self.steps_per_draw = steps_per_draw
        self.render_folder = render_folder
        self.W, self.H = sim.lx, sim.ly
        self.total_steps = 0
        self.save_images = save_images
        if force_method is None:
            self.force_method = self.sim.force_shan_chen
        else:
            self.force_method = force_method
        self.max_steps = max_steps

        app.Canvas.__init__(self, keys='interactive', size=((self.W * 5), (self.H * 5)))
        self.scaling_factor = scaling_factor
        self.max_magnitude = max_magnitude

        self.I = np.zeros((self.W, self.H), dtype=np.float32, order='F')
        self.I = self.sim.rho

        # A simple texture quad
        self.data = np.zeros(4, dtype=[('a_position', np.float32, 2),
                                  ('a_texcoord', np.float32, 2)])
        self.data['a_position'] = np.array([[0, 0], [self.W, 0], [0, self.H], [self.W, self.H]])
        self.data['a_texcoord'] = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])

        self.program = gloo.Program(VERT_SHADER, FRAG_SHADER)
        self.texture = gloo.Texture2D(self.sim.rho, interpolation='nearest',internalformat='r32f')

        self.program['u_texture'] = self.texture
        self.program.bind(gloo.VertexBuffer(self.data))

        self.view = np.eye(4, dtype=np.float32)
        self.model = np.eye(4, dtype=np.float32)
        self.projection = np.eye(4, dtype=np.float32)

        self.program['u_model'] = self.model
        self.program['u_view'] = self.view
        self.projection = ortho(0, self.W, 0, self.H, -1, 1)
        self.program['u_projection'] = self.projection
        self.program['scale_factor'] = self.scaling_factor
        self.program['max_magnitude'] = self.max_magnitude

        self.cmap = plt.cm.viridis
        norm = plt.Normalize(-self.max_magnitude, self.max_magnitude)
        self.num_colors = 1024
        possible_values = np.linspace(-self.max_magnitude, self.max_magnitude, 1024)
        self.colormap_array = self.cmap(norm(possible_values)).astype(np.float32)

        self.program['colormap_array'] = self.colormap_array
        self.program['colormap_array'].interpolation = 'nearest'

        gloo.set_clear_color('white')

        self._timer = app.Timer('auto', connect=self.update, start=True)

        self.show()

    def on_resize(self, event):
        width, height = event.physical_size
        gloo.set_viewport(0, 0, width, height)
        self.projection = ortho(0, width, 0, height, -100, 100)
        self.program['u_projection'] = self.projection

        # Compute thje new size of the quad
        r = width / float(height)
        R = self.W / float(self.H)
        if r < R:
            w, h = width, width / R
            x, y = 0, int((height - h) / 2)
        else:
            w, h = height * R, height
            x, y = int((width - w) / 2), 0
        self.data['a_position'] = np.array(
            [[x, y], [x + w, y], [x, y + h], [x + w, y + h]])
        self.program.bind(gloo.VertexBuffer(self.data))

    def on_draw(self, event):
        gloo.clear(color=True, depth=True)
        self.sim.step(self.steps_per_draw,force_method=self.force_method)
        # print(np.sum(self.sim.rho))
        self.texture.set_data(((self.sim.rho-np.min(self.sim.rho))/np.max(self.sim.rho-np.min(self.sim.rho))).astype(np.float32))
        self.program.draw('triangle_strip')
        # self.update()
        if self.save_images:
            if self.total_steps < self.max_steps:
                self.total_steps += 1
                screenshot = gloo.util._screenshot()
                io.write_png(self.render_folder + self.force_method.__name__+
                        str(self.total_steps).zfill(6) + '.png', screenshot)
            else:
                print("DONE SAVING :)")
                self.save_images = False


if __name__ == '__main__':
    c = Canvas()
    c.measure_fps()
    app.run()
