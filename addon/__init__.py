
import bpy

#unix:
#/usr/share/blender/2.77/scripts/addons/render_droplet
#/usr/lib64/python3.5/site-packages

from math import ceil
from mathutils import Vector
from bl_ui import properties_render

import nodeitems_utils
from bpy.props import PointerProperty

from bl_ui import properties_particle,properties_physics_common,properties_physics_field,properties_physics_smoke

#NOTE: addon directory name cannot coincide with pyd module name
import libdroplet #libblcloud #pyd/so filename
import numpy as np

if "bpy" in locals():
	import imp
	if "panel" in locals():
		imp.reload(panel);
	if "node" in locals():
		imp.reload(node);
	if "config" in locals():
		imp.reload(config);

from . import panel
from . import node
from . import config

bl_info = {
	"name":"Droplet Render",
	"author":"jaelpark",
	"version":(0,0,1),
	"blender":(2,77,0), #for python 3.5
	#"support":"COMMUNITY",
	"warning":"",
	"location":"Info header, render engine menu",
	"description":"Experimental volumetric cloud modeling and rendering engine",
	"category":"Render"
}

class CloudRenderEngine(bpy.types.RenderEngine):
	bl_idname = config.dre_engineid;
	bl_label = "Droplet Render";
	bl_use_preview = False;
	#outnode = None;

	#def update(self, data, scene):
		#group = data.node_groups[scene.blcloudrender.nodetree];
		#self.outnode = group.nodes["Material Output"]; #warning, unreliable: fix
		#print("INFO: scene update requested");

	def render(self, scene):
		#s = scene.render.resolution_percentage/100.0;
		#w = int(s*scene.render.resolution_x);
		#h = int(s*scene.render.resolution_y);
		w = max(int(scene.blcloudrender.res_p*scene.blcloudrender.res_x),1);
		h = max(int(scene.blcloudrender.res_p*scene.blcloudrender.res_y),1);
		scene.render.resolution_x = w;#scene.blcloudrender.res_x;
		scene.render.resolution_y = h;#scene.blcloudrender.res_y;
		scene.render.resolution_percentage = 100;#*scene.blcloudrender.res_p;

		rx = scene.blcloudperf.tilex;
		ry = scene.blcloudperf.tiley; #self.tile_x, tile_y
		st = scene.blcloudsampling.samples;
		sr = scene.blcloudperf.samples;
		sc = int(ceil(st/sr)); #external sample count

		ny = int(ceil(h/ry));
		nx = int(ceil(w/rx));

		self.update_stats("Droplet","Initializing");

		tl = [];
		for y in range(0,ny):
			for x in range(0,nx):
				tl.append((x*rx,y*ry)); #TODO: fix the alignment so that the first tile starts at the center

		#libdroplet.BeginRender(scene,self.outnode,rx,ry,w,h);
		libdroplet.BeginRender(scene,bpy.data,rx,ry,w,h);

		while len(tl) > 0:
			tc = min(tl,key=lambda tt: Vector((tt[0]+0.5*rx-0.5*w,tt[1]+0.5*ry-0.5*h)).length)
			tl.remove(tc);

			result = self.begin_result(tc[0],tc[1],rx,ry);
			rr = np.zeros((rx*ry,4));

			dd = 0;
			for i in range(0,sc):
				if self.test_break():
					break; #TODO: end_result to keep all the work
				self.update_stats("Path tracing tile ("+str((nx*ny)-len(tl))+"/"+str(nx*ny)+")",str(dd)+"/"+str(st)+" samples");
				d1 = min(st-i*sr,sr);

				dd += d1;
				rr += libdroplet.Render(tc[0],tc[1],d1);

				result.layers[0].passes[0].rect = rr/float(dd);
				self.update_result(result); #TODO: draw rectangle if i < sc-1
				self.update_progress(1.0-len(tl)/(nx*ny)+i/(sc*nx*ny));
			else:
				self.end_result(result);
				continue;
			break;

			#self.end_result(result);
			#self.update_progress((y*nx+x+1)/(ny*nx));
			#https://www.blender.org/api/blender_python_api_2_76_release/bpy.types.RenderEngine.html
			#self.update_memory_stats

		libdroplet.EndRender();

#bpy.utils.register_class(CloudRenderEngine)

def register():
	bpy.utils.register_module(__name__);
	#properties_render.RENDER_PT_render.COMPAT_ENGINES.add('droplet');
	properties_particle.PARTICLE_PT_context_particles.COMPAT_ENGINES.add(config.dre_engineid);
	properties_particle.PARTICLE_PT_emission.COMPAT_ENGINES.add(config.dre_engineid);
	properties_particle.PARTICLE_PT_draw.COMPAT_ENGINES.add(config.dre_engineid);

	properties_physics_common.PHYSICS_PT_add.COMPAT_ENGINES.add(config.dre_engineid);

	properties_physics_field.PHYSICS_PT_field.COMPAT_ENGINES.add(config.dre_engineid);

	properties_physics_smoke.PHYSICS_PT_smoke.COMPAT_ENGINES.add(config.dre_engineid);
	properties_physics_smoke.PHYSICS_PT_smoke_highres.COMPAT_ENGINES.add(config.dre_engineid);
	properties_physics_smoke.PHYSICS_PT_smoke_groups.COMPAT_ENGINES.add(config.dre_engineid);
	properties_physics_smoke.PHYSICS_PT_smoke_cache.COMPAT_ENGINES.add(config.dre_engineid);

	bpy.types.Scene.blcloudrender = PointerProperty(type=panel.ClRenderProperties);
	bpy.types.Scene.blcloudsampling = PointerProperty(type=panel.ClSamplingProperties);
	bpy.types.Scene.blcloudgrid = PointerProperty(type=panel.ClGridProperties);
	bpy.types.Scene.blcloudperf = PointerProperty(type=panel.ClPerformanceProperties);

	bpy.types.Object.droplet = PointerProperty(type=panel.ClObjectProperties);
	#bpy.types.ParticleSystem.droplet = PointerProperty(type=panel.ClParticleSystemProperties);
	bpy.types.Lamp.droplet = PointerProperty(type=panel.ClLampProperties);

	nodeitems_utils.register_node_categories("BLCLOUD_CATEGORIES",node.categories);

def unregister():
	bpy.utils.unregister_module(__name__);

	nodeitems_utils.unregister_node_categories("BLCLOUD_CATEGORIES");

if __name__ == "__main__":
	register();
