
import bpy

#unix:
#/usr/share/blender/2.77/scripts/addons/render_droplet
#/usr/lib64/python3.5/site-packages

from math import ceil
from mathutils import Vector
from bl_ui import properties_render

import nodeitems_utils
from bpy.props import PointerProperty

from bl_ui import properties_render_layer,properties_data_mesh,properties_data_camera,properties_particle,properties_physics_common,properties_physics_field,properties_physics_smoke

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

	def update(self, data, scene):
		self.samples_ext = scene.blcloudsampling.samples;
		self.samples_int = scene.blcloudperf.samples;
		self.width = scene.render.resolution_x;
		self.height = scene.render.resolution_y;
		self.tilew = scene.blcloudperf.tilex;
		self.tileh = scene.blcloudperf.tiley;
		self.layer = [x for x, m in enumerate(scene.render.layers) if scene.render.layers.active == m][0];
		self.smask = 0;
		for i,sclayer in enumerate(scene.layers):
			self.smask |= int(sclayer and scene.render.layers.active.layers[i])<<i;

		self.update_stats("Droplet","Initializing");
		libdroplet.BeginRender(scene,data,self.tilew,self.tileh,self.width,self.height,self.smask);
		while libdroplet.QueryStatus() != 0:
			pass; #TODO: query progress/memory usage etc

	# def DrawBorder(self, result, tilew, tileh, tilew1, tileh1):
	# 	bcolor = [0.9,0.5,0.0,1.0];
	# 	for x in range(0,tilew1):
	# 		result.layers[0].passes[0].rect[x] = bcolor;
	# 		result.layers[0].passes[0].rect[tileh*(tileh1-1)+x] = bcolor;
	# 	for y in range(0,tileh1):
	# 		result.layers[0].passes[0].rect[tileh*y] = bcolor;
	# 		result.layers[0].passes[0].rect[tileh*y+(tilew1-1)] = bcolor;

	def render(self, scene):
		sc = int(ceil(self.samples_ext/self.samples_int)); #external sample count
		nx = int(ceil(self.width/self.tilew));
		ny = int(ceil(self.height/self.tileh));

		tiles = [];
		for y in range(0,ny):
			for x in range(0,nx):
				xadj = x*self.tilew-0.5*(nx*self.tilew-self.width);
				yadj = y*self.tileh-0.5*(ny*self.tileh-self.height);
				tiles.append((xadj,yadj));

		while len(tiles) > 0:
			tile1 = min(tiles,key=lambda tileq: Vector((
				tileq[0]+0.5*self.tilew-0.5*self.width,
				tileq[1]+0.5*self.tileh-0.5*self.height)).length);
			tiles.remove(tile1);

			#crop the tile frame on edges
			tilew1 = self.tilew;
			if tile1[0]+self.tilew > self.width:
				tilew1 += int(self.width-(tile1[0]+self.tilew));
			tileh1 = self.tileh;
			if tile1[1]+self.tileh > self.height:
				tileh1 += int(self.height-(tile1[1]+self.tileh));

			tilew1 += int(min(tile1[0],0));
			tileh1 += int(min(tile1[1],0));

			tile1 = ((int(max(tile1[0],0)),int(max(tile1[1],0))));

			result = self.begin_result(tile1[0],tile1[1],tilew1,tileh1);
			cl = np.zeros((tilew1*tileh1,4)); #light sources
			cs = cl.copy(); #environment

			#self.DrawBorder(result,self.tilew,self.tileh,tilew1,tileh1);
			#self.update_result(result);

			dd = 0; #total accumulated sample count
			for i in range(0,sc):
				if self.test_break():
					self.end_result(result);
					break;
				self.update_stats("Path tracing tile ("+str((nx*ny)-len(tiles))+"/"+str(nx*ny)+")",str(dd)+"/"+str(self.samples_ext)+" samples");
				d1 = min(self.samples_ext-i*self.samples_int,self.samples_int);

				libdroplet.Render(tile1[0],tile1[1],tilew1,tileh1,d1);

				while True:
					qr = libdroplet.QueryResult(0);
					if qr is not None:
						break;
				dd += d1;
				cl += qr;
				cs += libdroplet.QueryResult(1);

				rpass = next((x for x in result.layers[self.layer].passes if x.type == 'COMBINED'),None);
				if rpass is not None:
					rpass.rect = (cl+cs)/float(dd);
				rpass = next((x for x in result.layers[self.layer].passes if x.type == 'DIFFUSE'),None);
				if rpass is not None:
					rpass.rect = np.delete(cl/float(dd),3,1);
				rpass = next((x for x in result.layers[self.layer].passes if x.type == 'ENVIRONMENT'),None);
				if rpass is not None:
					rpass.rect = np.delete(cs/float(dd),3,1);

				#if i < sc-1:
					#self.DrawBorder(result,self.tilew,self.tileh,tilew1,tileh1);

				self.update_result(result);
				self.update_progress(1.0-len(tiles)/(nx*ny)+i/(sc*nx*ny));
			else:
				self.end_result(result);
				continue;
			break;

			#self.end_result(result);
			#self.update_progress((y*nx+x+1)/(ny*nx));
			#https://www.blender.org/api/blender_python_api_2_76_release/bpy.types.RenderEngine.html
			#self.update_memory_stats

		libdroplet.EndRender();

def register():
	bpy.utils.register_module(__name__);
	properties_render_layer.RENDERLAYER_PT_layers.COMPAT_ENGINES.add(config.dre_engineid);

	properties_data_mesh.DATA_PT_context_mesh.COMPAT_ENGINES.add(config.dre_engineid);
	properties_data_mesh.DATA_PT_texture_space.COMPAT_ENGINES.add(config.dre_engineid);
	properties_data_mesh.DATA_PT_vertex_groups.COMPAT_ENGINES.add(config.dre_engineid);
	properties_data_mesh.DATA_PT_shape_keys.COMPAT_ENGINES.add(config.dre_engineid);
	properties_data_mesh.DATA_PT_uv_texture.COMPAT_ENGINES.add(config.dre_engineid);
	properties_data_mesh.DATA_PT_vertex_colors.COMPAT_ENGINES.add(config.dre_engineid);
	properties_data_mesh.DATA_PT_customdata.COMPAT_ENGINES.add(config.dre_engineid);

	properties_data_camera.DATA_PT_lens.COMPAT_ENGINES.add(config.dre_engineid);
	properties_data_camera.DATA_PT_camera_display.COMPAT_ENGINES.add(config.dre_engineid);

	properties_particle.PARTICLE_PT_context_particles.COMPAT_ENGINES.add(config.dre_engineid);
	properties_particle.PARTICLE_PT_emission.COMPAT_ENGINES.add(config.dre_engineid);
	properties_particle.PARTICLE_PT_draw.COMPAT_ENGINES.add(config.dre_engineid);
	properties_particle.PARTICLE_PT_velocity.COMPAT_ENGINES.add(config.dre_engineid);
	properties_particle.PARTICLE_PT_physics.COMPAT_ENGINES.add(config.dre_engineid);
	properties_particle.PARTICLE_PT_field_weights.COMPAT_ENGINES.add(config.dre_engineid);
	properties_particle.PARTICLE_PT_force_fields.COMPAT_ENGINES.add(config.dre_engineid);
	properties_particle.PARTICLE_PT_vertexgroups.COMPAT_ENGINES.add(config.dre_engineid);
	properties_particle.PARTICLE_PT_custom_props.COMPAT_ENGINES.add(config.dre_engineid);

	properties_physics_common.PHYSICS_PT_add.COMPAT_ENGINES.add(config.dre_engineid);

	properties_physics_field.PHYSICS_PT_field.COMPAT_ENGINES.add(config.dre_engineid);

	properties_physics_smoke.PHYSICS_PT_smoke.COMPAT_ENGINES.add(config.dre_engineid);
	properties_physics_smoke.PHYSICS_PT_smoke_highres.COMPAT_ENGINES.add(config.dre_engineid);
	properties_physics_smoke.PHYSICS_PT_smoke_groups.COMPAT_ENGINES.add(config.dre_engineid);
	properties_physics_smoke.PHYSICS_PT_smoke_cache.COMPAT_ENGINES.add(config.dre_engineid);
	properties_physics_smoke.PHYSICS_PT_smoke_field_weights.COMPAT_ENGINES.add(config.dre_engineid);

	bpy.types.Scene.blcloudrender = PointerProperty(type=panel.ClRenderProperties);
	bpy.types.Scene.blcloudsampling = PointerProperty(type=panel.ClSamplingProperties);
	bpy.types.Scene.blcloudgrid = PointerProperty(type=panel.ClGridProperties);
	bpy.types.Scene.blcloudperf = PointerProperty(type=panel.ClPerformanceProperties);
	bpy.types.World.droplet = PointerProperty(type=panel.ClCompositeProperties);

	bpy.types.Object.droplet = PointerProperty(type=panel.ClObjectProperties);
	bpy.types.ParticleSettings.droplet = PointerProperty(type=panel.ClParticleSystemProperties);
	bpy.types.Lamp.droplet = PointerProperty(type=panel.ClLampProperties);

	nodeitems_utils.register_node_categories("BLCLOUD_CATEGORIES",node.categories);

def unregister():
	bpy.utils.unregister_module(__name__);

	nodeitems_utils.unregister_node_categories("BLCLOUD_CATEGORIES");

if __name__ == "__main__":
	register();
