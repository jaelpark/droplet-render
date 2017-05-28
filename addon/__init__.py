
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
	bl_use_exclude_layers = True;

	def RenderTiles(self, f):
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
			tilew = self.tilew;
			if tile1[0]+self.tilew > self.width:
				tilew += int(self.width-(tile1[0]+self.tilew));
			tileh = self.tileh;
			if tile1[1]+self.tileh > self.height:
				tileh += int(self.height-(tile1[1]+self.tileh));

			tilew += int(min(tile1[0],0));
			tileh += int(min(tile1[1],0));

			tile1 = ((int(max(tile1[0],0)),int(max(tile1[1],0))));

			if not f((tile1[0],tile1[1],tilew,tileh),nx*ny-len(tiles),nx*ny):
				break;

	def update_render_passes(self, scene, srl):
		#intern/cycles/blender/addon/__init__.py
		#intern/cycles/blender/addon/engine.py
		self.register_pass(scene,srl,"Combined",4,"RGBA","COLOR");
		#self.register_pass(scene,srl,"Directional",3,"RGB","COLOR");
		#self.register_pass(scene,srl,"Environment",3,"RGB","COLOR");
		if srl.use_pass_transmission_direct:
			self.register_pass(scene,srl,"TransDir",3,"RGB","COLOR");
		if srl.use_pass_transmission_indirect:
			self.register_pass(scene,srl,"TransInd",3,"RGB","COLOR");
		if srl.use_pass_shadow:
			self.register_pass(scene,srl,"Shadow",3,"RGB","COLOR");

	def update(self, data, scene):
		self.samples_ext = scene.blcloudsampling.samples;
		self.samples_int = scene.blcloudperf.samples;
		self.width = scene.render.resolution_x;
		self.height = scene.render.resolution_y;
		self.tilew = scene.blcloudperf.tilex;
		self.tileh = scene.blcloudperf.tiley;
		self.layer = [x for x, m in enumerate(scene.render.layers) if scene.render.layers.active == m][0];
		self.smask = 0;
		self.primary = scene.render.layers.active.use_pass_combined or\
			scene.render.layers.active.use_pass_transmission_direct or scene.render.layers.active.use_pass_transmission_indirect;
		self.shadow = scene.render.layers.active.use_pass_shadow;
		for i,sclayer in enumerate(scene.layers):
			self.smask |= int(sclayer and scene.render.layers.active.layers[i])<<i;

		self.update_stats("Droplet","Initializing");
		libdroplet.BeginRender(scene,data,self.tilew,self.tileh,self.width,self.height,self.smask);
		while libdroplet.QueryStatus() != 0:
			pass; #TODO: query progress/memory usage etc

	def RenderScene(self, tile, index, total):
		sc = int(ceil(self.samples_ext/self.samples_int)); #external sample count

		result = self.begin_result(tile[0],tile[1],tile[2],tile[3]);
		cl = np.zeros((tile[2]*tile[3],4)); #light sources
		cs = cl.copy(); #environment

		dd = 0; #total accumulated sample count
		for i in range(0,sc):
			if self.test_break():
				self.end_result(result);
				return False;

			self.update_stats("Path tracing tile ("+str(index)+"/"+str(total)+")",str(dd)+"/"+str(self.samples_ext)+" samples");
			d1 = min(self.samples_ext-i*self.samples_int,self.samples_int);

			libdroplet.Render(tile[0],tile[1],tile[2],tile[3],d1);
			while True:
				qr = libdroplet.QueryResult(0);
				if qr is not None:
					break;

			dd += d1;
			cl += qr;
			cs += libdroplet.QueryResult(1);
			fd = 1.0/float(dd);

			rpass = result.layers[self.layer].passes.find_by_type("COMBINED",result.views[0].name); #find_by_name("Combined",result.views[0].name);
			if rpass is not None:
				rpass.rect = (cl+cs)*np.array([fd,fd,fd,0.5*fd]);
			rpass = result.layers[self.layer].passes.find_by_type("TRANSMISSION_DIRECT",result.views[0].name);
			if rpass is not None:
				rpass.rect = np.delete(cl,3,1)*fd;
			rpass = result.layers[self.layer].passes.find_by_type("TRANSMISSION_INDIRECT",result.views[0].name);
			if rpass is not None:
				rpass.rect = np.delete(cs,3,1)*fd;

			self.update_result(result);
			self.update_progress(1.0-(total-index)/(total)+i/(sc*total));

		self.end_result(result);
		return True;

	def RenderShadow(self, tile, index, total):
		#
		result = self.begin_result(tile[0],tile[1],tile[2],tile[3]);
		cl = np.zeros((tile[2]*tile[3],4)); #light sources
		cs = cl.copy(); #environment

		dd = 0;
		for i in range(0,1):
			if self.test_break():
				self.end_result(result);
				return False;

			self.update_stats("Shadowing tile ("+str(index)+"/"+str(total)+")","0/100 samples");
			d1 = 100;

			libdroplet.Shadow(tile[0],tile[1],tile[2],tile[3],d1);
			while True:
				qr = libdroplet.QueryResult(0);
				if qr is not None:
					break;

			dd += d1;
			fd = 1.0/float(dd);

			rpass = result.layers[self.layer].passes.find_by_type("SHADOW",result.views[0].name);
			if rpass is not None:
				rpass.rect = qr*np.array([fd,fd,fd,0.5*fd]); #TODO: channels?

			self.update_result(result);
			self.update_progress(1.0-(total-index)/(total)+i/(sc*total));

		self.end_result(result);
		return True;

	def render(self, scene):
		if self.primary:
			self.RenderTiles(self.RenderScene);
		if self.shadow:
			self.RenderTiles(self.RenderShadow);
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
	#bpy.types.Scene.blcloudpasses = PointerProperty(type=panel.ClPassProperties);
	bpy.types.World.droplet = PointerProperty(type=panel.ClEnvironmentProperties);

	bpy.types.Object.droplet = PointerProperty(type=panel.ClObjectProperties);
	bpy.types.ParticleSettings.droplet = PointerProperty(type=panel.ClParticleSystemProperties);
	bpy.types.Lamp.droplet = PointerProperty(type=panel.ClLampProperties);

	nodeitems_utils.register_node_categories("BLCLOUD_CATEGORIES",node.categories);

def unregister():
	bpy.utils.unregister_module(__name__);

	nodeitems_utils.unregister_node_categories("BLCLOUD_CATEGORIES");

if __name__ == "__main__":
	register();
