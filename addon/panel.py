
import bpy

from bpy.props import BoolProperty,IntProperty,FloatProperty,EnumProperty,PointerProperty,FloatVectorProperty,StringProperty

from . import config

# def SetResolution(self, context):
# 	self.width = max(int(context.scene.blcloudrender.res_p*context.scene.blcloudrender.res_x),1);
# 	self.height = max(int(context.scene.blcloudrender.res_p*context.scene.blcloudrender.res_y),1);
# 	context.scene.render.resolution_x = self.width;
# 	context.scene.render.resolution_y = self.height;
# 	context.scene.render.resolution_percentage = 100;

# class ClRenderProperties(bpy.types.PropertyGroup):
# 	res_x = IntProperty(name="Res.X",default=1920,min=1,description="Image width in pixels",update=SetResolution);
# 	res_y = IntProperty(name="Res.Y",default=1080,min=1,description="Image height in pixels",update=SetResolution);
# 	res_p = FloatProperty(name="Res.%",default=0.5,min=0.01,max=1,description="Resolution percentage",update=SetResolution);
#
# 	def draw(self, context, layout):
# 		layout.row().label("Dimensions:");
#
# 		s = layout.split();
# 		c = s.column();
# 		c.row().prop(context.scene.render,"resolution_x");
# 		c.row().prop(context.scene.render,"resolution_y");
#
# 		c = s.column();
# 		c.row().prop(context.scene.render,"resolution_percentage");

class ClRenderPanel(bpy.types.Panel):
	bl_idname = "ClRenderPanel";
	bl_label = "Render";
	bl_space_type = "PROPERTIES";
	bl_region_type = "WINDOW";
	bl_context = "render";

	@classmethod
	def poll(cls, context):
		return context.scene.render.engine == config.dre_engineid;

	def draw(self, context):
		#context.scene.blcloudrender.draw(context,self.layout);
		self.layout.row().label("Dimensions:");

		s = self.layout.split();
		c = s.column();
		c.row().prop(context.scene.render,"resolution_x");
		c.row().prop(context.scene.render,"resolution_y");

		c = s.column();
		c.row().prop(context.scene.render,"resolution_percentage");

class ClSamplingProperties(bpy.types.PropertyGroup):
	samples = IntProperty(name="Render",default=1000,min=1,description="Number of samples to be taken for each pixel.");
	scatterevs = IntProperty(name="Scattering",default=20,min=0,max=32,description="Maximum volume scattering events.");
	msigmas = FloatProperty(name="Sigma.S",default=18.0,min=0.001,description="Macroscopic scattering cross section for maximum density.");
	msigmaa = FloatProperty(name="Sigma.A",default=0.001,min=0.001,description="Macroscopic absorption cross section for maximum density.");
	phasef = EnumProperty(name="Phase function",default="M",items=(
		("H","Henyey-Greenstein","Henyey-Greenstein phase function. A fast approximation with plausible results."),
		("M","Mie","Precomputed RGB Mie phase function for typical cloud droplets. Being the most accurate this is also the most inefficient due to partly unvectorized table lookups. Note that spectral rendering is required to correctly sample for different wavelengths, although in case of Mie the dispersion is small enough to be approximated without separating the RGB channels.")));
	phasea = FloatProperty(name="Anisotropy",default=0.75,description="Anisotropy parameter 'g' for the Henyey-Greenstein phase function.");

	def draw(self, context, layout):
		s = layout.split();
		c = s.column();
		c.row().label("Samples:");
		c.row().prop(self,"samples");
		#seed, default 1000
		c.row().label("Light transport:");
		c.row().prop(self,"msigmas");
		c.row().prop(self,"msigmaa");

		c = s.column();
		c.row().label("Path tracing:");
		c.row().prop(self,"scatterevs");

		c.row().label("Phase function:");
		c.row().prop(self,"phasef");
		if self.phasef == "H":
			c.row().prop(self,"phasea");

class ClSamplingPanel(bpy.types.Panel):
	bl_idname = "ClSamplingPanel";
	bl_label = "Sampling";
	bl_space_type = "PROPERTIES";
	bl_region_type = "WINDOW";
	bl_context = "render";

	@classmethod
	def poll(cls, context):
		return context.scene.render.engine == config.dre_engineid;

	def draw(self, context):
		context.scene.blcloudsampling.draw(context,self.layout);

class ClGridProperties(bpy.types.PropertyGroup):
	detailsize = FloatProperty(name="Detail size",default=0.01,min=0.0001,precision=4,description="Smallest detail size during scene costruction in blender units.");
	maxdepth = IntProperty(name="Max Depth",default=12,min=1,description="Maximum octree depth. Limiting depth to smaller values increases render performance, but at the cost of less sparse data and higher memory requirements.");
	qfbandw = FloatProperty(name="Band",default=1.0,min=0.01,precision=2,description="Outer narrow-band width of the low-resolution distance query field. This field is only constructed when the 'distance' output of the SceneInfo-node is used. A separate low-resolution field is created to allow approximate distance evaluation in larger global domains, as opposed to tight and local surface-surrounding field of the high-resolution field.");

	def draw(self, context, layout):
		s = layout.split();
		c = s.column();
		c.row().label("Resolution:");
		c.row().prop(self,"detailsize");
		c.row().label("Octree:");
		c.row().prop(self,"maxdepth");

		c = s.column();
		c.row().label("SceneInfo Query:");
		c.row().prop(self,"qfbandw");

class ClGridPanel(bpy.types.Panel):
	bl_idname = "ClGridPanel";
	bl_label = "Grid";
	bl_space_type = "PROPERTIES";
	bl_region_type = "WINDOW";
	bl_context = "render";

	@classmethod
	def poll(cls, context):
		return context.scene.render.engine == config.dre_engineid;

	def draw(self, context):
		context.scene.blcloudgrid.draw(context,self.layout);

class ClPerformanceProperties(bpy.types.PropertyGroup):
	tilex = IntProperty(name="X",default=128,min=4,description="Horizontal tile size. By design all threads contribute to one tile simultaneously."); #step=2
	tiley = IntProperty(name="Y",default=128,description="Vertical tile size. By design all threads contribute to one tile simultaneously.");
	#deprecated ######
	cache = BoolProperty(name="Enable",default=False,description="Enable the grid disk caching for individual objects. Until the object cache is reconstructed, the object is unaffected by any changes to it or its nodes.");
	cachelayer = IntProperty(name="Layer",default=10,min=0,max=19,description="Objects in this scene layer are read from the cache, or written to it if the cache doesn't exist. Remove the object from this layer to reconstruct the cache, or manually delete the cache files.");
	cachedir = StringProperty(name="Path",subtype="DIR_PATH",default="/tmp/",description="Location for the VDB cache.");
	samples = IntProperty(name="Int.Samples",default=100,min=1,description="Maximum number of samples taken internally by the render engine before returning to update the render result. Higher number of internal samples results in slightly faster render times, but also increases the interval between visual updates.");

	def draw(self, context, layout):
		s = layout.split();
		c = s.column();
		c.row().label("Tiles:");
		c.row().prop(self,"tilex");
		c.row().prop(self,"tiley");
		c.row().label("Internal sampling:");
		c.row().prop(self,"samples");

		c = s.column();
		c.row().label("Caching:");#,icon="FILE");
		c.row().prop(self,"cache");
		c.row().prop(self,"cachelayer");
		c.row().label("Location:");
		c.row().prop(self,"cachedir");

class ClPerformancePanel(bpy.types.Panel):
	bl_idname = "ClPerformancePanel";
	bl_label = "Performance";
	bl_space_type = "PROPERTIES";
	bl_region_type = "WINDOW";
	bl_context = "render";

	@classmethod
	def poll(cls, context):
		return context.scene.render.engine == config.dre_engineid;

	def draw(self, context):
		context.scene.blcloudperf.draw(context,self.layout);

class ClLayerPanel(bpy.types.Panel):
	bl_idname = "ClLayerPanel";
	bl_label = "Layer";
	bl_space_type = "PROPERTIES";
	bl_region_type = "WINDOW";
	bl_context = "render_layer";

	@classmethod
	def poll(cls, context):
		return context.scene.render.engine == config.dre_engineid;

	def draw(self, context):
		#self.layout.row().prop(context.scene.blcloudrender,"layer");
		#pass #context.scene.blcloudperf.draw(context,self.layout);
		s = self.layout.split();
		s.column().prop(context.scene,"layers",text="Scene");
		s.column().prop(context.scene.render.layers.active,"layers",text="Layer");

# #proxy to enable render results outputs and rename them
# def EnableDirectional(self, context):
# 	context.scene.render.layers.active.use_pass_transmission_direct = context.scene.blcloudpasses.dirlight;
# 	for i in context.scene.node_tree.nodes:
# 		if i.type != "R_LAYERS":
# 			continue;
# 		n = i.outputs.find("TransDir");
# 		#i.outputs[n].name = "Directional";
# 		i.outputs[n].enabled = context.scene.blcloudpasses.dirlight;
#
# def EnableEnvironment(self, context):
# 	context.scene.render.layers.active.use_pass_transmission_indirect = context.scene.blcloudpasses.envlight;
# 	for i in context.scene.node_tree.nodes:
# 		if i.type != "R_LAYERS":
# 			continue;
# 		n = i.outputs.find("TransInd");
# 		#i.outputs[n].name = "Environment";
# 		i.outputs[n].enabled = context.scene.blcloudpasses.envlight;
#
# class ClPassProperties(bpy.types.PropertyGroup):
# 	dirlight = BoolProperty(name="Directional",default=False,update=EnableDirectional);
# 	envlight = BoolProperty(name="Environment",default=False,update=EnableEnvironment);
#
# 	def draw(self, context, layout):
# 		layout.row().prop(self,"dirlight");
# 		layout.row().prop(self,"envlight");

class ClPassPanel(bpy.types.Panel):
	bl_idname = "ClPassPanel";
	bl_label = "Passes";
	bl_space_type = "PROPERTIES";
	bl_region_type = "WINDOW";
	bl_context = "render_layer";

	@classmethod
	def poll(cls, context):
		return context.scene.render.engine == config.dre_engineid;

	def draw(self, context):
		self.layout.row().label("Directional and environment lighting");
		self.layout.row().prop(context.scene.render.layers.active,"use_pass_combined"); #update(), render()
		self.layout.row().prop(context.scene.render.layers.active,"use_pass_transmission_direct");
		self.layout.row().prop(context.scene.render.layers.active,"use_pass_transmission_indirect");
		#context.scene.blcloudpasses.draw(context,self.layout);
		self.layout.row().label("Shadow mask for external depth buffer");
		self.layout.row().prop(context.scene.render.layers.active,"use_pass_shadow");

def TextureSelectionRGB(self, context):
	return [("(droplet.nan)","( Not Used )","Map disabled","X",0)]\
		+list(filter(lambda t: bpy.data.images[t[0]].channels == 4 and t[0] != "Render Result" and t[0] != "Viewer Node",[(m.name,m.name,m.name,"IMAGE_COL",x+1) for x, m in enumerate(bpy.data.images)])); #IMAGE_RGP

def TextureSelectionR(self, context):
	return [("(droplet.nan)","( Not Used )","Map disabled","X",0)]\
		+list(filter(lambda t: bpy.data.images[t[0]].channels == 4 and t[0] != "Render Result" and t[0] != "Viewer Node",[(m.name,m.name,m.name,"TEXTURE",x+1) for x, m in enumerate(bpy.data.images)]));

class ClEnvironmentProperties(bpy.types.PropertyGroup):
	envtex = EnumProperty(name="Environment Map",items=TextureSelectionRGB,description="Equirectangularly mapped environment texture to be used (sky and groud). The image should be low-frequency and should not include the sun. Droplet ignores the zeroth order environment lighting, so that optionally different background may be used for alpha compositing.");
	depthtex = EnumProperty(name="Depth Map",items=TextureSelectionR,description="The optionial depth texture from the primary render engine. This can be used to calculate local shadowing for the reconstructed locations (shadow pass). The camera properties should match the primary render (view, projection) and the map should be unnormalized and linear. Cycles depth output is readily usable, assuming that the camera is shared between the two renders. Filtering may also have to be disabled (Gaussian, width = ~0) for correct results.");
	depthcomp = BoolProperty(name="Depth Composition",default=False,description="Enable depth map occlusion testing. Use the depth map to limit the camera ray traversal.");
	occlusion = BoolProperty(name="Occlusion Geometry",default=False,description="Enable holdout geometry occlusion testing. Every object marked as holdout will occlude rays creating shadows and lightshafts. Holdout object itself is not visible. This feature may have a significant performance impact, and requires Droplet to be built with Intel Embree.");

	def draw(self, context, layout):
		layout.row().label("Environment Lighting:");
		layout.row().prop(self,"envtex");
		layout.row().label("Render Compositing:");
		layout.row().prop(self,"depthtex"); #TODO: allow only 1-channel textures
		layout.row().prop(self,"depthcomp");
		layout.row().prop(self,"occlusion");

class ClEnvironmentPanel(bpy.types.Panel):
	bl_idname = "ClEnvironmentPanel";
	bl_label = "Environment"
	bl_space_type = "PROPERTIES";
	bl_region_type = "WINDOW";
	bl_context = "world";

	@classmethod
	def poll(cls, context):
		return context.scene.render.engine == config.dre_engineid;

	def draw(self, context):
		context.scene.world.droplet.draw(context,self.layout);

#https://wiki.blender.org/index.php/Linking_Custom_Node_Tree_to_DataBlock
def NodeGroupSelection(self, context):
	#ml = sorted(context.scene.timeline_markers,key = lambda m: m.frame,reverse=True);
	#return [(m.name, m.name, m.name, x) for x, m in enumerate(ml)];
	#return [(m.name,m.name,m.name,"NODETREE",x) for x, m in enumerate(bpy.data.node_groups)];
	return list(filter(lambda t: bpy.data.node_groups[t[0]].bl_idname == "ClNodeTree",[(m.name,m.name,m.name,"NODETREE",x) for x, m in enumerate(bpy.data.node_groups)]));

class ClObjectProperties(bpy.types.PropertyGroup):
	holdout = BoolProperty(name="Holdout Mesh",default=False,description="Tell Droplet that this is a holdout mesh. Holdouts will occlude rays and create shadowing among clouds. This is also required when compositing with results from other render engines. Available only if \"occlusion geometry\" option is enabled and Droplet was built with Intel Embree support.");
	nodetree = EnumProperty(name="Node group",items=NodeGroupSelection,description="Node group to be used for this object");
	#cache
	#vdbdir = StringProperty(name="Directory",subtype="DIR_PATH",default="/tmp/");
	# vdbsdf = BoolProperty(name="Surface",default=False,description="Enable surface cache for this object. Changes to the scene or node setup won't affect the result, until the cache is recomputed.");
	# vdbfog = BoolProperty(name="Fog",default=False,description="Enable fog cache for this object. Changes to the scene or node setup won't affect the result, until the cache is recomputed.");
	#smoke
	vdbcache = StringProperty(name="File",subtype="FILE_PATH",description="Path to the OpenVDB .vdb cache. Required if the node tree makes use of the SmokeCache. Can be set to point to the Blender produced .vdb cache of desired frame (smoke simulations), for example. Loaded density and/or velocity grids will be upsampled to match current grid resolution.");
	vdbrho = StringProperty(name="Density",default="density",description="Density grid name. For Blender smoke caches, default value can be used. Leave empty if unavailable.");
	vdbvel = StringProperty(name="Velocity",default="velocity",description="Velocity grid name. For Blender smoke caches, default value can be used. Leave empty if unavailable.");

class ClMaterialPanel(bpy.types.Panel):
	bl_idname = "ClMaterialPanel";
	bl_label = "Surface";
	bl_space_type = "PROPERTIES";
	bl_region_type = "WINDOW";
	bl_context = "material";
	bl_options = {'HIDE_HEADER'};

	@classmethod
	def poll(cls, context):
		return context.scene.render.engine == config.dre_engineid and context.active_object.type == 'MESH';

	def draw(self, context):
		self.layout.row().prop(context.object.droplet,"holdout");
		if not context.object.droplet.holdout:
			self.layout.row().prop(context.object.droplet,"nodetree");

# class ClCachePanel(bpy.types.Panel):
# 	bl_idname = "ClCachePanel";
# 	bl_label = "OpenVDB disk cache";
# 	bl_space_type = "PROPERTIES";
# 	bl_region_type = "WINDOW";
# 	bl_context = "material";
#
# 	@classmethod
# 	def poll(cls, context):
# 		return context.scene.render.engine == config.dre_engineid and context.active_object.type == 'MESH';
#
# 	def draw(self, context):
# 		#TODO: cache entries should be grid resolution dependent
# 		#self.layout.row().prop(context.object.droplet,"vdbdir");
# 		self.layout.row().prop(context.object.droplet,"vdbsdf");
# 		self.layout.row().prop(context.object.droplet,"vdbfog");

class ClSmokePanel(bpy.types.Panel):
	bl_idname = "ClSmokePanel";
	bl_label = "OpenVDB smoke cache";
	bl_space_type = "PROPERTIES";
	bl_region_type = "WINDOW";
	bl_context = "material";

	@classmethod
	def poll(cls, context):
		return context.scene.render.engine == config.dre_engineid and context.active_object.type == 'MESH';

	def draw(self, context):
		self.layout.row().label("VDB source for the SmokeCache node");
		self.layout.row().prop(context.object.droplet,"vdbcache");
		self.layout.row().prop(context.object.droplet,"vdbrho");
		self.layout.row().prop(context.object.droplet,"vdbvel");


class ClParticleSystemProperties(bpy.types.PropertyGroup):
	nodetree = EnumProperty(name="Node group",items=NodeGroupSelection,description="Node group to be used for this particle system");

	def draw(self, context, layout):
		layout.row().prop(self,"nodetree");

class ClParticleSystemPanel(bpy.types.Panel):
	bl_idname = "ClParticleSystemPanel";
	bl_label = "Render";
	bl_space_type = "PROPERTIES";
	bl_region_type = "WINDOW";
	bl_context = "particle";

	@classmethod
	def poll(cls, context):
		return context.scene.render.engine == config.dre_engineid and context.particle_system != None;

	def draw(self, context):
		context.particle_system.settings.droplet.draw(context,self.layout);

class ClLampProperties(bpy.types.PropertyGroup):
	intensity = FloatProperty(name="Intensity",default=1.0,min=0.0);
	color = FloatVectorProperty(name="Color",default=[1,1,1],subtype='COLOR',size=3);
	angle = FloatProperty(name="Angle",default=0.010,min=0.0,max=1.0,precision=3);

	def draw(self, context, layout):
		s = layout.split();
		c = s.column();
		c.row().prop(self,"intensity");
		c.row().prop(self,"angle");

		c = s.column();
		c.row().prop(self,"color");

class ClLampPanel(bpy.types.Panel):
	bl_idname = "ClLampPanel";
	bl_label = "Lamp";
	bl_space_type = "PROPERTIES";
	bl_region_type = "WINDOW";
	bl_context = "data";

	@classmethod
	def poll(cls, context):
		return context.scene.render.engine == config.dre_engineid and context.active_object.type == 'LAMP';

	def draw(self, context):
		context.object.data.droplet.draw(context,self.layout);
