
import bpy

from bpy.props import BoolProperty, IntProperty, FloatProperty, EnumProperty, PointerProperty, FloatVectorProperty

from . import config

class ClRenderProperties(bpy.types.PropertyGroup):
	res_x = IntProperty(name="Res.X",default=1920,min=1,description="Image width in pixels");
	res_y = IntProperty(name="Res.Y",default=1080,min=1,description="Image height in pixels");
	res_p = FloatProperty(name="Res.%",default=0.5,min=0.01,max=1,description="Resolution percentage");
	transparent = BoolProperty(name="Transparent",default=False,description="Enable alpha channel and ignore background for 0th order scattering.");
	#nodetree = EnumProperty(name="Node tree",items=NodeTreeSelection,description="Node tree to be used globally");

	def draw(self, context, layout):
		layout.row().label("Dimensions:");

		s = layout.split();
		c = s.column();
		c.row().prop(self,"res_x");
		c.row().prop(self,"res_y");

		c = s.column();
		c.row().prop(self,"res_p");
		c.row().prop(self,"transparent");

		#layout.row().prop(self,"nodetree");

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
		context.scene.blcloudrender.draw(context,self.layout);

class ClSamplingProperties(bpy.types.PropertyGroup):
	samples = IntProperty(name="Render",default=500,min=1,description="Number of samples to be taken");
	scatterevs = IntProperty(name="Scattering",default=12,min=0,max=32,description="Maximum volume scattering events.");
	msigmas = FloatProperty(name="Sigma.S",default=7.11,min=0.0,description="Macroscopic scattering cross section for maximum density.");
	msigmaa = FloatProperty(name="Sigma.A",default=0.03,min=0.0,description="Macroscopic absorption cross section for maximum density.");
	phasef = EnumProperty(name="Phase function",default="M",items=(
		("H","Henyey-Greenstein","Henyey-Greenstein phase function with anisotropy g=0.35. A fast approximation with plausible results."),
		("M","Mie","Precomputed RGB Mie phase function for typical cloud droplets. Being the most accurate phase function this is also the most inefficient due to unvectorized table lookups. Note that spectral rendering is required to correctly sample for different wavelengths.")));

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
	detailsize = FloatProperty(name="Detail size",default=0.02,min=0.0001,precision=4,description="Smallest detail size in blender units.");
	adaptive = BoolProperty(name="Adaptive resolution",default=False,description="Distance and FOV based resolution level of detail");

	def draw(self, context, layout):
		s = layout.split();
		c = s.column();
		c.row().prop(self,"detailsize");

		c = s.column();
		c.row().prop(self,"adaptive");
		#layout.row().prop(self,"voxelsize");
		#TODO: calculate and show resolution and octree depth

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

#def ClPerformanceCacheTick(self, context):
#	#TODO: need to set also for new nodes
#	for n in bpy.data.node_groups[context.scene.blcloudrender.nodetree].nodes:
#		if n.outputs and n.outputs[0].bl_idname == "ClNodeSurfaceSocket":
#			n.color = (1.0,0.7,0.0);
#			n.use_custom_color = self.cache;

class ClPerformanceProperties(bpy.types.PropertyGroup):
	tilex = IntProperty(name="X",default=128,min=4,description="Horizontal tile size. By design all threads contribute to one tile simultaneously."); #step=2
	tiley = IntProperty(name="Y",default=128,description="Vertical tile size. By design all threads contribute to one tile simultaneously.");
	cache = EnumProperty(name="Cache mode",default="0",items=(
		("0","Off","Caching disabled. Grid is always recreated."),
		("R","RW","Read the grid from a cache. A new cache is created from the current scene if unavailable. Currently manual cache management is required since it's not possible to track all changes made to the scene (data, surface nodes, textures etc.)"),
		("W","W","Write always and overwrite any previous caches. Manual overwriting is required when changes have been made.")));
	samples = IntProperty(name="Int.Samples",default=100,min=1,description="Maximum number of samples taken internally by the render engine before returning to update the render result. Higher number of internal samples results in slightly faster render times, but also increases the interval between visual updates and may limit application responsiveness.");

	def draw(self, context, layout):
		s = layout.split();
		c = s.column();
		c.row().label("Tiles:");
		c.row().prop(self,"tilex");
		c.row().prop(self,"tiley");
		c.row().label("Internal sampling:");
		c.row().prop(self,"samples");

		c = s.column();
		c.row().label("Caching:",icon="FILE");
		c.row().prop(self,"cache",expand=True);

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

#class ClRenderPanel(bpy.types.Panel):
#	bl_idname = "ClRenderPanel";
#	bl_label = "Render";
#	bl_space_type = "PROPERTIES";
#	bl_region_type = "WINDOW";
#	bl_context = "render";
#
#	@classmethod
#	def poll(cls, context):
#		return context.scene.render.engine == config.dre_engineid;
#
#	#def draw(

def NodeGroupSelection(self, context):
	#ml = sorted(context.scene.timeline_markers,key = lambda m: m.frame,reverse=True);
	#return [(m.name, m.name, m.name, x) for x, m in enumerate(ml)];
	return [(m.name,m.name,m.name,"NODETREE",x) for x, m in enumerate(bpy.data.node_groups)];

#class ClMaterialProperties(bpy.types.PropertyGroup):
	#nodetree = EnumProperty(name="Node tree",items=NodeTreeSelection,description="Node tree to be used globally");
	#
	#def draw(self, context, layout):
		#layout.row().prop(self,"nodetree");

class ClObjectProperties(bpy.types.PropertyGroup):
	nodetree = EnumProperty(name="Node group",items=NodeGroupSelection,description="Node group to be used for this object");

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
		#self.layout.row().prop(context.object,"active_material_index");
		#context.active_object.data.droplet.draw(context,self.layout);
		#context.object.data.droplet.draw(context,self.layout);
		self.layout.row().prop(context.object.droplet,"nodetree");

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
		#self.layout.row().prop(context.particle_system.droplet,"nodetree");
		context.particle_system.settings.droplet.draw(context,self.layout);
		#context.particle_system.droplet[1]['type'].draw(context,self.layout); #wtf is this tuple nonsense

def TextureSelection(self, context):
	return [(m.name,m.name,m.name,"TEXTURE",x) for x, m in enumerate(bpy.data.images)];

class ClLampProperties(bpy.types.PropertyGroup):
	intensity = FloatProperty(name="Intensity",default=50.0,min=0.0);
	color = FloatVectorProperty(name="Color",default=[1,1,1],subtype='COLOR',size=3);
	angle = FloatProperty(name="Angle",default=0.95,min=0.0,max=1.0,precision=3);
	#color = FloatVectorProperty(name="Color",default=[1,1,1]);
	#cr = FloatProperty(name="Red",default=1.0,min=0.0,max=1.0);
	#phasef = EnumProperty(name="Phase function",default="H",items=(
	#	("H","Henyey-Greenstein","Henyey-Greenstein phase function"),
	#	("M","Mie","Precomputed Mie phase function for typical cloud droplets")));
		#("c","Curve","Custom curve created with the curve editor [0,pi]"),
	#phasetex = EnumProperty(name="Phase Texture",items=TextureSelection,description="A four-channel texture describing a custom phase function (normalized probability distribution function). Texture resolution must be 1024x1 which is the mapped to [0,pi].");
	#phasesax = BoolProperty(name="Spectral approximation",default=False,description="The phase function evaluates all color channels of the texture instead of just the first one. However, since spectral rendering isn't explicitly supported, light paths are still sampled using only the first channel. This is done by the assumption that the PDF variation between channels is small.");
	#hgparam = FloatProperty(name="Anisotropy",min=0.001,max=1,default=0.35);

	def draw(self, context, layout):
		s = layout.split();
		c = s.column();
		c.row().prop(self,"intensity");
		c.row().prop(self,"angle");

		c = s.column();
		c.row().prop(self,"color");

	#def draw2(self, context, layout):
		#layout.row().prop(self,"phasef",expand=False);
		#if self.phasefunc == 'H':
			#layout.row().prop(self,"hgparam");

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

#class ClPhaseProperties(bpy.types.PropertyGroup):

# class ClPhasePanel(bpy.types.Panel):
# 	bl_idname = "ClPhasePanel";
# 	bl_label = "Phase";
# 	bl_space_type = "PROPERTIES";
# 	bl_region_type = "WINDOW";
# 	bl_context = "data";
#
# 	@classmethod
# 	def poll(cls, context):
# 		return context.scene.render.engine == config.dre_engineid and context.active_object.type == 'LAMP';
#
# 	def draw(self, context):
# 		context.object.data.droplet.draw2(context,self.layout);
