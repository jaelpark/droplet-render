
import bpy

from nodeitems_utils import NodeItem, NodeCategory
from bpy.props import BoolProperty, IntProperty, FloatProperty, EnumProperty, FloatVectorProperty, PointerProperty

from . import config

#https://wiki.blender.org/index.php/Dev:Py/Scripts/Guidelines/Layouts

###############################################################

class ClNodeTree(bpy.types.NodeTree):
	bl_idname = "ClNodeTree";
	bl_label = "Droplet Render nodes";
	bl_icon = "MATSPHERE"; #MATERIAL_DATA
	bl_description = bl_label;

	@classmethod
	def poll(cls, context):
		return context.scene.render.engine == config.dre_engineid;

class ClNodeIntSocket(bpy.types.NodeSocket):
	bl_idname = "ClNodeIntSocket";
	bl_label = "Int socket";

	value = IntProperty(name="",default=0,min=0);
	type = "INT";

	def draw(self, context, layout, node, x):
		if self.is_linked or self.is_output:
			layout.label(self.name);
		else:
			layout.prop(self,"value",text=self.name);

	def draw_color(self, context, node):
		#return (1.0,0.8,0.8,1);
		return (0.5,0.5,0.5,1);

class ClNodeFloatSocket(bpy.types.NodeSocket):
	bl_idname = "ClNodeFloatSocket";
	bl_label = "Float socket";

	value = FloatProperty(name="",default=0,precision=3);
	type = 'VALUE';

	def draw(self, context, layout, node, x):
		if self.is_linked or self.is_output:
			layout.label(self.name);
		else:
			layout.prop(self,"value",text=self.name);

	def draw_color(self, context, node):
		return (0.5,0.5,0.5,1);

class ClNodeVectorSocket(bpy.types.NodeSocket):
	bl_idname = "ClNodeVectorSocket";
	bl_label = "Vector socket";

	value = FloatProperty(name="",default=0); #TODO: vector property
	type = 'VALUE';

	def draw(self, context, layout, node, x):
		layout.label(self.name);

	def draw_color(self, context, node):
		return (0.3,0.3,0.3,1);

class ClNodeShaderSocket(bpy.types.NodeSocket):
	bl_idname = "ClNodeShaderSocket";
	bl_label = "Shader socket";

	value = FloatProperty(name="",default=0);
	type = 'SHADER';

	def draw(self, context, layout, node, x):
		layout.label(self.name);

	def draw_color(self, context, node):
		return (0.4,1,0.4,1);
		#return (0.9,0.9,0.2,1);

class ClNodeFogSocket(bpy.types.NodeSocket):
	bl_idname = "ClNodeFogSocket";
	bl_label = "Fog socket";

	value = FloatProperty(name="",default=0);
	type = 'CUSTOM';

	def draw(self, context, layout, node, x):
		layout.label(self.name);

	def draw_color(self, context, node):
		return (0.9,0.9,0.2,1);

class ClNodeSurfaceSocket(bpy.types.NodeSocket):
	bl_idname = "ClNodeSurfaceSocket";
	bl_label = "Displacement socket";

	value = FloatProperty(name="",default=0);
	type = 'CUSTOM';

	def draw(self, context, layout, node, x):
		layout.label(self.name);

	def draw_color(self, context, node):
		return (1,1,1,1);

class ClNodeVectorFieldSocket(bpy.types.NodeSocket):
	bl_idname = "ClNodeVectorFieldSocket";
	bl_label = "Vector field socket";

	value = FloatProperty(name="",default=0);
	type = 'CUSTOM';

	def draw(self, context, layout, node, x):
		layout.label(self.name);

	def draw_color(self, context, node):
		return (0.9,0.5,0.2,1);

class ClNodeFieldSocket(bpy.types.NodeSocket):
	bl_idname = "ClNodeFieldSocket";
	bl_label = "Field socket";

	value = FloatProperty(name="",default=0);
	type = 'CUSTOM';

	def draw(self, context, layout, node, x):
		layout.label(self.name);

	def draw_color(self, context, node):
		return (0.6,0,1,1);

class ClNodeVoxelInfo(bpy.types.Node):
	bl_idname = "ClNodeVoxelInfo";
	bl_label = "Voxel Info";

	def init(self, context):
		self.outputs.new("ClNodeVectorSocket","Vox.world");
		self.outputs.new("ClNodeVectorSocket","CPT.world");
		self.outputs.new("ClNodeFloatSocket","Distance");
		self.outputs.new("ClNodeFloatSocket","Density");

class ClNodeSceneInfo(bpy.types.Node):
	bl_idname = "ClNodeSceneInfo";
	bl_label = "Scene Info";

	def init(self, context):
		self.inputs.new("ClNodeVectorSocket","World");
		self.outputs.new("ClNodeFloatSocket","Distance");
		self.outputs.new("ClNodeFloatSocket","Density");
		self.outputs.new("ClNodeFloatSocket","Density.Final");

class ClNodeSurfaceInput(bpy.types.Node):
	bl_idname = "ClNodeSurfaceInput";
	bl_label = "Surface";

	def init(self, context):
		self.outputs.new("ClNodeSurfaceSocket","Surface");

class ClNodeParticleInput(bpy.types.Node):
	bl_idname = "ClNodeParticleInput";
	bl_label = "ParticleSystem";

	def init(self, context):
		#self.inputs.new("ClNodeFloatSocket","Raster.res");
		#self.inputs.new("ClNodeFloatSocket","Weight");
		self.inputs.new("ClNodeFloatSocket","Size");
		self.inputs.new("ClNodeFloatSocket","Cutoff");
		self.outputs.new("ClNodeFogSocket","Fog");
		#self.outputs.new("ClNodeVectorFieldSocket","Velocity");

class ClNodeSmokeCache(bpy.types.Node):
	bl_idname = "ClNodeSmokeCache";
	bl_label = "SmokeCache";

	def init(self, context):
		self.outputs.new("ClNodeFogSocket","Fog");

class ClNodeFogPostInput(bpy.types.Node):
	bl_idname = "ClNodeFogPostInput";
	bl_label = "FogPostInput";

	def init(self, context):
		 self.outputs.new("ClNodeFogSocket","Fog");

class ClNodeFloatAdd(bpy.types.Node):
	bl_idname = "ClNodeFloatAdd";
	bl_label = "Add";

	def init(self, context):
		self.inputs.new("ClNodeFloatSocket","a");
		self.inputs.new("ClNodeFloatSocket","b");
		self.outputs.new("ClNodeFloatSocket","Out");

class ClNodeFloatSub(bpy.types.Node):
	bl_idname = "ClNodeFloatSub";
	bl_label = "Subtract";

	def init(self, context):
		self.inputs.new("ClNodeFloatSocket","a");
		self.inputs.new("ClNodeFloatSocket","b");
		self.outputs.new("ClNodeFloatSocket","Out");

class ClNodeFloatMul(bpy.types.Node):
	bl_idname = "ClNodeFloatMul";
	bl_label = "Multiply";

	def init(self, context):
		self.inputs.new("ClNodeFloatSocket","a");
		self.inputs.new("ClNodeFloatSocket","b");
		self.outputs.new("ClNodeFloatSocket","Out");

class ClNodeFloatDiv(bpy.types.Node):
	bl_idname = "ClNodeFloatDiv";
	bl_label = "Divide";

	def init(self, context):
		self.inputs.new("ClNodeFloatSocket","a");
		self.inputs.new("ClNodeFloatSocket","b");
		self.outputs.new("ClNodeFloatSocket","Out");

class ClNodeFloatPow(bpy.types.Node):
	bl_idname = "ClNodeFloatPow";
	bl_label = "Power";

	def init(self, context):
		self.inputs.new("ClNodeFloatSocket","a");
		self.inputs.new("ClNodeFloatSocket","b");
		self.outputs.new("ClNodeFloatSocket","Out");

class ClNodeFloatMin(bpy.types.Node):
	bl_idname = "ClNodeFloatMin";
	bl_label = "Min";

	def init(self, context):
		self.inputs.new("ClNodeFloatSocket","a");
		self.inputs.new("ClNodeFloatSocket","b");
		self.outputs.new("ClNodeFloatSocket","Out");

class ClNodeFloatMax(bpy.types.Node):
	bl_idname = "ClNodeFloatMax";
	bl_label = "Max";

	def init(self, context):
		self.inputs.new("ClNodeFloatSocket","a");
		self.inputs.new("ClNodeFloatSocket","b");
		self.outputs.new("ClNodeFloatSocket","Out");

class ClNodeFloatInput(bpy.types.Node):
	bl_idname = "ClNodeFloatInput";
	bl_label = "Value";

	def init(self, context):
		self.inputs.new("ClNodeFloatSocket","Value");
		self.outputs.new("ClNodeFloatSocket","Out");

class ClNodeVectorInput(bpy.types.Node):
	bl_idname = "ClNodeVectorInput";
	bl_label = "Vector";

	def init(self, context):
		self.inputs.new("ClNodeFloatSocket","x");
		self.inputs.new("ClNodeFloatSocket","y");
		self.inputs.new("ClNodeFloatSocket","z");
		self.outputs.new("ClNodeVectorSocket","Out");

class ClNodeVectorMath(bpy.types.Node):
	bl_idname = "ClNodeVectorMath";
	bl_label = "Vector Math";

	op = EnumProperty(name="",default="+",items=(
		("+","Add",""),
		("-","Subtract",""),
		("*","Multiply",""),
		("/","Divide",""),
		("X","Cross","")));

	def init(self, context):
		self.inputs.new("ClNodeVectorSocket","a");
		self.inputs.new("ClNodeVectorSocket","b");
		self.outputs.new("ClNodeVectorSocket","Out");

	def draw_buttons(self, context, layout):
		layout.prop(self,"op");

#class ClPropertyEmpty(bpy.types.PropertyGroup):
	#pass

class ClNodeFbmNoise(bpy.types.Node):
	bl_idname = "ClNodeFbmNoise";
	bl_label = "FbmNoise";

	def init(self, context):
		self.inputs.new("ClNodeIntSocket","octaves");
		self.inputs.new("ClNodeFloatSocket","freq");
		self.inputs.new("ClNodeFloatSocket","amp");
		self.inputs.new("ClNodeFloatSocket","fjump");
		self.inputs.new("ClNodeFloatSocket","gain");
		self.inputs.new("ClNodeVectorSocket","World");
		self.outputs.new("ClNodeFloatSocket","Out.Scalar");
		self.outputs.new("ClNodeVectorSocket","Out.Vector");
		self.outputs.new("ClNodeFloatSocket","Max");

class ClNodeSurfaceOutput(bpy.types.Node):
	bl_idname = "ClNodeSurfaceOutput";
	bl_label = "Surface Output";

	def init(self, context):
		#TODO: 3d space to surface node: using sdf gradients, get the closest surface point and convert to texc
		#self.inputs.new("ClNodeShaderSocket","Shader");
		self.inputs.new("ClNodeFogSocket","Fog");
		self.inputs.new("ClNodeFogSocket","Fog.Post");
		self.inputs.new("ClNodeSurfaceSocket","Surface");

		#self.color = (0.7,0.7,0.8);
		#self.use_custom_color = True;

class ClNodeComposite(bpy.types.Node):
	bl_idname = "ClNodeComposite";
	bl_label = "Composition";

	def init(self, context):
		self.inputs.new("ClNodeFloatSocket","In");
		self.inputs.new("ClNodeFogSocket","Fog");
		self.outputs.new("ClNodeFogSocket","Fog");

class ClNodeAdvection(bpy.types.Node):
	bl_idname = "ClNodeAdvection";
	bl_label = "Advection";

	#full_iter = BoolProperty(name="Full Iteration",default=False,description="Always step the specified distance regardless whether density above the threshold was encountered before reaching the iteration limit. Enabling this will in some cases form bodies separated from the advection source while additionally also creating a better density variation when advecting from a heterogenous for volume.");
	break_iter = BoolProperty(name="Break Iteration",default=False,description="Break the iteration when encountering a density above specificied threshold. Density variation will be more subtle and the advection will no longer form bodies separated from the advection source. Advection performance is significantly improved for low threshold values and dense sources.");

	sample_local = BoolProperty(name="Provide Local Data",default=False,description="When enabled, local density information is provided for the VoxelInfo node (for each advection step when not using SceneInfo during post-processing). For performance reasons and by the assumption that advection is performed during post-processing, this is turned off by default. By disabling local sampling, the VoxelInfo density output will remain constant (the value for the voxel currently being processed), which may also be desired in some cases.");

	def init(self, context):
		self.inputs.new("ClNodeFloatSocket","Threshold");
		self.inputs.new("ClNodeFloatSocket","Distance");
		self.inputs.new("ClNodeIntSocket","Iterations");
		#self.inputs.new("ClNodeVectorFieldSocket","Velocity");
		self.inputs.new("ClNodeFloatSocket","Density");
		self.inputs.new("ClNodeVectorSocket","Velocity");
		self.inputs.new("ClNodeFogSocket","Fog");
		self.outputs.new("ClNodeFogSocket","Fog");

	def draw_buttons(self, context, layout):
		layout.prop(self,"break_iter");
		layout.prop(self,"sample_local");

class ClNodeDisplacement(bpy.types.Node):
	bl_idname = "ClNodeDisplacement";
	bl_label = "Displacement";

	#props = PointerProperty(type=ClPropertyDisplacement);
	#props = PointerProperty(type=ClPropertyEmpty);

	def init(self, context):
		self.inputs.new("ClNodeFloatSocket","Distance");
		self.inputs.new("ClNodeFloatSocket","Max");
		self.inputs.new("ClNodeFloatSocket","Billow");
		self.inputs.new("ClNodeSurfaceSocket","Surface");
		self.outputs.new("ClNodeSurfaceSocket","Surface");

#^^this alternative is still unoptimal, as all (noise) nodes of the tree would be evaluated (above level+1), even if only was was needed for current displacement node.
# class ClNodefBmPerlinNoise(bpy.types.Node):
# 	bl_idname = "ClNodefBmPerlinNoise";
# 	bl_label = "fBm Perlin";
#
# 	#props = PointerProperty(type=ClPropertyfBmPerlinNoise);
#
# 	def init(self, context):
# 		self.inputs.new("ClNodeIntSocket","octaves");
# 		self.inputs.new("ClNodeFloatSocket","freq");
# 		self.inputs.new("ClNodeFloatSocket","amp");
# 		self.inputs.new("ClNodeFloatSocket","fjump");
# 		self.inputs.new("ClNodeFloatSocket","gain");
# 		self.inputs.new("ClNodeFloatSocket","billow");
# 		self.inputs.new("ClNodeSurfaceSocket","Surface");
# 		self.outputs.new("ClNodeSurfaceSocket","Surface");
# 		#self.outputs.new("ClNodeGridSocket","Distance");
#
# 	#def draw_buttons(self, context, layout):
# 		#self.props.draw(context,layout);

class ClNodeVectorFieldSampler(bpy.types.Node):
	bl_idname = "ClNodeVectorFieldSampler";
	bl_label = "VectorField sampler";

	def init(self, context):
		self.inputs.new("ClNodeVectorFieldSocket","Field");
		self.inputs.new("ClNodeVectorSocket","World");
		self.outputs.new("ClNodeVectorSocket","Out");

class ClNodeCategory(NodeCategory):
	@classmethod
	def poll(cls, context):
		return context.space_data.tree_type == "ClNodeTree";

categories = [
	ClNodeCategory("INPUT_CATEGORY","Input",items = [
		NodeItem("ClNodeSurfaceInput"),
		NodeItem("ClNodeParticleInput"),
		NodeItem("ClNodeSmokeCache"),
		NodeItem("ClNodeFogPostInput"),
		NodeItem("ClNodeFloatInput"),
		NodeItem("ClNodeVectorInput"),
		NodeItem("ClNodeVoxelInfo"),
		NodeItem("ClNodeSceneInfo"),
	]),
	ClNodeCategory("OUTPUT_CATEGORY","Output",items = [
		NodeItem("ClNodeSurfaceOutput"),
		#NodeItem("ClNodeFieldOutput"),
	]),
	ClNodeCategory("MATH_CATEGORY","Math",items = [
		NodeItem("ClNodeFloatAdd"),
		NodeItem("ClNodeFloatSub"),
		NodeItem("ClNodeFloatMul"),
		NodeItem("ClNodeFloatDiv"),
		NodeItem("ClNodeFloatPow"),
		NodeItem("ClNodeFloatMin"),
		NodeItem("ClNodeFloatMax"),
		NodeItem("ClNodeVectorMath"),
	]),
	ClNodeCategory("CONVERSION_CATEGORY","Conversion",items = [
		NodeItem("ClNodeVectorFieldSampler"),
	]),
	ClNodeCategory("NOISE_CATEGORY","Noise",items = [
		NodeItem("ClNodeFbmNoise"),
	]),
	ClNodeCategory("DENSITY_CATEGORY","Fog",items = [
		NodeItem("ClNodeComposite"),
		NodeItem("ClNodeAdvection"),
	]),
	ClNodeCategory("SURFACE_CATEGORY","Surface",items = [
		NodeItem("ClNodeDisplacement"),
		#NodeItem("ClNodefBmPerlinNoise"),
	]),
];
