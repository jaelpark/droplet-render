
import bpy

from nodeitems_utils import NodeItem, NodeCategory
from bpy.props import BoolProperty, IntProperty, FloatProperty, EnumProperty, FloatVectorProperty, PointerProperty

from . import config

#https://wiki.blender.org/index.php/Dev:Py/Scripts/Guidelines/Layouts

###############################################################

class ClNodeTree(bpy.types.NodeTree):
	bl_idname = "ClNodeTree";
	#bl_label = dre_name+" nodes";
	bl_label = "Droplet Render nodes";
	bl_icon = "MATSPHERE"; #MATERIAL_DATA
	#bl_description = dre_name+" nodes";
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

	value = FloatProperty(name="",default=0,min=0);
	type = 'VALUE';

	def draw(self, context, layout, node, x):
		if self.is_linked or self.is_output:
			layout.label(self.name);
		else:
			layout.prop(self,"value",text=self.name);

	def draw_color(self, context, node):
		return (0.5,0.5,0.5,1);

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

#class ClNodeGridSocket(bpy.types.NodeSocket):
#	bl_idname = "ClNodeGridSocket";
#	bl_label = "Grid socket";
#
#	value = FloatProperty(name="",default=0);
#	type = 'CUSTOM';
#
#	def draw(self, context, layout, node, x):
#		layout.label(self.name);
#
#	def draw_color(self, context, node):
#		#return (0.9,0.9,0.2,1);
#		return (1.0,0.8,0.8,1);

class ClNodeSurfaceSocket(bpy.types.NodeSocket):
	bl_idname = "ClNodeSurfaceSocket";
	bl_label = "Displacement socket";

	value = FloatProperty(name="",default=0);
	type = 'CUSTOM';

	def draw(self, context, layout, node, x):
		layout.label(self.name);

	def draw_color(self, context, node):
		return (1,1,1,1);

# class ClNodeVectorFieldSocket(bpy.types.NodeSocket):
# 	bl_idname = "ClNodeVectorFieldSocket";
# 	bl_label = "Vector field socket";
#
# 	value = FloatProperty(name="",default=0);
# 	type = 'CUSTOM';
#
# 	def draw(self, context, layout, node, x):
# 		layout.label(self.name);
#
# 	def draw_color(self, context, node):
# 		return (0.6,0,1,1);

class ClNodeFieldSocket(bpy.types.NodeSocket):
	bl_idname = "ClNodeFieldSocket";
	bl_label = "Field socket";

	value = FloatProperty(name="",default=0);
	type = 'CUSTOM';

	def draw(self, context, layout, node, x):
		layout.label(self.name);

	def draw_color(self, context, node):
		return (0.6,0,1,1);

#class ClNodeFloatInput(bpy.types.Node):
#	bl_idname = "ClNodeFloatInput";
#	bl_label = "Float";

class ClNodeSurfaceInput(bpy.types.Node):
	bl_idname = "ClNodeSurfaceInput";
	bl_label = "Surface";

	def init(self, context):
		self.outputs.new("ClNodeSurfaceSocket","Surface");

class ClNodeParticleInput(bpy.types.Node):
	bl_idname = "ClNodeParticleInput";
	bl_label = "ParticleSystem";

	#bool: velocity

	def init(self, context):
		self.inputs.new("ClNodeFloatSocket","Raster.res");
		self.inputs.new("ClNodeFloatSocket","Weight");
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

#class ClPropertyEmpty(bpy.types.PropertyGroup):
	#pass

class ClNodeSurfaceOutput(bpy.types.Node):
	bl_idname = "ClNodeSurfaceOutput";
	bl_label = "Surface Output";

	#props = PointerProperty(type=ClPropertyEmpty);

	def init(self, context):
		#TODO: 3d space to surface node: using sdf gradients, get the closest surface point and convert to texc
		#self.inputs.new("ClNodeShaderSocket","Shader");
		self.inputs.new("ClNodeFogSocket","Fog");
		self.inputs.new("ClNodeSurfaceSocket","Surface.Field");
		self.inputs.new("ClNodeSurfaceSocket","Surface");

		#self.color = (0.7,0.7,0.8);
		#self.use_custom_color = True;

	#def draw_buttons(self, context, layout):
		#layout.row().label("wefwefwef",icon="LOCKED");

#class ClPropertyDisplacement(bpy.types.PropertyGroup):
#	maxd = FloatProperty(name="Max",default=1.0,min=0,description="Maximum displacement value. Inputs exceeding this will be clamped");
#
#	def draw(self, context, layout):
#		layout.row().prop(self,"maxd");

#class ClNodeFieldOutput(bpy.types.Node):
#	bl_idname = "ClNodeFieldOutput";
#	bl_label = "Field Output";
#
#	def init(self, context):
#		self.inputs.new("ClNodeFieldSocket","Field");

#class ClNodeFogVolume(bpy.types.Node):
#	bl_idname = "ClNodeFogVolume";
#	bl_label = "SDF to Fog";
#
#	def init(self, context):
#		self.inputs.new("ClNodeSurfaceSocket","Surface");
#		self.outputs.new("ClNodeFogSocket","Fog");

class ClNodeAdvection(bpy.types.Node):
	bl_idname = "ClNodeAdvection";
	bl_label = "Advection";

	def init(self, context):
		self.inputs.new("ClNodeIntSocket","Iterations");
		self.inputs.new("ClNodeFloatSocket","Step size");
		self.inputs.new("ClNodeFogSocket","Fog");
		self.outputs.new("ClNodeFogSocket","Fog");

	#def draw_buttons(self, context, layout):
		#layout.row().label("Baked advection");

class ClNodeDisplacement(bpy.types.Node):
	bl_idname = "ClNodeDisplacement";
	bl_label = "Displacement";

	#props = PointerProperty(type=ClPropertyDisplacement);
	#props = PointerProperty(type=ClPropertyEmpty);

	def init(self, context):
		self.inputs.new("ClNodeFloatSocket","Distance");
		self.inputs.new("ClNodeSurfaceSocket","Surface");
		self.outputs.new("ClNodeSurfaceSocket","Surface");

	#def draw_buttons(self, context, layout):
		#self.props.draw(context,layout);

#class ClPropertyfBmPerlinNoise(bpy.types.PropertyGroup):
#	#TODO: allow node-driven inputs
#	octaves = IntProperty(name="Octaves",default=2,min=1,max=12);
#	freq = FloatProperty(name="Frequency",default=1.0,min=0.0,max=10.0);
#	amp = FloatProperty(name="Amplitude",default=1.0,min=0.0);
#	fjump = FloatProperty(name="FJump",default=1.5,min=0.0);
#	gain = FloatProperty(name="Gain",default=0.5,min=0.0);
#	billow = FloatProperty(name="Billow",default=0.0,min=0.0,max=1.0);
#	qscale = FloatProperty(name="QScale",default=1.0,min=0.0);
#
#	def draw(self, context, layout):
#		layout.row().prop(self,"octaves");
#		layout.row().prop(self,"freq");
#		layout.row().prop(self,"amp");
#		layout.row().prop(self,"fjump");
#		layout.row().prop(self,"gain");
#		layout.row().prop(self,"billow");
#		layout.row().prop(self,"qscale");

#class ClPropertyfBmPerlinNoise(bpy.types.PropertyGroup):
#	#maxd = FloatProperty(name="Max",default=1.0,min=0,description="Maximum displacement value. Inputs exceeding this will be clamped");
#	#lckd = BoolProperty(name="",default=True); #ICON=LINKED/UNLINKED, LOCKED/UNLOCKED
#	octaves = IntProperty(name="Octaves",default=2,min=1,max=12);
#
#	def draw(self, context, layout):
#		#r = layout.row();
#		#r.prop(self,"maxd");
#		#r.prop(self,"lckd");
#		layout.row().prop(self,"octaves");

#TODO: deprecate this in favor of displacement + separate value perlin noise?
#^^this alternative is still unoptimal, as all (noise) nodes of the tree would be evaluated (above level+1), even if only was was needed for current displacement node.
class ClNodefBmPerlinNoise(bpy.types.Node):
	bl_idname = "ClNodefBmPerlinNoise";
	bl_label = "fBm Perlin";

	#props = PointerProperty(type=ClPropertyfBmPerlinNoise);

	def init(self, context):
		self.inputs.new("ClNodeIntSocket","octaves");
		self.inputs.new("ClNodeFloatSocket","freq");
		self.inputs.new("ClNodeFloatSocket","amp");
		self.inputs.new("ClNodeFloatSocket","fjump");
		self.inputs.new("ClNodeFloatSocket","gain");
		self.inputs.new("ClNodeFloatSocket","billow");
		#self.inputs.new("ClNodeFloatSocket","qscale");
		self.inputs.new("ClNodeSurfaceSocket","Surface");
		#self.inputs.new("ClNodeGridSocket","Billowing");

		self.outputs.new("ClNodeSurfaceSocket","Surface");
		#self.outputs.new("ClNodeGridSocket","Distance");

	#def draw_buttons(self, context, layout):
		#self.props.draw(context,layout);

class ClNodeCategory(NodeCategory):
	@classmethod
	def poll(cls, context):
		return context.space_data.tree_type == "ClNodeTree";

categories = [
	ClNodeCategory("INPUT_CATEGORY","Input",items = [
		NodeItem("ClNodeSurfaceInput"),
		NodeItem("ClNodeParticleInput"),
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
	]),
	ClNodeCategory("DENSITY_CATEGORY","Fog",items = [
		NodeItem("ClNodeAdvection"),
	]),
	ClNodeCategory("SURFACE_CATEGORY","Surface",items = [
		NodeItem("ClNodeDisplacement"),
		NodeItem("ClNodefBmPerlinNoise"),
	]),
];
