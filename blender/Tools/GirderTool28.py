
import bpy, bmesh

#print( "============= SCTIPT START ============ \n\n" )

def getSomeOrtho(fw):
    up=fw*0
    if abs(fw.x) < abs(fw.y):
        if abs(fw.z) < abs(fw.x):
            up.z=1
        else:
            up.x=1
    else:
        if abs(fw.z) < abs(fw.y):
            up.z=1
        else:
            up.y=1
    up = up - fw*up.dot(fw)
    return up.normalized()    

def makeGirder( bm, v1, v2, up=None, dl=0.4, w=0.05 ):
    v12 = v2 - v1 
    l   = v12.length
    fw = v12.normalized()
    if up is None:
        up = getSomeOrtho(fw)
    lf = up.cross(v12).normalized()
    n   = int( l/dl )
    dl  = l/n
    print(l,n)
    nv = len(bm.verts)
    print( "nverts ", nv)
    print( "fw ",     fw)
    print( "up ",     up)
    print( "lf ",     lf)
    ovs = None
    for i in range( n ):
        vi1 = v1+up*w + fw*(dl* i     );  vi1= bm.verts.new( vi1 )
        vi2 = v1+lf*w + fw*(dl*(i+0.5));  vi2= bm.verts.new( vi2 )
        vi3 = v1-up*w + fw*(dl* i     );  vi3= bm.verts.new( vi3 )
        vi4 = v1-lf*w + fw*(dl*(i+0.5));  vi4= bm.verts.new( vi4 )
        nv+=4;
        
        bm.edges.new( (vi1, vi3) )
        bm.edges.new( (vi2, vi4) )
        
        bm.edges.new( (vi1, vi2) )
        bm.edges.new( (vi2, vi3) )
        bm.edges.new( (vi3, vi4) )
        bm.edges.new( (vi4, vi1) )
        if ovs is not None:
            bm.edges.new( (vi1, ovs[0]) )
            bm.edges.new( (vi2, ovs[1]) )
            bm.edges.new( (vi3, ovs[2]) )
            bm.edges.new( (vi4, ovs[3]) )
            
            bm.edges.new( (vi1, ovs[1]) )
            bm.edges.new( (vi3, ovs[1]) )
            bm.edges.new( (vi1, ovs[3]) )
            bm.edges.new( (vi3, ovs[3]) )
        ovs = (vi1,vi2,vi3,vi4)
    vi1 = v1+up*w + fw*(dl*n);  vi1= bm.verts.new( vi1 )
    vi3 = v1-up*w + fw*(dl*n);  vi3= bm.verts.new( vi3 )
    bm.edges.new( (vi1, vi3) )
    bm.edges.new( (vi1, ovs[0]) )
    bm.edges.new( (vi3, ovs[2]) )
    bm.edges.new( (vi1, ovs[1]) )
    bm.edges.new( (vi3, ovs[1]) )
    bm.edges.new( (vi1, ovs[3]) )
    bm.edges.new( (vi3, ovs[3]) )


def makeTube( bm, v1, v2, up=None, w=0.01 ):
    v12 = v2 - v1 
    l   = v12.length
    fw = v12.normalized()
    if up is None:
        up = getSomeOrtho(fw)
    lf = up.cross(v12).normalized()

    #nv = len(bm.verts)
    #print( "nverts ", nv)
    print( "fw ", fw)
    print( "up ", up)
    print( "lf ", lf)
    #ovs = None
    
    v1a = v1 + up*w
    v2a = v1 + lf*w
    v3a = v1 - up*w
    v4a = v1 - lf*w 
    v1b=v1a+v12
    v2b=v2a+v12
    v3b=v3a+v12
    v4b=v4a+v12
    
    v1a= bm.verts.new( v1a ); v1b= bm.verts.new( v1b );
    v2a= bm.verts.new( v2a ); v2b= bm.verts.new( v2b );
    v3a= bm.verts.new( v3a ); v3b= bm.verts.new( v3b );
    v4a= bm.verts.new( v4a ); v4b= bm.verts.new( v4b );
    
    bm.faces.new( (v1a, v2a, v2b, v1b) )
    bm.faces.new( (v2a, v3a, v3b, v2b) )
    bm.faces.new( (v3a, v4a, v4b, v3b) )
    bm.faces.new( (v4a, v1a, v1b, v4b) )  

def edges2tubes( w=0.01 ):
    obj_from = bpy.context.object
    me = bpy.data.meshes .new('TubeMesh') 
    ob = bpy.data.objects.new('TubeObject', me)   
    ob.show_name = True
    bpy.context.collection.objects.link(ob)
    bm = bmesh.new()
    bm.from_mesh(me)
    ovs = obj_from.data.vertices
    print( ovs )
    for e in obj_from.data.edges:
        vs = e.vertices
        makeTube( bm, ovs[vs[0]].co, ovs[vs[1]].co, w=w )
    bm.to_mesh(me)

def edges2girders( dl=0.4, w=0.05 ):
    obj_from = bpy.context.object
    me = bpy.data.meshes .new('GirderMesh') 
    ob = bpy.data.objects.new('GirderObject', me)   
    ob.show_name = True
    bpy.context.collection.objects.link(ob)
    bm = bmesh.new()
    bm.from_mesh(me)
    ovs = obj_from.data.vertices
    print( ovs )
    for e in obj_from.data.edges:
        vs = e.vertices
        makeGirder( bm, ovs[vs[0]].co, ovs[vs[1]].co, dl=dl, w=w )
    bm.to_mesh(me)

bl_info = {
    "name": "Prokop Girder Tool",
    "description": "",
    "author": "Prokop Hapala",
    "version": (0, 0, 3),
    "blender": (2, 82, 0),
    "location": "3D View > Tools",
    "warning": "", # used for warning icon and text in addons panel
    "wiki_url": "",
    "tracker_url": "",
    "category": "Development"
}

from bpy.props import (StringProperty,
                       BoolProperty,
                       IntProperty,
                       FloatProperty,
                       FloatVectorProperty,
                       EnumProperty,
                       PointerProperty,
                       )
from bpy.types import (Panel,
                       Menu,
                       Operator,
                       PropertyGroup,
                       )


# ------------------------------------------------------------------------
#    Scene Properties
# ------------------------------------------------------------------------

class GirderProperties(PropertyGroup):

    fStep: FloatProperty(
        name = "step length",
        description = "length of girder element",
        default = 0.2,
        min = 0.01,
        max = 1.0
        )

    fWidth: FloatProperty(
        name = "width",
        description = "girder width",
        default = 0.05,
        min = 0.01,
        max = 1.0
        )

    fTubeWidth: FloatProperty(
        name = "tube width",
        description = "tube width",
        default = 0.01,
        min = 0.001,
        max = 0.1
        )

# ------------------------------------------------------------------------
#    Operators
# ------------------------------------------------------------------------

class WM_OT_edges2girders(Operator):
    bl_label  = "edges2girders"
    bl_idname = "wm.edges2girders"

    def execute(self, context):
        scene = context.scene
        girdertool = scene.girder_tool
        
        edges2girders( dl=girdertool.fStep, w=girdertool.fWidth )

        return {'FINISHED'}

class WM_OT_edges2tubes(Operator):
    bl_label  = "edges2tubes"
    bl_idname = "wm.edges2tubes"

    def execute(self, context):
        scene = context.scene
        girdertool = scene.girder_tool
        
        edges2tubes( w=girdertool.fTubeWidth )

        return {'FINISHED'}

# ------------------------------------------------------------------------
#    Menus
# ------------------------------------------------------------------------

#class OBJECT_MT_CustomMenu(bpy.types.Menu):
class OBJECT_MT_GirderMenu(bpy.types.Menu):
    bl_label = "Select"
    bl_idname = "OBJECT_MT_girder_menu"

    def draw(self, context):
        layout = self.layout
        # Built-in operators
        layout.operator("object.select_all",    text="Select/Deselect All").action = 'TOGGLE'
        layout.operator("object.select_all",    text="Inverse").action = 'INVERT'
        layout.operator("object.select_random", text="Random")

# ------------------------------------------------------------------------
#    Panel in Object Mode
# ------------------------------------------------------------------------

#class OBJECT_PT_CustomPanel(Panel):
class OBJECT_PT_GirderPanel(Panel):
    bl_label = "Girderizer"
    bl_idname = "OBJECT_PT_girder_panel"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "Tools"
    bl_context = "objectmode"


    @classmethod
    def poll(self,context):
        return context.object is not None

    def draw(self, context):
        layout = self.layout
        scene  = context.scene
        girdertool = scene.girder_tool

        layout.prop(girdertool, "fStep")
        layout.prop(girdertool, "fWidth")
        layout.prop(girdertool, "fTubeWidth")

        layout.operator("wm.edges2girders")
        layout.operator("wm.edges2tubes")
        #layout.menu(OBJECT_MT_GirderMenu.bl_idname, text="Presets", icon="SCENE")
        layout.separator()

# ------------------------------------------------------------------------
#    Registration
# ------------------------------------------------------------------------

classes = (
    GirderProperties,
    WM_OT_edges2girders,
    WM_OT_edges2tubes,
    OBJECT_MT_GirderMenu,
    OBJECT_PT_GirderPanel
)

def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
    bpy.types.Scene.girder_tool = PointerProperty(type=GirderProperties)

def unregister():
    from bpy.utils import unregister_class
    for cls in reversed(classes):
        unregister_class(cls)
    del bpy.types.Scene.girder_tool


if __name__ == "__main__":
    register()

