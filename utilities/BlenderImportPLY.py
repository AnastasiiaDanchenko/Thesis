import bpy
import os

def add_particle_spheres(object_name="Fluid", sphere_radius=0.1):  
    obj = bpy.data.objects.get(object_name)
    if obj is None:
        print(f"Object '{object_name}' not found.")
        return

    sphere_collection = bpy.data.collections.new("Particle_Spheres")
    bpy.context.scene.collection.children.link(sphere_collection)

    bpy.ops.mesh.primitive_uv_sphere_add(radius=sphere_radius, segments=16, ring_count=8)
    base_sphere = bpy.context.object
    base_sphere.name = "Particle_Sphere"

    spheres = []
    for vertex in obj.data.vertices:
        sphere = base_sphere.copy()
        sphere.data = base_sphere.data.copy()
        sphere.name = f"Particle_{vertex.index}"

        sphere_collection.objects.link(sphere)
        sphere.location = vertex.co
        sphere.parent = obj
        spheres.append(sphere)

    bpy.data.objects.remove(base_sphere, do_unlink=True)


def import_ply(object_name="Fluid", nbFrames=100):
    directory = "C:\\dev\\MasterProject\\Thesis\\output\\blender_frames"

    firstFrame = os.path.join(directory, "0.ply")
    if not os.path.exists(firstFrame):
        print(f"File {firstFrame} does not exist")
        return

    bpy.ops.import_mesh.ply(filepath=firstFrame)
    obj = bpy.context.selected_objects[0]
    obj.name = object_name
    bpy.context.view_layer.objects.active = obj

    bpy.ops.object.shape_key_add(from_mix=False)
    obj.data.shape_keys.key_blocks[-1].name = "Base"

    for frame in range(1, nbFrames):
        filename = os.path.join(directory, f"{frame}.ply")

        if not os.path.exists(filename):
            print(f"File {filename} does not exist")
            continue

        with open(filename, 'r') as file:
            lines = file.readlines()

        vertex_data_start = lines.index("end_header\n") + 1
        vertex_positions = [list(map(float, line.strip().split())) for line in lines[vertex_data_start:]]

        if len(vertex_positions) != len(obj.data.vertices):
            raise ValueError(f"Frame {frame} has a different number of vertices than the initial frame!")
        
        bpy.ops.object.shape_key_add(from_mix=True)
        key_block = obj.data.shape_keys.key_blocks[-1]
        key_block.name = f"Frame_{frame}"

        for i, vertex in enumerate(obj.data.vertices):
            key_block.data[i].co = vertex_positions[i]

        key_block.value = 0.0
        key_block.keyframe_insert(data_path="value", frame=frame-1)
        key_block.value = 1.0
        key_block.keyframe_insert(data_path="value", frame=frame)
        key_block.value = 0.0
        key_block.keyframe_insert(data_path="value", frame=frame+1)

if __name__ == "__main__":
    object_name = "Dam Break"
    import_ply(object_name, 100)
    add_particle_spheres(object_name=object_name, sphere_radius=0.1)
