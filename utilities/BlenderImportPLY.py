import bpy
import os

import bpy

def create_and_animate_spheres(object_name, sphere_radius=0.1):
    # Get the mesh object containing the vertices
    obj = bpy.data.objects.get(object_name)
    if obj is None:
        print(f"Object '{object_name}' not found.")
        return
    
    # Ensure the object is a mesh
    if obj.type != 'MESH':
        print(f"Object '{object_name}' is not a mesh.")
        return
    
    # Create a new collection for the spheres
    sphere_collection = bpy.data.collections.new("Particle_Spheres")
    bpy.context.scene.collection.children.link(sphere_collection)

    # Store spheres and vertex index mapping
    spheres = []
    vertex_to_sphere = {}

    # Create spheres and position them at the vertices
    for vertex in obj.data.vertices:
        # Create a new sphere
        bpy.ops.mesh.primitive_uv_sphere_add(radius=sphere_radius, segments=16, ring_count=8)
        sphere = bpy.context.object
        sphere.name = f"Particle_Sphere_{vertex.index}"
        
        # Link the sphere to the collection
        sphere_collection.objects.link(sphere)
        bpy.context.collection.objects.unlink(sphere)
        
        # Position the sphere at the vertex location
        sphere.location = obj.matrix_world @ vertex.co
        
        # Store the sphere and vertex index mapping
        spheres.append(sphere)
        vertex_to_sphere[vertex.index] = sphere

    # Add vertex groups to the mesh object
    mesh = obj.data
    for vertex_index in vertex_to_sphere.keys():
        group_name = f"VertexGroup_{vertex_index}"
        if group_name not in mesh.vertex_groups:
            mesh.vertex_groups.new(name=group_name)
        group = mesh.vertex_groups[group_name]
        group.add([vertex_index], 1.0, 'REPLACE')

    # Create shape keys for animation
    if not mesh.shape_keys:
        bpy.ops.object.shape_key_add(from_mix=False)
    base_key = mesh.shape_keys.key_blocks[0]
    
    # Ensure there is a shape key for each frame
    for i in range(1, 101):  # Example: 100 frames
        if i >= len(mesh.shape_keys.key_blocks):
            bpy.ops.object.shape_key_add(from_mix=False)
        shape_key = mesh.shape_keys.key_blocks[i]
        shape_key.name = f"Frame_{i}"

        # Set the influence of the shape key to the corresponding sphere
        for vertex_index, sphere in vertex_to_sphere.items():
            group = mesh.vertex_groups.get(f"VertexGroup_{vertex_index}")
            if group:
                group.add([vertex_index], shape_key.value, 'REPLACE')

    # Set the initial frame to ensure proper animation
    bpy.context.scene.frame_set(1)

    print(f"Created and animated {len(spheres)} spheres for particles.")

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

def import_ply_sequence(nbFrames=100):
    directory = "C:\\dev\\MasterProject\\Thesis\\output\\blender_frames"

    imported_objects = []

    for frame in range(nbFrames):
        filename = os.path.join(directory, f"{frame}_r.ply")
        
        if not os.path.exists(filename):
            print(f"File {filename} does not exist")
            continue
        
        # Import the PLY file
        bpy.ops.import_mesh.ply(filepath=filename)
        obj = bpy.context.selected_objects[0]
        obj.name = f"{object_name}_{frame}"
        
        # Hide the object by default
        obj.hide_viewport = True
        obj.hide_render = True
        obj.keyframe_insert(data_path="hide_viewport", frame=frame-1 if frame > 0 else 0)
        obj.keyframe_insert(data_path="hide_render", frame=frame-1 if frame > 0 else 0)
        
        # Unhide the object on its corresponding frame
        obj.hide_viewport = False
        obj.hide_render = False
        obj.keyframe_insert(data_path="hide_viewport", frame=frame)
        obj.keyframe_insert(data_path="hide_render", frame=frame)
        
        # Hide it again immediately after (next frame)
        obj.hide_viewport = True
        obj.hide_render = True
        obj.keyframe_insert(data_path="hide_viewport", frame=frame+1)
        obj.keyframe_insert(data_path="hide_render", frame=frame+1)

        imported_objects.append(obj)

    # Optionally, you can parent all these objects to an empty or another object for easier manipulation
    if imported_objects:
        bpy.ops.object.empty_add(type='PLAIN_AXES', location=(0, 0, 0))
        parent_obj = bpy.context.object
        parent_obj.name = f"{object_name}_Controller"
        
        for obj in imported_objects:
            obj.parent = parent_obj


if __name__ == "__main__":
    object_name = "Dam Break"
    frames_num = 100
    import_ply_sequence(nbFrames=frames_num)
    # import_ply(object_name, frames_num)
    # create_and_animate_spheres(object_name="Fluid", sphere_radius=0.1)

