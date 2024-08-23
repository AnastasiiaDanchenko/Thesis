import bpy
import os

def import_ply(nbFrames=100):
    directory = "C:\\dev\\MasterProject\\Thesis\\output\\blender_frames"

    firstFrame = os.path.join(directory, "0.ply")
    if not os.path.exists(firstFrame):
        print(f"File {firstFrame} does not exist")
        return

    bpy.ops.import_mesh.ply(filepath=firstFrame)
    obj = bpy.context.selected_objects[0]
    obj.name = "Fluid"
    # obj.keyframe_insert(data_path="location", frame=1)

    # Ensure the object is active
    bpy.context.view_layer.objects.active = obj

    # Create a shape key for the base frame (rest position)
    bpy.ops.object.shape_key_add(from_mix=False)
    obj.data.shape_keys.key_blocks[-1].name = "Base"

    for frame in range(1, nbFrames):
        filename = os.path.join(directory, f"{frame}.ply")

        if not os.path.exists(filename):
            print(f"File {filename} does not exist")
            continue

        # bpy.ops.import_mesh.ply(filepath=filename)
        
        # # Set the imported object to the current frame
        # obj = bpy.context.selected_objects[0]
        # obj.name = f"frame_{frame}"
        # obj.animation_data_create()
        # obj.animation_data.action = bpy.data.actions.new(name="Action")
        # obj.location = (0, 0, 0)
        
        # # Insert keyframes to show/hide the object at the correct frame
        # obj.hide_viewport = True
        # obj.keyframe_insert(data_path="hide_viewport", frame=frame)
        # obj.hide_viewport = False
        # obj.keyframe_insert(data_path="hide_viewport", frame=frame + 1)

        with open(filename, 'r') as file:
            lines = file.readlines()

        vertex_data_start = lines.index("end_header\n") + 1
        vertex_positions = [list(map(float, line.strip().split())) for line in lines[vertex_data_start:]]

        # Ensure the mesh has the correct number of vertices
        if len(vertex_positions) != len(obj.data.vertices):
            raise ValueError(f"Frame {frame} has a different number of vertices than the initial frame!")

        # # Update the mesh vertices
        # for i, vertex in enumerate(obj.data.vertices):
        #     vertex.co.x = vertex_positions[i][0]
        #     vertex.co.y = vertex_positions[i][1]
        #     vertex.co.z = vertex_positions[i][2]

        # # Insert keyframes to animate vertex positions
        # obj.data.update()
        # obj.keyframe_insert(data_path="data.vertices[:].co", frame=frame)

        # Create a new shape key for this frame
        bpy.ops.object.shape_key_add(from_mix=False)
        key_block = obj.data.shape_keys.key_blocks[-1]
        key_block.name = f"Frame_{frame}"

        # Update the mesh vertices for the shape key
        for i, vertex in enumerate(obj.data.vertices):
            key_block.data[i].co = vertex_positions[i]

        # Insert keyframe for shape key value
        key_block.value = 0.0
        key_block.keyframe_insert(data_path="value", frame=frame-1)
        key_block.value = 1.0
        key_block.keyframe_insert(data_path="value", frame=frame)

if __name__ == "__main__":
    import_ply(100)