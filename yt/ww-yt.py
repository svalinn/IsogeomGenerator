import yt

mesh = yt.load("../magic/mesh-ww-tags.h5m")

min_ww, min_loc = mesh.find_min("ww_n")
max_ww, max_loc = mesh.find_max("ww_n")

print(mesh.find_field_values_at_point("ww_n", [0.,0.,0.]))

dd = mesh.all_data()
dd.extract_isocontours("ww_n", 0.002, "test.obj", True)
#print(dd.quantities)
#print(mesh.field_list)

#mesh.surface(mesh, 'ww_n', 0.1)
