#include "moab/Interface.hpp"
#include "moab/Core.hpp"
#include "moab/Types.hpp"
#include "MBTagConventions.hpp"
#include "moab/GeomTopoTool.hpp"
#include <iostream>
#include <map>
#include <set>

using namespace moab;

Tag category_tag;
Tag geom_tag;
Tag name_tag;
Tag obj_name_tag;
Tag dim_tag, id_tag;

ErrorCode build_cube( Interface *mbi,
                      std::vector<double> scale_vec,
                      std::vector<double> trans_vec,
                      int    object_id,
                      EntityHandle &volume  )
{
  GeomTopoTool *GTT = new GeomTopoTool(mbi);

  ErrorCode rval;

  // Define a 1x1x1 cube centered at orgin

  // coordinates of each corner
  const double coords[] = {
    0.5, -0.5, -0.5,
    0.5,  0.5, -0.5,
   -0.5,  0.5, -0.5,
   -0.5, -0.5, -0.5,
    0.5, -0.5,  0.5,
    0.5,  0.5,  0.5,
   -0.5,  0.5,  0.5,
   -0.5, -0.5,  0.5 };

  // connectivity of 2 triangles per
  //  each face of the cube
  const int connectivity[] = {
    0, 3, 1,  3, 2, 1, // -Z
    0, 1, 4,  5, 4, 1, // +X
    1, 2, 6,  6, 5, 1, // +Y
    6, 2, 3,  7, 6, 3, // -X
    0, 4, 3,  7, 3, 4, // -Y
    4, 5, 7,  5, 6, 7 // +Z
  };


  // Create the geometry
  const int num_verts = 8;
  const int num_tris = 12;
  EntityHandle verts[num_verts], tris[num_tris], surf;

  rval = mbi->create_meshset( MESHSET_SET, surf ); MB_CHK_ERR(rval);
/*
  // scale coords
  int i;
  double scaled_coords[24];
  for ( i = 0; i < num_verts; i++ )
    {
      scaled_coords[3*i]   = coords[3*i]*scale_vec[0];
      scaled_coords[3*i+1] = coords[3*i+1]*scale_vec[1];
      scaled_coords[3*i+2] = coords[3*i+2]*scale_vec[2];
    }

  // translate coords
  double trans_coords[24];
  for ( i = 0; i < num_verts; i++ )
    {
      trans_coords[3*i]   = scaled_coords[3*i] + trans_vec[0];
      trans_coords[3*i+1] = scaled_coords[3*i+1] + trans_vec[1];
      trans_coords[3*i+2] = scaled_coords[3*i+2] + trans_vec[2];
    }
*/
  // transform coords-- scale and translate
  // double trans_coords[24];
  // for ( int i = 0; i < num_verts; i++ )
  //   {
  //     trans_coords[3*i]   = coords[3*i]*scale_vec[0]   + trans_vec[0];
  //     trans_coords[3*i+1] = coords[3*i+1]*scale_vec[1] + trans_vec[1];
  //     trans_coords[3*i+2] = coords[3*i+2]*scale_vec[2] + trans_vec[2];
  //   }
  //
  // // create vertices and add to meshset
  // for ( int i = 0; i < num_verts; ++i)
  //   {
  //     rval = mbi->create_vertex( trans_coords + 3*i, verts[i] ); MB_CHK_ERR(rval);
  //
  //     rval = mbi->add_entities( surf, &verts[i], 1 ); MB_CHK_ERR(rval);
  //
  //   }
  //
  // // create triangles and add to meshset
  // for ( int i = 0; i < num_tris; ++i)
  //   {
  //     const EntityHandle conn[] = { verts[connectivity[3*i  ]],
  //                                   verts[connectivity[3*i+1]],
  //                                   verts[connectivity[3*i+2]] };
  //     rval = mbi->create_element( MBTRI, conn, 3, tris[i] ); MB_CHK_ERR(rval);
  //
  //     rval = mbi->add_entities( surf, &tris[i], 1 ); MB_CHK_ERR(rval);
  //   }
  //
  //
  // // set name, id, geom, and category tags for SURFACE
  // rval = mbi->tag_set_data( name_tag, &surf, 1, "Surface\0" ); MB_CHK_ERR(rval);
  // std::string object_name;
  // rval = mbi->tag_set_data( obj_name_tag, &surf, 1, object_name.c_str() ); MB_CHK_ERR(rval);
  // rval = mbi->tag_set_data( id_tag, &surf, 1, &object_id ); MB_CHK_ERR(rval);
  // int two = 2;
  // rval = mbi->tag_set_data( geom_tag, &surf, 1, &(two) ); MB_CHK_ERR(rval);
  // rval = mbi->tag_set_data( category_tag, &surf, 1, "Surface\0" ); MB_CHK_ERR(rval);
  //
  // // create volume meshset associated with surface meshset
  // // EntityHandle volume;
  // rval = mbi->create_meshset( MESHSET_SET, volume ); MB_CHK_ERR(rval);
  //
  // // set name, id, geom, and category tags for VOLUME
  // rval = mbi->tag_set_data( name_tag, &volume, 1, "Volume\0" ); MB_CHK_ERR(rval);
  // rval = mbi->tag_set_data( obj_name_tag, &surf, 1, object_name.c_str() ); MB_CHK_ERR(rval);
  // rval = mbi->tag_set_data( id_tag, &volume, 1, &(object_id) ); MB_CHK_ERR(rval);
  // int three = 3;
  // rval = mbi->tag_set_data( geom_tag, &volume, 1, &(three) ); MB_CHK_ERR(rval);
  // rval = mbi->tag_set_data( category_tag, &volume, 1, "Volume\0" ); MB_CHK_ERR(rval);
  //
  //
  // // set surface as child of volume
  // rval = mbi->add_parent_child( volume, surf ); MB_CHK_ERR(rval);
  //
  // // set sense tag
  // rval = GTT->set_sense(surf, volume, SENSE_FORWARD); MB_CHK_ERR(rval);

  delete GTT;

  return MB_SUCCESS;
}

int main()
{
    return 0;
}
