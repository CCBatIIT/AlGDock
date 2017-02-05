Timer unit: 1e-06 s

Total time: 125.031 s
File: /Users/dminh/.ipython/cython/_cython_magic_66ad0cdc396f48457a5191d065a19ef4.pyx
Function: calc_desolvationGrid at line 389

Line #      Hits         Time  Per Hit   % Time  Line Contents
==============================================================
   389                                           def calc_desolvationGrid(int_t[:,:,:] receptor_MS_grid, float_t[:] spacing, int_t[:] counts, \
   390                                                                      float_t[:,:] receptor_SAS_points, float_t[:,:] receptor_coordinates, \
   391                                                                      float_t[:,:] SAS_sphere_pts, float_t[:] LJ_r2, \
   392                                                                      float_t ligand_atom_radius, float_t probe_radius, float_t integration_cutoff):
   393                                             cdef int_t i, j, k, n
   394                                             cdef int_t nreceptor_SAS_points
   395                                             cdef float_t x_min, x_max, y_min, y_max, z_min, z_max
   396                                             cdef float_t ligand_atom_radius2
   397                                             cdef float_t grid_point_x, grid_point_y, grid_point_z
   398                                             cdef np.ndarray[np.double_t, ndim=3] desolvationGrid
   399                                           
   400                                             # To use for new SAS points
   401                                             cdef int_t clash
   402                                             cdef int_t nreceptor_atoms, nsphere_points, n_newly_inaccessible_SAS_points
   403                                             cdef int_t atom_j, sphere_i, d
   404                                             cdef float_t SAS_point_x, SAS_point_y, SAS_point_z
   405                                             cdef float_t dx, dy, dz
   406                                             cdef np.ndarray[int_t, ndim=3] grid_c
   407                                               
   408         1            2      2.0      0.0    ligand_atom_radius2 = ligand_atom_radius*ligand_atom_radius
   409         1            1      1.0      0.0    nreceptor_SAS_points = receptor_SAS_points.shape[0]
   410         1           16     16.0      0.0    nreceptor_atoms = len(LJ_r2)
   411         1            6      6.0      0.0    nsphere_points = len(SAS_sphere_pts)  
   412                                             
   413         1          296    296.0      0.0    desolvationGrid = np.zeros(shape=tuple(counts), dtype=np.float)
   414                                             
   415                                           #   i = 0
   416                                           #   while i < counts[0]:
   417                                           #     grid_point_x = i*spacing[0]
   418                                           #     j = 0
   419                                           #     while j < counts[1]:
   420                                           #       grid_point_y = j*spacing[1]
   421                                           #       k = 45 # TODO: Loop over k
   422                                           #       grid_point_z = k*spacing[2]
   423         1            1      1.0      0.0    i = 26
   424         1            0      0.0      0.0    while i < 27:
   425         1            1      1.0      0.0      grid_point_x = i*spacing[0]
   426         1            0      0.0      0.0      j = 0
   427         1            1      1.0      0.0      while j < counts[1]:
   428        73           49      0.7      0.0        grid_point_y = j*spacing[1]
   429        73           42      0.6      0.0        k = 45 # TODO: Loop over k
   430        73           40      0.5      0.0        grid_point_z = k*spacing[2]
   431                                           
   432        73          522      7.2      0.0        newly_inaccessible_SAS_points = []
   433                                           
   434        73           40      0.5      0.0        x_min = grid_point_x - ligand_atom_radius
   435        73           52      0.7      0.0        x_max = grid_point_x + ligand_atom_radius
   436        73           35      0.5      0.0        y_min = grid_point_y - ligand_atom_radius
   437        73           59      0.8      0.0        y_max = grid_point_y + ligand_atom_radius
   438        73           39      0.5      0.0        z_min = grid_point_z - ligand_atom_radius
   439        73           48      0.7      0.0        z_max = grid_point_z + ligand_atom_radius
   440                                           
   441        73           37      0.5      0.0        n = 0
   442        73           45      0.6      0.0        while n < nreceptor_SAS_points:
   443  12176254      7058454      0.6      5.6          SAS_point_x = receptor_SAS_points[n,0]
   444  12176254      7061936      0.6      5.6          SAS_point_y = receptor_SAS_points[n,1]
   445  12176254      7066737      0.6      5.7          SAS_point_z = receptor_SAS_points[n,2]
   446  24352508     14093985      0.6     11.3          if (SAS_point_x>x_min) and (SAS_point_y>y_min) and \
   447   3809636      2207269      0.6      1.8             (SAS_point_z>z_min) and \
   448   1814184      1054232      0.6      0.8             (SAS_point_x<x_max) and (SAS_point_y<y_max) and \
   449     43606        25666      0.6      0.0             (SAS_point_z<z_max):
   450      7055         4031      0.6      0.0            dx = SAS_point_x-grid_point_x
   451      7055         4080      0.6      0.0            dy = SAS_point_y-grid_point_y
   452      7055         4052      0.6      0.0            dz = SAS_point_z-grid_point_z
   453      7055         4132      0.6      0.0            if (dx*dx + dy*dy + dz*dz)<ligand_atom_radius2:
   454      3577         2505      0.7      0.0              newly_inaccessible_SAS_points.append(\
   455      3577         2823      0.8      0.0                (SAS_point_x,SAS_point_y,SAS_point_z))
   456  12176254      7050012      0.6      5.6          n += 1
   457                                                 
   458        73           46      0.6      0.0        n_newly_inaccessible_SAS_points = len(newly_inaccessible_SAS_points)
   459        73           49      0.7      0.0        if n_newly_inaccessible_SAS_points==0:
   460        12       295235  24602.9      0.2          desolvationGrid[i,j,k] = fraction_r4inv_low_dielectric(\
   461                                                       receptor_MS_grid, spacing, counts, \
   462                                                       grid_point_x, grid_point_y, grid_point_z, \
   463                                                       ligand_atom_radius, integration_cutoff)
   464                                                 else:
   465        61        34230    561.1      0.0          grid_c = np.copy(receptor_MS_grid)
   466                                           
   467                                                   # Find new SAS points around the ligand atom and 
   468                                                   # increment the marks of the grid points within a probe radius
   469        61          125      2.0      0.0          sphere_i = 0
   470        61           56      0.9      0.0          while sphere_i < nsphere_points:
   471                                                     # Propose a point at the SAS of the atom
   472     15250         9430      0.6      0.0            SAS_point_x = SAS_sphere_pts[sphere_i,0] + grid_point_x
   473     15250         9070      0.6      0.0            SAS_point_y = SAS_sphere_pts[sphere_i,1] + grid_point_y
   474     15250         9478      0.6      0.0            SAS_point_z = SAS_sphere_pts[sphere_i,2] + grid_point_z
   475     15250         8879      0.6      0.0            clash = 0
   476     15250         8901      0.6      0.0            atom_j = 0
   477     15250         8713      0.6      0.0            while atom_j < nreceptor_atoms:
   478  26785867     15446422      0.6     12.4              dx = receptor_coordinates[atom_j,0] - SAS_point_x
   479  26785867     15398550      0.6     12.3              dy = receptor_coordinates[atom_j,1] - SAS_point_y
   480  26785867     15390550      0.6     12.3              dz = receptor_coordinates[atom_j,2] - SAS_point_z
   481  26785867     15442672      0.6     12.4              if (dx*dx + dy*dy + dz*dz) < LJ_r2[atom_j]:
   482      7358         4263      0.6      0.0                clash = 1
   483      7358         4650      0.6      0.0                atom_j = nreceptor_atoms
   484                                                       else:
   485  26778509     15371695      0.6     12.3                atom_j += 1
   486     15250         8914      0.6      0.0            if clash==0:
   487                                                       # Increment the marks of the grid points within a probe radius
   488      7892        57178      7.2      0.0              increment_inside_sphere(grid_c, spacing, counts, \
   489      7892       544770     69.0      0.4                SAS_point_x, SAS_point_y, SAS_point_z, probe_radius)
   490     15250         9557      0.6      0.0            sphere_i += 1
   491                                           
   492                                                   # Decrement of the marks of grid points within the probe radius
   493        61           55      0.9      0.0          while n < n_newly_inaccessible_SAS_points:
   494                                                     decrement_inside_sphere(grid_c, spacing, counts, \
   495                                                       newly_inaccessible_SAS_points[n][0], \
   496                                                       newly_inaccessible_SAS_points[n][1], \
   497                                                       newly_inaccessible_SAS_points[n][2], probe_radius)
   498                                                   # Blot the region inside the ligand vdW as low dielectric
   499        61          340      5.6      0.0          set_inside_sphere_to(grid_c, spacing, counts, \
   500                                                     grid_point_x, grid_point_y, grid_point_z, \
   501        61         4076     66.8      0.0            ligand_atom_radius, 0)
   502                                           
   503       122          389      3.2      0.0          desolvationGrid[i,j,k] = fraction_r4inv_low_dielectric(grid_c, \
   504                                                     spacing, counts, grid_point_x, grid_point_y, grid_point_z, \
   505        61      1321086  21657.1      1.1            ligand_atom_radius, integration_cutoff)
   506                                           
   507        73           73      1.0      0.0        j += 1
   508         1            1      1.0      0.0      i += 1
   509                                           
   510         1          404    404.0      0.0    return desolvationGrid
