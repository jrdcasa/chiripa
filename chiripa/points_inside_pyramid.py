import numpy as np
import random
import os
import chiripa as chi

# #########################################################################
def calculateNormal(p1, p2, p3):

    """
    Calculate the normal vector of a plane from three points in it.

    From three vertices, a triangle plane is constructed. The vertice p1 must be the begin point
    of the two edges vectors.

    .. image:: ../imgs/normal_vector.png

    Args:
        p1 (ndarray): 1D array containing the Point 1. Common point of the two edges vectors.
        p2 (ndarray): 1D array containing the Point 2
        p3 (ndarray): 1D array containing the Point 3

    Returns:
        A 1D array containing the normal vector at the origin (0,0,0)

    """

    v1 = p1 - p3
    v2 = p1 - p2
    n = np.cross(v1, v2)
    n = n/np.linalg.norm(n)

    return n

# #########################################################################
def calculate_all_normals_cube(Lx, vertices=None):
    """
    This only works with cubic boxes, with the origin in [0, 0, 0].
    The coordinate reference is placed as:

    Reference system:

    .. image:: ../imgs/referencesystem.png


    Args:
        Lx (float): Length of the cube
        vertices (list): [alpha, beta, gamma]


    Returns:
        A tuple containing the **normals** and the **dots** (D) coeficients of the planes.

        * normals (ndarray): Shape [6,5,3]. Each row contains all normal vectors (5) of each \
        of the pyramids to divide in 6 regions the cube.
        * dots (ndarray): Shape [6,5]. Each file contains the value D of each plane equation (nx*x+ny*y+nz*z+D=0)

    """

    if vertices is None:
        # Vertices p1 to p8 for a cube with origin (0,0,0) and perfecly aligned with the reference axis
        c = Lx / 2.0
        # Vertices and origin
        p1 = np.array([ c, -c, c])
        p2 = np.array([-c, -c, c])
        p3 = np.array([-c,  c, c])
        p4 = np.array([ c,  c, c])

        p5 = np.array([ c, -c, -c])
        p6 = np.array([-c, -c, -c])
        p7 = np.array([-c,  c, -c])
        p8 = np.array([ c,  c, -c])
        o = np.array([0.0, 0.0, 0.0])
    else:
        p1 = vertices[0]
        p2 = vertices[1]
        p3 = vertices[2]
        p4 = vertices[3]

        p5 = vertices[4]
        p6 = vertices[5]
        p7 = vertices[6]
        p8 = vertices[7]
        o = np.array([0.0, 0.0, 0.0])

    # Top pyramid Z ==========================================
    # normals
    n1 = list()
    a = calculateNormal(p1, p2, p3)
    b = calculateNormal( o, p2, p1)  # Normal vector outside
    c = calculateNormal( o, p3, p2)  # Normal vector outside
    d = calculateNormal( o, p4, p3)  # Normal vector outside
    e = calculateNormal( o, p1, p4)  # Normal vector outside
    n1 = [a, b, c, d, e]

    # dots
    # Plane equation [nx*x+ny*y+nz*z+D = 0] where D = - (nx*p1x + ny*p1y + nz*p1z)  --> -np.dot(n,p1)
    d1 = list()
    a = - np.dot(n1[0], p1)
    b = - np.dot(n1[1], o)
    c = - np.dot(n1[2], o)
    d = - np.dot(n1[3], o)
    e = - np.dot(n1[4], o)
    d1 = [a, b, c, d, e]


    # Bottom pyramid Z=======================================
    n2 = list()
    a = calculateNormal(p6, p5, p7)
    b = calculateNormal( o, p5, p6)
    c = calculateNormal( o, p6, p7)
    d = calculateNormal( o, p7, p8)
    e = calculateNormal( o, p8, p5)
    n2 = [a, b, c, d, e]

    # dots
    # Plane equation [nx*x+ny*y+nz*z+D = 0] where D = - (nx*p1x + ny*p1y + nz*p1z)  --> -np.dot(n,p1)
    d2 = list()
    a = - np.dot(n2[0], p6)
    b = - np.dot(n2[1], o)
    c = - np.dot(n2[2], o)
    d = - np.dot(n2[3], o)
    e = - np.dot(n2[4], o)
    d2 = [a, b, c, d, e]

    # Front pyramid X========================================
    n3 = list()
    a = calculateNormal(p4, p5, p1)
    b = calculateNormal(o, p4, p1)
    c = calculateNormal(o, p1, p5)
    d = calculateNormal(o, p5, p8)
    e = calculateNormal(o, p8, p4)
    n3 = [a, b, c, d, e]

    # dots
    # Plane equation [nx*x+ny*y+nz*z+D = 0] where D = - (nx*p1x + ny*p1y + nz*p1z) --> -np.dot(n,p1)
    d3 = list()
    a = - np.dot(n3[0], p4)
    b = - np.dot(n3[1], o)
    c = - np.dot(n3[2], o)
    d = - np.dot(n3[3], o)
    e = - np.dot(n3[4], o)
    d3 = [a, b, c, d, e]

    # Bottom pyramid X========================================
    n4 = list()
    a = calculateNormal(p2, p6, p7)
    b = calculateNormal( o, p6, p2)
    c = calculateNormal( o, p7, p6)
    d = calculateNormal( o, p3, p7)
    e = calculateNormal( o, p2, p3)
    n4 = [a, b, c, d, e]

    # dots
    # Plane equation [nx*x+ny*y+nz*z+D = 0] where D = - (nx*p1x + ny*p1y + nz*p1z) --> -np.dot(n,p1)
    d4 = list()
    a = - np.dot(n4[0], p2)
    b = - np.dot(n4[1], o)
    c = - np.dot(n4[2], o)
    d = - np.dot(n4[3], o)
    e = - np.dot(n4[4], o)
    d4 = [a, b, c, d, e]

    # Front pyramid Y========================================
    n5 = list()
    a = calculateNormal(p3, p8, p4)
    b = calculateNormal( o, p4, p8)
    c = calculateNormal( o, p8, p7)
    d = calculateNormal( o, p7, p3)
    e = calculateNormal( o, p3, p4)
    n5 = [a, b, c, d, e]

    # dots
    # Plane equation [nx*x+ny*y+nz*z+D = 0] where D = - (nx*p1x + ny*p1y + nz*p1z) --> -np.dot(n,p1)
    d5 = list()
    a = - np.dot(n5[0],p3)
    b = - np.dot(n5[1], o)
    c = - np.dot(n5[2], o)
    d = - np.dot(n5[3], o)
    e = - np.dot(n5[4], o)
    d5 = [a, b, c, d, e]

    # Bottom pyramid Y========================================
    n6 = list()
    a = calculateNormal(p2, p1, p5)
    b = calculateNormal( o, p1, p2)
    c = calculateNormal( o, p5, p1)
    d = calculateNormal( o, p6, p5)
    e = calculateNormal( o, p2, p6)
    n6 = [a, b, c, d, e]

    # dots
    # Plane equation [nx*x+ny*y+nz*z+D = 0] where D = - (nx*p1x + ny*p1y + nz*p1z) --> -np.dot(n,p1)
    d6 = list()
    a = - np.dot(n6[0], p2)
    b = - np.dot(n6[1], o)
    c = - np.dot(n6[2], o)
    d = - np.dot(n6[3], o)
    e = - np.dot(n6[4], o)
    d6 = [a, b, c, d, e]

    # Stack normals in rows. Each row is a pyramid (axis = 0), each column (axis=1)
    # is one of the normals in that pyramid and each axis=2 is a component of the
    # normal vector
    normals = np.stack((n1, n2, n3, n4, n5, n6))
    dots = np.stack((d1, d2, d3, d4, d5, d6))
    vertices = np.stack((p1, p2, p3, p4, p5, p6, p7, p8))

    return normals, dots, vertices

# #########################################################################
def testPointInsidePyramid(boxl, list_of_points, nrotations=100,
                           debug=False, iseed=None):

    """
    This function divides a cubic region in 6 subregions defined by the six squared-pyramids centered in the cube
    origin. The origin of the cube must be [0.0, 0.0, 0.0].
    Then, it calculates the index of anisotropy (S) ``nrotations`` times with different divisions of the space.
    Subsequenly, both the averaged and the standard deviation are returned.
    The way to calculate the index of anisotropy (S) is closed to the method proposed by
    `Okuwaki et al. <https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.7b08461>`_

    .. image:: ../imgs/pyramids.png

    The algorithm is:

    1. Repeat the calculations for ``nrotations`` of the cube.

        a. In the first iteration (irot=0), the cube is aligned with the axis.
        Get the p1-p8 vertices (8 vertices) of the cube.

        b. The rest of iterations (irot != 0) the cube is generated by random rotation of the vertices
        using a tuple of euler angles
        (:py:meth:`chiripa.internal_coordinates.generate_random_euler_angles` and
        :py:meth:`chiripa.internal_coordinates.euler_rotation_matrix` functions).
        Get the p1-p8 vertices (8 vertices) of the cube.

        .. image:: ../imgs/cube_rotations.png

        c. Calculate all outward normal vectors and D coeffients of the five (4 triangles and 1 square) faces
        in each pyramid. Use the :py:meth:`calculate_all_normals_cube` function.

        d. For a point p(px, py, pz) and a plane(nx, ny, nz, D), the sign of the Dist value
        is able to determinate which half space the point is in.

        .. image:: ../imgs/distanceplane.png

        e. If all signs of the np.dot(p,n_i)+D are negative then the point is inside the region.



    Args:
        boxl (float): Cubic length in angstroms.
        list_of_points (list of ndarray(3)): All points to be evaluated in the regions.
            Usually, these points are the center of geometry of the generated confirmations.
        nrotations (int): Number of rotations to be performed to calculate the value
            of asymmetry parameter. Default=100
        debug (bool): If True write debug information in a directory called ROTATIONS. Default=False
        iseed (int): Seed for the random number generator. If is ``None`` then the iseed is
            randomly generated. Default=None

    Returns:
        tuple: The mean and the standard deviation for the index of anisotropy (S)

    .. warning:: This function only works with cubic boxes of length `boxl` centered in [0., 0., 0.].

    """

    # Seed for the random generator ========
    if iseed is None:
        iseed = int(982628732*random.random())

    if os.path.isfile("vertex.xyz"):
        os.remove("vertex.xyz")
    if debug:
        if os.path.isdir("ROTATIONS"):
            os.system("rm -r ./ROTATIONS")
        os.mkdir("ROTATIONS")

    N_regions = 6
    S_list = []

    vertex = None
    for irot in range(0, nrotations):

        S = 0.0
        d_factor = np.zeros(6, dtype=int)

        if irot == 0:
            normals, dots, vertex = calculate_all_normals_cube(boxl)
        else:
            iseed += 1
            euler = chi.generate_random_euler_angles(seed=iseed)
            S = chi.euler_rotation_matrix(euler)
            a = np.transpose(np.array(vertex))
            v = np.transpose(np.matmul(S, a))
            normals, dots, vertex = calculate_all_normals_cube(boxl, vertices=v)

        if debug:
            irow,icol = vertex.shape
            with open("vertex.xyz", 'a') as f:
                f.writelines("{}\n".format(irow))
                f.writelines("Frame\n")
                for j in range(0, irow):
                    l = "P  {} {} {}\n".format(vertex[j][0], vertex[j][1], vertex[j][2])
                    f.writelines(l)

            debug_write_tcl_pyramid(irot, normals, vertex)
            i, j, k = normals.shape

        for ipoint in list_of_points:
            for iregion in range(normals.shape[0]):
                index = 0
                while True:
                    distPointPlane = np.dot(ipoint, normals[iregion, index, :]) + dots[iregion, index]
                    if distPointPlane > 0.0: break
                    if index == 4:
                        d_factor[iregion] += 1
                        break
                    index += 1

        #print(d_factor, S)
        d_factor = d_factor/d_factor.max()
        d_factor_f = ['%.2f' % item for item in d_factor]
        S = sum(d_factor/N_regions)
        S_list.append(S)
        #print(d_factor, S, S_list)

    S_avg = np.mean(S_list)
    S_std = np.std(S_list)

    return S_avg, S_std

# #########################################################################
def debug_write_tcl_pyramid(irot, normals, vertex):

    npyr, nnorm, dim = normals.shape

    # ===================================================
    # Write pyramid points
    # PYRAMID 1 -- Z above
    with open("ROTATIONS/pyramid_points_Zabove_{0:03d}.xyz".format(irot), 'w') as f:
        f.writelines("5\n")
        f.writelines("Pyramid Z above\n")
        f.writelines("O 0.0 0.0 0.0\n")
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[0][0], vertex[0][1], vertex[0][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[1][0], vertex[1][1], vertex[1][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[2][0], vertex[2][1], vertex[2][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[3][0], vertex[3][1], vertex[3][2]))

    # PYRAMID 2 -- Z below
    with open("ROTATIONS/pyramid_points_Zbelow_{0:03d}.xyz".format(irot), 'w') as f:
        f.writelines("5\n")
        f.writelines("Pyramid Z below\n")
        f.writelines("O 0.0 0.0 0.0\n")
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[4][0], vertex[4][1], vertex[4][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[5][0], vertex[5][1], vertex[5][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[6][0], vertex[6][1], vertex[6][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[7][0], vertex[7][1], vertex[7][2]))

    # PYRAMID 3 -- Y Above
    with open("ROTATIONS/pyramid_points_Yabove_{0:03d}.xyz".format(irot), 'w') as f:
        f.writelines("5\n")
        f.writelines("Pyramid Y above\n")
        f.writelines("O 0.0 0.0 0.0\n")
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[2][0], vertex[2][1], vertex[2][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[3][0], vertex[3][1], vertex[3][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[6][0], vertex[6][1], vertex[6][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[7][0], vertex[7][1], vertex[7][2]))

    # PYRAMID 4 -- Y Below
    with open("ROTATIONS/pyramid_points_Ybelow_{0:03d}.xyz".format(irot), 'w') as f:
        f.writelines("5\n")
        f.writelines("Pyramid Y below\n")
        f.writelines("O 0.0 0.0 0.0\n")
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[0][0], vertex[0][1], vertex[0][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[1][0], vertex[1][1], vertex[1][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[4][0], vertex[4][1], vertex[4][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[5][0], vertex[5][1], vertex[5][2]))

    # PYRAMID 5 -- X Above
    with open("ROTATIONS/pyramid_points_Xabove_{0:03d}.xyz".format(irot), 'w') as f:
        f.writelines("5\n")
        f.writelines("Pyramid X above\n")
        f.writelines("O 0.0 0.0 0.0\n")
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[0][0], vertex[0][1], vertex[0][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[3][0], vertex[3][1], vertex[3][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[4][0], vertex[4][1], vertex[4][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[7][0], vertex[7][1], vertex[7][2]))

    # PYRAMID 6 -- X Below
    with open("ROTATIONS/pyramid_points_Xbelow_{0:03d}.xyz".format(irot), 'w') as f:
        f.writelines("5\n")
        f.writelines("Pyramid Y below\n")
        f.writelines("O 0.0 0.0 0.0\n")
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[1][0], vertex[1][1], vertex[1][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[2][0], vertex[2][1], vertex[2][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[5][0], vertex[5][1], vertex[5][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[6][0], vertex[6][1], vertex[6][2]))

    # Write vertices
    with open("ROTATIONS/vertex_{0:03d}.xyz".format(irot), 'w') as f:
        f.writelines("9\n")
        f.writelines("Vertex\n")
        f.writelines("O 0.0 0.0 0.0\n")
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[0][0], vertex[0][1], vertex[0][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[1][0], vertex[1][1], vertex[1][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[2][0], vertex[2][1], vertex[2][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[3][0], vertex[3][1], vertex[3][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[4][0], vertex[4][1], vertex[4][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[5][0], vertex[5][1], vertex[5][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[6][0], vertex[6][1], vertex[6][2]))
        f.writelines("N {0:.3f} {1:.3f} {2:.3f}\n".format(vertex[7][0], vertex[7][1], vertex[7][2]))

    # ===================================================
    # Write run_pyramid.tcl

    txt="""
proc newRep { sel type color rep imol} {

    mol selection $sel
    mol representation $type
    mol addrep $imol
    mol showrep $imol $rep on
    mol modcolor $rep $imol $color

}

#Draw Arrow
proc vmd_draw_arrow {mol start end} {
  # an arrow is made of a cylinder and a cone
  set middle [vecadd $start [vecscale 0.9 [vecsub $end $start]]]
  graphics $mol color yellow
  graphics $mol cylinder $start $middle radius 0.10
  graphics $mol color blue
  graphics $mol cone $middle $end radius 0.25
  unset middle
}

# Sum two lists
proc sum_lists {l1 l2} {

    set result {}
    foreach x $l1 y $l2 {
        lappend result [expr {$x + $y}]
    }
    return $result
}

# Center of geometry
proc geom_center {selection} {
        # set the geometrical center to 0
        set gc [veczero]
        # [$selection get {x y z}] returns a list of {x y z}
        #    values (one per atoms) so get each term one by one
        foreach coord [$selection get {x y z}] {
           # sum up the coordinates
           set gc [vecadd $gc $coord]
        }
        # and scale by the inverse of the number of atoms
        return [vecscale [expr 1.0 /[$selection num]] $gc]
}

# Draw normal vector
proc draw_normal_vector {sel n1 index} {
   set c1 [ geom_center $sel ] 
   set c2 [sum_lists $c1 $n1]
   draw arrow $c1 $c2
   draw color black
   draw text $c1 " Normal $index" size 1 
}
    
    """

    txt += """
proc pyramid_Z_above {imol islabel} {
   # Normal 1
   set s1 [atomselect $imol "index 1 to 4"]
   # Normal 2
   set s2 [atomselect $imol "index 0 1 2"]
    # Normal 3
   set s3 [atomselect $imol "index 0 2 3"]
   # Normal 4
   set s4 [atomselect $imol "index 0 3 4"]
   # Normal 5
   set s5 [atomselect $imol "index 0 1 4"]
   """
    txt += "\n"
    txt += "   set n1 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[0][0][0], normals[0][0][1], normals[0][0][2])
    txt += "   set n2 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[0][1][0], normals[0][1][1], normals[0][1][2])
    txt += "   set n3 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[0][2][0], normals[0][2][1], normals[0][2][2])
    txt += "   set n4 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[0][3][0], normals[0][3][1], normals[0][3][2])
    txt += "   set n5 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[0][4][0], normals[0][4][1], normals[0][4][2])

    txt += """
   if {$islabel == 1} {
       draw_normal_vector $s1 $n1 "z1"
       draw_normal_vector $s2 $n2 "z2"
       draw_normal_vector $s3 $n3 "z3"
       draw_normal_vector $s4 $n4 "z4"
       draw_normal_vector $s5 $n5 "z5"
   }

}
"""

    txt += """
proc pyramid_Z_below {imol islabel} {
   # Normal 1
   set s1 [atomselect $imol "index 1 to 4"]
   # Normal 2
   set s2 [atomselect $imol "index 0 1 2"]
   # Normal 3
   set s3 [atomselect $imol "index 0 2 3"]
   # Normal 4
   set s4 [atomselect $imol "index 0 3 4"]
   # Normal 5
   set s5 [atomselect $imol "index 0 1 4"]
    """
    txt += "\n"
    txt += "   set n1 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[1][0][0], normals[1][0][1], normals[1][0][2])
    txt += "   set n2 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[1][1][0], normals[1][1][1], normals[1][1][2])
    txt += "   set n3 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[1][2][0], normals[1][2][1], normals[1][2][2])
    txt += "   set n4 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[1][3][0], normals[1][3][1], normals[1][3][2])
    txt += "   set n5 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[1][4][0], normals[1][4][1], normals[1][4][2])

    txt += """
   if {$islabel == 1} {
      draw_normal_vector $s1 $n1 "z6"
      draw_normal_vector $s2 $n2 "z7"
      draw_normal_vector $s3 $n3 "z8"
      draw_normal_vector $s4 $n4 "z9"
      draw_normal_vector $s5 $n5 "z10"
   }

}
    """

    txt += """
proc pyramid_X_above {imol islabel} {
   # Normal 1
   set s1 [atomselect $imol "index 1 to 4"]
   # Normal 2
   set s2 [atomselect $imol "index 0 1 2"]
   # Normal 3
   set s3 [atomselect $imol "index 0 1 3"]
   # Normal 4
   set s4 [atomselect $imol "index 0 3 4"]
   # Normal 5
   set s5 [atomselect $imol "index 0 2 4"]
    """
    txt += "\n"
    txt += "   set n1 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[2][0][0], normals[2][0][1], normals[2][0][2])
    txt += "   set n2 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[2][1][0], normals[2][1][1], normals[2][1][2])
    txt += "   set n3 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[2][2][0], normals[2][2][1], normals[2][2][2])
    txt += "   set n4 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[2][3][0], normals[2][3][1], normals[2][3][2])
    txt += "   set n5 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[2][4][0], normals[2][4][1], normals[2][4][2])

    txt += """
   if {$islabel == 1} {
      draw_normal_vector $s1 $n1 "x1"
      draw_normal_vector $s2 $n2 "x2"
      draw_normal_vector $s3 $n3 "x3"
      draw_normal_vector $s4 $n4 "x4"
      draw_normal_vector $s5 $n5 "x5"
   }

}
    """

    txt += """
proc pyramid_X_below {imol islabel} {
   # Normal 1
   set s1 [atomselect $imol "index 1 to 4"]
   # Normal 2
   set s2 [atomselect $imol "index 0 1 3"]
   # Normal 3
   set s3 [atomselect $imol "index 0 3 4"]
   # Normal 4
   set s4 [atomselect $imol "index 0 2 4"]
   # Normal 5
   set s5 [atomselect $imol "index 0 1 2"]
    """
    txt += "\n"
    txt += "   set n1 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[3][0][0], normals[3][0][1], normals[3][0][2])
    txt += "   set n2 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[3][1][0], normals[3][1][1], normals[3][1][2])
    txt += "   set n3 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[3][2][0], normals[3][2][1], normals[3][2][2])
    txt += "   set n4 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[3][3][0], normals[3][3][1], normals[3][3][2])
    txt += "   set n5 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[3][4][0], normals[3][4][1], normals[3][4][2])

    txt += """
   if {$islabel == 1} {
      draw_normal_vector $s1 $n1 "x6"
      draw_normal_vector $s2 $n2 "x7"
      draw_normal_vector $s3 $n3 "x8"
      draw_normal_vector $s4 $n4 "x9"
      draw_normal_vector $s5 $n5 "x10"
   }

}
    """

    txt += """
proc pyramid_Y_above {imol islabel} {
    # Normal 1
    set s1  [atomselect $imol "index 1 to 4"]
    # Normal 2
    set s2  [atomselect $imol "index 0 2 4"]
    # Normal 3
    set s3  [atomselect $imol "index 0 3 4"]
    # Normal 4
    set s4  [atomselect $imol "index 0 1 3"]
    # Normal 5
    set s5  [atomselect $imol "index 0 1 2"]
    """

    txt += "\n"
    txt += "   set n1 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[4][0][0], normals[4][0][1], normals[4][0][2])
    txt += "   set n2 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[4][1][0], normals[4][1][1], normals[4][1][2])
    txt += "   set n3 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[4][2][0], normals[4][2][1], normals[4][2][2])
    txt += "   set n4 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[4][3][0], normals[4][3][1], normals[4][3][2])
    txt += "   set n5 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[4][4][0], normals[4][4][1], normals[4][4][2])


    txt += """
   if {$islabel == 1} {
    draw_normal_vector $s1 $n1 "y1"
    draw_normal_vector $s2 $n2 "y2"
    draw_normal_vector $s3 $n3 "y3"
    draw_normal_vector $s4 $n4 "y4"
    draw_normal_vector $s5 $n5 "y5"
   }
}
    """

    txt += """
    proc pyramid_Y_below {imol islabel} {
    # Normal 1
    set s1  [atomselect $imol "index 1 to 4"]
    # Normal 2
    set s2  [atomselect $imol "index 0 1 2"]
    # Normal 3
    set s3  [atomselect $imol "index 0 1 3"]
    # Normal 4
    set s4  [atomselect $imol "index 0 3 4"]
    # Normal 5
    set s5  [atomselect $imol "index 0 2 4"]
    """

    txt += "\n"
    txt += "   set n1 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[5][0][0], normals[5][0][1], normals[5][0][2])
    txt += "   set n2 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[5][1][0], normals[5][1][1], normals[5][1][2])
    txt += "   set n3 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[5][2][0], normals[5][2][1], normals[5][2][2])
    txt += "   set n4 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[5][3][0], normals[5][3][1], normals[5][3][2])
    txt += "   set n5 [ list {0:.5f} {1:.5f} {2:.5f}]\n".format(normals[5][4][0], normals[5][4][1], normals[5][4][2])


    txt += """
   if {$islabel == 1} {
    draw_normal_vector $s1 $n1 "y6"
    draw_normal_vector $s2 $n2 "y7"
    draw_normal_vector $s3 $n3 "y8"
    draw_normal_vector $s4 $n4 "y9"
    draw_normal_vector $s5 $n5 "y10"
   }
}
    """

    txt += """\n
display projection orthographic
#axes location off
color Display Background white
display depthcue off

# Load Centers of Geometry ======================
set acog_ref ../acog_debug.xyz
mol new $acog_ref
set imol_ref [molinfo top]
mol delrep 0 $imol_ref
set rep1 0
newRep "all" "Points 10.0" "Name" $rep1 $imol_ref

# Load Trajectory ===============================
set trj_pdb ../anisotropy_full.pdb
mol new $trj_pdb
set imol_trj [molinfo top]
mol delrep 0 $imol_trj
set rep2 0
newRep "all" "CPK" "Name" $rep2 $imol_trj
#pbc box -center origin 
    """

    txt += """
 # Load Pyramid points ===============================
set pyr_xyz_za pyramid_points_Zabove_{0:03d}.xyz
mol new $pyr_xyz_za
set imol_pyr_za [molinfo top]
mol delrep 0 $imol_pyr_za
set rep3 0
newRep "all" "CPK" "Name" $rep3 $imol_pyr_za
mol modmaterial $rep3 $imol_pyr_za Transparent
topo addbond 0 1
topo addbond 0 2
topo addbond 0 3
topo addbond 0 4

# Load Pyramid points ===============================
set pyr_xyz_zb pyramid_points_Zbelow_{0:03d}.xyz
mol new $pyr_xyz_zb
set imol_pyr_zb [molinfo top]
mol delrep 0 $imol_pyr_zb
set rep4 0
newRep "all" "CPK" "Name" $rep4 $imol_pyr_zb
mol modmaterial $rep4 $imol_pyr_zb Transparent
topo addbond 0 1
topo addbond 0 2
topo addbond 0 3
topo addbond 0 4

# Load Pyramid points ===============================
set pyr_xyz_xa pyramid_points_Xabove_{0:03d}.xyz
mol new $pyr_xyz_xa
set imol_pyr_xa [molinfo top]
mol delrep 0 $imol_pyr_xa
set rep5 0
newRep "all" "CPK" "Name" $rep5 $imol_pyr_xa
mol modmaterial $rep5 $imol_pyr_xa Transparent
topo addbond 0 1
topo addbond 0 2
topo addbond 0 3
topo addbond 0 4

# Load Pyramid points ===============================
set pyr_xyz_xb pyramid_points_Xbelow_{0:03d}.xyz
mol new $pyr_xyz_xb
set imol_pyr_xb [molinfo top]
mol delrep 0 $imol_pyr_xb
set rep6 0
newRep "all" "CPK" "Name" $rep6 $imol_pyr_xb
mol modmaterial $rep6 $imol_pyr_xb Transparent
topo addbond 0 1
topo addbond 0 2
topo addbond 0 3
topo addbond 0 4

# Load Pyramid points ===============================
set pyr_xyz_ya pyramid_points_Yabove_{0:03d}.xyz
mol new $pyr_xyz_ya
set imol_pyr_ya [molinfo top]
mol delrep 0 $imol_pyr_ya
set rep7 0
newRep "all" "CPK" "Name" $rep7 $imol_pyr_ya
mol modmaterial $rep7 $imol_pyr_ya Transparent
topo addbond 0 1
topo addbond 0 2
topo addbond 0 3
topo addbond 0 4

# Load Pyramid points ===============================
set pyr_xyz_yb pyramid_points_Ybelow_{0:03d}.xyz
mol new $pyr_xyz_yb
set imol_pyr_yb [molinfo top]
mol delrep 0 $imol_pyr_yb
set rep7 0
newRep "all" "CPK" "Name" $rep7 $imol_pyr_yb
mol modmaterial $rep7 $imol_pyr_yb Transparent
topo addbond 0 1
topo addbond 0 2
topo addbond 0 3
topo addbond 0 4   
    """.format(irot)

    txt += """
# Load Vertex points ===============================
set vert vertex_{0:03d}.xyz
mol new $vert
set imol_vert [molinfo top]
mol delrep 0 $imol_vert
set rep8 0
newRep "all" "CPK" "Name" $rep8 $imol_vert
topo addbond 1 2
topo addbond 1 4
topo addbond 1 5
topo addbond 2 3
topo addbond 2 6
topo addbond 3 4
topo addbond 3 7
topo addbond 4 8
topo addbond 5 6
topo addbond 5 8
topo addbond 6 7
topo addbond 7 8
    """.format(irot)

    txt += """
# Pyramid Z above ================================
#pyramid_Z_above $imol_pyr_za 1
#pyramid_Z_below $imol_pyr_zb 1
#pyramid_X_above $imol_pyr_xa 1
#pyramid_X_below $imol_pyr_xb 1
pyramid_Y_above $imol_pyr_ya 1
pyramid_Y_below $imol_pyr_yb 1
"""

    with open("ROTATIONS/run_pyramid_{0:03d}.tcl".format(irot), 'w') as f:
        f.writelines(txt)

