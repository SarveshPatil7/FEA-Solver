# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 12:05:07 2025

@author: sarve
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.tri as tri

#DO NOT MODIFY THIS FUNCTION
def plot_fea_mesh(fig, n_elem, elem_node, x, y, linecolor='k'):
    ax = fig.gca()
    for ELEMENT_ID in range(n_elem):
        node_list = elem_node[ELEMENT_ID,:]
        XEL = x[node_list-1] # Element x coordinates
        YEL = y[node_list-1] # Element y coordinates
        X = np.append(XEL,x[node_list[0]-1]) # Close the drawing loop by appending first node at end
        Y = np.append(YEL,y[node_list[0]-1]) # Close the drawing loop by appending first node at end
        # Plot the  element
        plt.plot(X,Y,linecolor,linewidth=0.5) # Check, was ax.plot

    plt.draw()

#DO NOT MODIFY THIS FUNCTION
def plot_Heat2D_results(tur, n_elem, elem_node, coords, field='temperature', show_plot=False,filename='temperature_contour.png'):

    '''
    Functions assumes plotting for Q4 elements
    Inputs are as follows:
    * tur: (n_node, 1) array of temperatures
    * n_elem: number of elements in the mesh
    * elem_node: (n_elem, 4) array of element nodes
    * coords; (n_node, 2) array of nodal coordinates in the mesh
    * field: 'temperature' or 'flux' or 'gradT'
    '''

    def quads_to_tris(elem_node):
        '''
        Cuts each Q4 element into 2 triangular element

        Returns list of triangle elements and their node numbers
        '''
        tris = [[None for j in range(3)] for i in range(2*len(elem_node))]
        for i in range(len(elem_node)):
            j = 2*i # For for access
            # List of the node numbers of the Q4 element
            # --> would work for straight-sided Q8 as well
            n0 = elem_node[i][0]-1
            n1 = elem_node[i][1]-1
            n2 = elem_node[i][2]-1
            n3 = elem_node[i][3]-1
            # First triangle
            tris[j][0] = n0
            tris[j][1] = n1
            tris[j][2] = n2
            # Second triangle
            tris[j + 1][0] = n2
            tris[j + 1][1] = n3
            tris[j + 1][2] = n0
        return tris

    # Nodal coordinates
    x=coords[:,0]
    y=coords[:,1]

    # 'z' is the parameter to plot
    if(field=='temperature'):
        z=tur
    elif(field=='flux'):
        z=...
    elif(field=='gradT'):
        z=...
    else:
        print('Error: this plotting component is not supported')
        return

    # Initial figure setup
    plotHeat = plt.figure()
    if(field=='temperature'):
        plottitle = 'Temperature solution'
    elif(field=='flux'):
        plottitle = 'Heat flux solution'
    elif(field=='flux'):
        plottitle = 'Temperature gradient solution'
    plt.axes().set_aspect('equal')
    plt.title(plottitle)
    plt.xlabel('x')
    plt.ylabel('y')

    # Plot undeformed mesh
    plot_fea_mesh(plotHeat, n_elem, elem_node, x, y, linecolor='k')

    # Triangulate the domain
    tris = quads_to_tris(elem_node) # Convert elements to triangles

    # Unstructured triangular grid (undeformed mesh)
    triangulation = tri.Triangulation(x, y, tris) 

    # Plot the contours
    levels = np.linspace(min(z), max(z), 13) # '13' matches Abaqus colorbar
    plt.tricontourf(triangulation, z, cmap='rainbow', levels=levels)
    cbar = plt.colorbar()
    if(field=='temperature'):
        cbar.ax.set_ylabel('Temperature', rotation=90)
    elif(field=='flux'):
        cbar.ax.set_ylabel('Heat Flux', rotation=90)
    # Use scientific notation
    cbar.ax.ticklabel_format(scilimits=[0,0]) 

    # Show the plot
    plt.savefig(filename)
    if show_plot:
        plt.show()

    return(plotHeat)



