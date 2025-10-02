# -*- coding: utf-8 -*-
"""
Created on Mon Apr 28 12:03:01 2025

@author: sarve
"""

import numpy as np
import pandas as pd
#DO NOT MODIFY THIS FUNCTION
def parse_abaqus_input(file_path, print_mesh_data=False):

    '''
    Read an Abaqus input file for 
    1. a 2D static linear elasticity problem
    2. a 2D thermal problem
    Known restrictions (may not be exhaustive):
    * Only one element type allowed in the mesh
    * Only one part/part instance allowed in the assembly
    * Loads supported for linear elasticty:
      - Concentrated nodal load (*Cload)
      - Constant traction (*Dsload, load_type = TRVEC)
      - Constant body load (*Dload)
    * Loads supported for thermal problems:
      - Constant surface flux (*Dsflux)
      - Constant body load (*Dflux)
      - Convection with constant properties (*Sfilm)
    * Material parameters supported for linear elasticity:
      - Modulus and Poisson ratio (elastic)
      - Plane stress/plane strain thickness
    * Material properties supported for thermal problems:
      - Thermal conductivity
      - Thickness in 3rd dimension
    * Assumptions to work with FEA code:
      - Nodes numbered starting with 1, no skipping
      - Elements numbered starting with 1, no skipping
    '''

    # Initialize data structures to store parsed information
    instances = []
    nodes = {'numbers': [], 'coordinates': []}
    elements = {'element_numbers': [], 'nodes': []}
    part_node_sets = {}
    assembly_node_sets = {}
    part_element_sets = {}
    assembly_element_sets = {}
    surfaces = {}
    materials = {}
    boundaries = [] 
    concentrated_loads = [] 
    surface_loads = [] 
    surface_fluxes = []
    body_loads = {}
    body_thermal_loads = []
    convection_fluxes = []
    sections = [] 
    element_type = None

    # Initialize variables to track the current parsing context
    current_section = None
    material_name = None
    current_set_name = None
    is_generate = False
    reading_elastic = False
    reading_thermal = False
    reading_solid_section = False
    reading_instance = False
    reading_part = False
    reading_assembly = False

    # Open the input file and read it line by line
    with open(file_path, 'r') as file:
        for line in file:

            '''
            Reminder of common string operations
            * line.strip(): Remove spaces at beginning and end of string
            * line.split(','): Split string into lists using comma 
                               (in this case) as separator
            '''

            # Remove any space at beginning and end of string
            line = line.strip()
            # Skip empty lines and comments
            if not line or line.startswith('**'):
                continue
            # Check for new section headers
            if line.startswith('*'):
                line_data = line.split(',')
                header = line_data[0]

                '''
                In this first set of conditionals, the header line is identified
                and useful info is extracted. Relevant data will be read
                in the subsequent loop iteration, which is handled in the second
                block of conditional statements.

                You can add additional input file read capability by defining
                the behavior for a given header here
                '''

                # Identify the type of section and update the current context
                # Need a few special treatments when reading data in the Instance section
                # So make a flag that indiates when we begin and stop reading info on an instance
                if header == '*Instance':
                    reading_instance = True
                elif header =='*End Instance':
                    reading_instance = False
                elif header == '*Part':
                    reading_part = True
                elif header =='*End Part':
                    reading_part = False
                elif header == '*Assembly':
                    reading_assembly = True
                elif header =='*End Assembly':
                    reading_assembly = False
                
                # Reading node data
                elif header == '*Node':
                    current_section = 'nodes'
                # Reading element data
                elif header == '*Element':
                    element_type = line_data[1].split('=')[1]
                    current_section = 'elements'
                # Reading a node set
                elif header.startswith('*Nset'):
                    current_set_name = line_data[1].split('=')[1]
                    # Part or assembly set
                    if(reading_instance or reading_part):
                        part_node_sets[current_set_name] = []
                    else:
                        assembly_node_sets[current_set_name] = []                        
                    is_generate = 'generate' in line_data[-1]
                    current_section = 'nset'
                # Reading an element set
                elif header.startswith('*Elset'):
                    current_set_name = line_data[1].split('=')[1]
                    # Part or assembly set
                    if(reading_instance or reading_part):
                        part_element_sets[current_set_name] = []
                    else:
                        assembly_element_sets[current_set_name] = []
                    is_generate = 'generate' in line_data[-1]
                    current_section = 'elset'
                # Reading a surface data
                elif header.startswith('*Surface'):
                    surface_name = line_data[2].split('=')[1]
                    surface_type = line_data[1].split('=')[1]
                    surfaces[surface_name] = {'type': surface_type, 
                                              'element_sets': [], 
                                              'element_faces': []}
                    current_section = 'surface'
                # Reading material data
                elif header == '*Material':
                    material_name = line_data[1].split('=')[1]
                    materials[material_name] = {}
                    current_section = 'material'
                # Reading elastic material properties
                elif header == '*Elastic':
                    current_section = 'elastic'
                    reading_elastic = True
                # Reading thermal material properties
                elif header == '*Conductivity':
                    current_section = 'conductivity'
                    reading_thermal = True
                # Reading density
                elif header == '*Density':
                    current_section = 'density'
                    reading_thermal = True
                # Reading specific heat
                elif header == '*Specific Heat':
                    current_section = 'specific_heat'
                    reading_thermal = True
                # Reading essential BC info
                elif header == '*Boundary':
                    current_section = 'boundary'
                # Concentrated force
                elif header == '*Cload':
                    current_section = 'cload'
                # Constant traction
                elif header == '*Dsload':
                    current_section = 'dsload'
                # Constant flux
                elif header == '*Dsflux':
                    current_section = 'dsflux'
                # Solid section data
                elif header == '*Solid Section':
                    elset = line_data[1].split('=')[1]
                    material = line_data[2].split('=')[1]
                    # Next line sets thickness to a default of 1.0 if not specified
                    sections.append({'part_element_set': elset, 'material_name': material, 'thickness': 1.0})
                    current_section = 'solid_section'
                    reading_solid_section = True
                # Body load for structural problem
                elif header == '*Dload':
                    current_section = 'dload'
                # Body load for thermal problem
                elif header == '*Dflux':
                    current_section = 'dflux'
                # Convection BC for thermal problem
                elif header == '*Sfilm':
                    current_section = 'sfilm'
                else:
                    current_section = None
                continue

            # Split the line into line_data
            line_data = [p.strip() for p in line.split(',') if p.strip()]

            '''
            This 2nd block of conditional statements is set up to read all of the
            data associated with a particular section in the code. The header line
            was read in the previous iteration of the loop
            '''

            if current_section == 'instance':
                instances.append['instance_name']
            
            # Process data according to the current section
            if current_section == 'nodes':
                
                # Node number
                node_id = int(line_data[0]) 

                # For the coordinates code below:
                # * line_data[1:] is nodal coordinates in string format
                # * map(float, ...) converts each string to a float, stored in an iterable
                # list(...) convert the iterable to a list
                coordinates = list(map(float, line_data[1:]))

                # The node numbers and coordinates are appended to their lists
                # in the dictionary "nodes"
                nodes['numbers'].append(node_id)
                nodes['coordinates'].append(coordinates)

            elif current_section == 'elements':
                
                # Convert the line into a list of integers
                line_data = list(map(int, line_data))

                # Element number is the first entry in the line
                element_number = line_data[0]

                # Remaining numbers are the node numbers
                node_numbers = line_data[1:]

                # Append to relevant lists in the dictionary "elements"
                elements['element_numbers'].append(element_number)
                elements['nodes'].append(node_numbers)

            elif current_section == 'nset':

                # Convert line data into integers
                line_data = list(map(int, line_data))

                # If line ends in 'generate', create node list using 'range' (start, stop, increment)
                if is_generate and len(line_data) == 3:
                    if(reading_instance or reading_part):
                        part_node_sets[current_set_name].extend(range(line_data[0], line_data[1] + 1, line_data[2]))
                    else:
                        assembly_node_sets[current_set_name].extend(range(line_data[0], line_data[1] + 1, line_data[2]))
                # If line does not end in 'generate', the data is just a list of all the nodes
                else:
                    if(reading_instance or reading_part):
                        part_node_sets[current_set_name].extend(line_data)
                    else:
                        assembly_node_sets[current_set_name].extend(line_data)

            elif current_section == 'elset':
                
                # Convert line data into integers
                line_data = list(map(int, line_data))

                # If line ends in 'generate', create node list using 'range' (start, stop, increment)
                if is_generate and len(line_data) == 3:
                    if(reading_instance or reading_part):
                        part_element_sets[current_set_name].extend(range(line_data[0], line_data[1] + 1, line_data[2]))
                    else:
                        assembly_element_sets[current_set_name].extend(range(line_data[0], line_data[1] + 1, line_data[2]))
                       
                # If line does not end in 'generate', the data is just a list of all the nodes
                else:
                    if(reading_instance or reading_part):
                        part_element_sets[current_set_name].extend(line_data)
                    else:
                        assembly_element_sets[current_set_name].extend(line_data)
                    
            elif current_section == 'surface':
                
                # Store the name of the element set
                surfaces[surface_name]['element_sets'].append(line_data[0])
                # Store the name of the face for all elements in this set
                surfaces[surface_name]['element_faces'].append(line_data[1].strip())

                # if line_data[0].isdigit():
                #     surfaces[surface_name]['elements'].append(int(line_data[0]))
                #     surfaces[surface_name]['face_number'].append(line_data[1].strip())
                # else:
                #     surfaces[surface_name]['elements'].append(line_data[0])
                #     surfaces[surface_name]['face_number'].append(line_data[1].strip())


            elif current_section == 'material':
                continue  # Material type is read in the next line

            elif current_section == 'elastic' and reading_elastic:
                line_data = list(map(float, line_data))
                materials[material_name]['E']  = line_data[0]  # Young's modulus
                materials[material_name]['nu'] = line_data[1] # Poisson's ratio
                reading_elastic = False  # Reset the flag after reading

            elif current_section == 'conductivity' and reading_thermal:
                line_data = list(map(float, line_data))
                materials[material_name]['k'] = line_data[0]  # Thermal conductivity
                reading_thermal = False  # Reset the flag after reading

            elif current_section == 'density' and reading_thermal:
                line_data = list(map(float, line_data))
                materials[material_name]['rho'] = line_data[0]  # Density
                reading_thermal = False  # Reset the flag after reading

            elif current_section == 'specific_heat' and reading_thermal:
                line_data = list(map(float, line_data))
                materials[material_name]['c'] = line_data[0]  # Specific heat
                reading_thermal = False  # Reset the flag after reading

            elif current_section == 'boundary':
                assembly_node_set = line_data[0]
                # For x-symmetry, U1 = 0
                if(line_data[1] == 'XSYMM'):
                    dof_id_start = 1
                    dof_id_end   = 1
                # For y-symmetry, U2 = 0
                elif(line_data[1] == 'YSYMM'):
                    dof_id_start = 2
                    dof_id_end   = 2
                # For general listing of fixed DOF ids
                else:
                    dof_id_start = int(line_data[1])
                    dof_id_end   = int(line_data[2])
                # Allow for non-zero prescribed displacment
                if len(line_data) == 4:
                    dof_val = float(line_data[3])
                else:
                    # DOF is zero-valued if not prescirbed
                    dof_val = 0.
                boundaries.append({'assembly_node_set': assembly_node_set, 
                                   'dof_id_start': dof_id_start, 
                                   'dof_id_end': dof_id_end,
                                   'dof_val': dof_val})
                
            # Concentrated loads
            elif current_section == 'cload':
                assembly_node_set = line_data[0]
                dof_id = int(line_data[1])
                load_value = float(line_data[2])
                # concentrated_loads.append([node_set, dof_id, load_value])
                concentrated_loads.append({'assembly_node_set': assembly_node_set, 
                                           'dof_id': dof_id, 
                                           'load_value': load_value})
                
            # Constant surface traction for structural problems
            elif current_section == 'dsload':
                surface_name = line_data[0]
                load_type = line_data[1]
                load_magnitude = float(line_data[2])
                direction_vector = list(map(float, line_data[3:]))
                surface_loads.append({'surface_name': surface_name,
                                      'load_type': load_type,
                                      'load_magnitude': load_magnitude,
                                      'direction_vector': direction_vector})
                
            # Constant surface flux for thermal problems
            elif current_section == 'dsflux':
                surface_name = line_data[0]
                flux_type = line_data[1]
                flux_val = float(line_data[2])
                surface_fluxes.append({'surface_name': surface_name,
                                       'flux_type': flux_type,
                                       'flux_val': flux_val})  
                
            # Convection flux for thermal problems
            elif current_section == 'sfilm':
                surface_name = line_data[0]
                load_type = line_data[1]
                Tinf = float(line_data[2])
                h    = float(line_data[3])
                convection_fluxes.append({'surface_name': surface_name,
                                          'flux_type': load_type,
                                          'Tinf': Tinf,
                                          'h': h}) 

            # Solid section thickness in 3rd dimension
            elif current_section == 'solid_section' and reading_solid_section:
                if line_data:
                    try:
                        thickness = float(line_data[0])
                    except ValueError:
                        thickness = 1.0
                    sections[-1]['thickness'] = thickness
                reading_solid_section = False

            # Body load for structural problems
            elif current_section == 'dload':
                element_set = line_data[0] # Name of element set that gets the load
                component = line_data[1]   # BX or BY
                load = float(line_data[2]) # Load value (per unit volume, e.g. rho*g)

                # # Add this dictionary to the list of body loads
                # body_loads.append({'element_set': element_set, 'load': [BX, BY]})

                if element_set not in body_loads:
                    body_loads[element_set] = {'BX': 0.0, 'BY': 0.0}
                if component == 'BX':
                    body_loads[element_set]['BX'] = load
                elif component == 'BY':
                    body_loads[element_set]['BY'] = load

            # Body load for thermal problems
            elif current_section == 'dflux':
                element_set = line_data[0] # Name of element set that gets the load
                load_type = line_data[1]   # Type of load, must be BF
                load = float(line_data[2]) # Load value (per unit volume)

                body_thermal_loads.append({'element_set': element_set,
                                           'load_type': load_type,
                                           'load': load}) 


    # Convert lists to numpy arrays for easier numerical processing
    nodes['coordinates'] = np.array(nodes['coordinates'])
    elements['nodes'] = np.array(elements['nodes'])

    # Print the parsed data for verification
    if(print_mesh_data):
        print('=========================================================')
        print('Abaqus mesh data from input file: ' + file_path)
        print('=========================================================')
        print("Nodes:\n", nodes)
        print("Elements:\n", elements)
        print("Part Node Sets:\n", part_node_sets)
        print("Assembly Node Sets:\n", assembly_node_sets)
        print("Part Element Sets:\n", part_element_sets)
        print("Assembly Element Sets:\n", assembly_element_sets)
        print("Surfaces:\n", surfaces)
        print("Materials:\n", materials)
        print("Boundaries:\n", boundaries)
        print("Concentrated Loads:\n", concentrated_loads)
        print("Surface Loads:\n", surface_loads)
        print("Surface Fluxes:\n", surface_fluxes)
        print("Convection Fluxes:\n", convection_fluxes)
        print("Body Loads:\n", body_loads)
        print("Body Thermal Loads:\n", body_thermal_loads)
        print("Sections:\n", sections)
        print('=========================================================')
        print('End of mesh data')
        print('=========================================================')

    return {
        'nodes': nodes,
        'elements': {
            'type': element_type,
            'element_numbers': np.array(elements['element_numbers']),
            'nodes': elements['nodes']
        },
        'part_node_sets': part_node_sets,
        'assembly_node_sets': assembly_node_sets,
        'part_element_sets': part_element_sets,
        'assembly_element_sets': assembly_element_sets,
        'surfaces': surfaces,
        'materials': materials,
        'boundaries': boundaries,
        'concentrated_loads': concentrated_loads,
        'surface_loads': surface_loads,
        'surface_fluxes': surface_fluxes,
        'convection_fluxes': convection_fluxes,
        'body_loads': body_loads,
        'body_thermal_loads': body_thermal_loads,
        'sections': sections
    }

#DO NOT MODIFY THIS FUNCTION
def make_thermal_problem_data(filename):

    mesh_data = parse_abaqus_input(filename)

    # Number of nodes and elements
    N_NODE = len(mesh_data['nodes']['numbers'])
    N_ELEM = len(mesh_data['elements']['element_numbers'])

    # Determine parameters that depend on element type
    element_type = mesh_data['elements']['type']

    # Q4 heat transfer element with full integration
    if(element_type == 'DC2D4'):
        NNODE_ELE = 4
        ShapeOrder = 1
        Ng = 2
    # This element type is not currently supported
    else:
        print(element_type + ' is not currently supported')
        return

    # COORDS: (N_NODE, 2) numpy array of coordinates
    COORDS = np.array(mesh_data['nodes']['coordinates'])

    # Element number (id) for each element: Must be numbered 1 to N_ELEM
    ELEM_ID = np.array(mesh_data['elements']['element_numbers'], dtype=int)

    # Node list for each element
    ELEM_NODE = np.array(mesh_data['elements']['nodes'], dtype=int)

    # ELEM_PROPS contains [k, thickness, rho, c]
    ELEM_PROPS = np.zeros((N_ELEM, 4))

    # ELEM_BLOAD contains [Q]
    ELEM_BLOAD = np.zeros(N_ELEM)

    # Loop on elements
    for id in ELEM_ID:
        # Loop on sections to assign material properties and thickness
        for section in mesh_data['sections']:
            element_set = section['part_element_set'] # Name of element set associated with this section
            # If the elemnet number (id) is in the elemnet set, then assign these section properties
            if id in mesh_data['part_element_sets'][element_set]:
                material = section['material_name']
                ELEM_PROPS[id-1,0] = mesh_data['materials'][material]['k']
                ELEM_PROPS[id-1,1] = section['thickness']
                ELEM_PROPS[id-1,2] = mesh_data['materials'][material]['rho']
                ELEM_PROPS[id-1,3] = mesh_data['materials'][material]['c']
        
        # Loop on body loads to assign body load values
        for bload in mesh_data['body_thermal_loads']:

            print(bload)

            element_set = bload['element_set']
            Q = bload['load']

            if id in mesh_data['assembly_element_sets'][element_set]:
                ELEM_BLOAD[id-1] = Q


    # T_NODE contains a list of nodes and DOF ids that are fixed
    T_NODE = []
    T_VAL  = []
    # Loop on boundaries
    for boundary in mesh_data['boundaries']:

        # For thermal problems, the dof_id is always 1
        # Note it is stored as 11 in Abaqus

        node_set_name = boundary['assembly_node_set']
        dof_val       = boundary['dof_val']

        # If node_set_name begins with "Part", then assume it is a single node
        # The node number is at the end
        # This assumes a model with 1 "Part"
        # Should confirm this is a consistent behavior/format
        if(node_set_name[:4] == 'Part'):
            node_number = int(node_set_name.split('.')[1]) # Node number is 'x' in Part-1-1.x
            T_NODE.append(node_number)
            T_VAL.append(dof_val)

        # Otherwise, node_set_name indicates a node set with an accompanying list of nodes
        else:
            # Loop through the nodes in the node set
            for node_number in mesh_data['assembly_node_sets'][node_set_name]:
                T_NODE.append(node_number)
                T_VAL.append(dof_val)

    # Make T_NODE a numpy array of integers
    T_NODE = np.array(T_NODE, dtype=int)
    T_VAL  = np.array(T_VAL, dtype=float)
    # Sort T_NODE on nodes, then on DOF id
    if(len(T_NODE) > 0):
        indices = np.argsort(T_NODE)
        T_NODE = T_NODE[indices]
        T_VAL  = T_VAL[indices]
    N_PRE_T = len(T_NODE) # Returns number of rows

    # # Create FORCE_NODE and FORCE_VAL arrays
    # FORCE_NODE = []
    # FORCE_VAL  = []
    # for cload in mesh_data['concentrated_loads']:
        
    #     node_set_name = cload['node_set']
    #     dof_id = cload['dof_id']
    #     load_value  = cload['load_value']

    #     # Loop on nodes in the node set
    #     for node_number in mesh_data['node_sets'][node_set_name]:
    #         FORCE_NODE.append([node_number, dof_id])
    #         FORCE_VAL.append(load_value)

    # # Make FORCE_NODE and FORCE_VAL numpy arrays
    # FORCE_NODE = np.array(FORCE_NODE, dtype=int)
    # FORCE_VAL = np.array(FORCE_VAL, dtype=float)
    # # Sort FORCE_NODE and FORCE_VAL on nodes first, dof id second
    # if(len(FORCE_NODE > 0)):
    #     indices = np.lexsort((FORCE_NODE[:, 1], FORCE_NODE[:, 0]))
    #     FORCE_NODE = FORCE_NODE[indices]
    #     FORCE_VAL= FORCE_VAL[indices]
    # N_CLOAD = len(FORCE_NODE) # Returns number of rows


    # Loop on surface and convection fluxes
    FLUX_q_ELE = []    
    FLUX_q_VAL = []
    for surface_flux in mesh_data['surface_fluxes']:
        # Get the surface name
        surface_name = surface_flux['surface_name']
        # Get the surface associated with the load
        surface = mesh_data['surfaces'][surface_name]
        # Reader assumes constant value
        flux_val = surface_flux['flux_val']
        # Loop on the element sets associated with that surface
        for index, element_set in enumerate(surface['element_sets']):
            # Get all the elements for the current element set
            elements = mesh_data['assembly_element_sets'][element_set]
            # Get the face number: same for all elements in this set
            # Store face as an integer, using the 2nd character in the face name
            face_label = surface['element_faces'][index]
            face = int(face_label[1]) # Last Character in face label; assumes label is formatted e.g. 'S2'
            # Append the data to FLUX_q_ELE
            for element in elements:
                FLUX_q_ELE.append([element, face])
                FLUX_q_VAL.append(flux_val)
    FLUX_c_ELE = []
    FLUX_c_VAL = []
    for convection_flux in mesh_data['convection_fluxes']:
        # Get the surface name
        surface_name = convection_flux['surface_name']
        # Get the surface associated with the load
        surface = mesh_data['surfaces'][surface_name]
        # Reader assumes constant values
        Tinf = convection_flux['Tinf']
        h    = convection_flux['h']
        # Loop on the element sets associated with that surface
        for index, element_set in enumerate(surface['element_sets']):
            # Get all the elements for the current element set
            elements = mesh_data['assembly_element_sets'][element_set]
            # Get the face number: same for all elements in this set
            # Store face as an integer, using the 2nd character in the face name
            face_label = surface['element_faces'][index]
            face = int(face_label[1]) # Last Character in face label; assumes label is formatted e.g. 'S2'
            # Append the data to FLUX_q_ELE
            for element in elements:
                FLUX_c_ELE.append([element, face])
                FLUX_c_VAL.append([h, Tinf])

    # Convert the lists into numpy arrays
    FLUX_q_ELE = np.array(FLUX_q_ELE, dtype=int)
    FLUX_q_VAL = np.array(FLUX_q_VAL, dtype=float)
    N_FLUX_q   = len(FLUX_q_ELE) # Returns number of rows
    FLUX_c_ELE = np.array(FLUX_c_ELE, dtype=int)
    FLUX_c_VAL = np.array(FLUX_c_VAL, dtype=float)
    N_FLUX_c   = len(FLUX_c_ELE) # Returns number of rows

    # Print all outputs
    #print('N_NODE = ', N_NODE)
    #print('N_ELEM = ', N_ELEM)
    #print('NNODE_ELE = ', NNODE_ELE)
    #print('ShapeOrder = ', ShapeOrder)
    #print('Ng = ', Ng)
    #print('N_CLOAD = ', N_CLOAD)
    #print('N_PRE_T = ', N_PRE_T)
    #print('N_FLUX_q = ', N_FLUX_q)
    #print('N_FLUX_c = ', N_FLUX_c)

    #print('COORDS = \n', COORDS)
    #print('ELEM_NODE = \n', ELEM_NODE)

    #print('ELEM_PROPS = \n', ELEM_PROPS)
    #print('ELEM_BLOAD = \n', ELEM_BLOAD)

    #print('T_NODE = \n', T_NODE)
    #print('T_VAL = \n', T_VAL)

    #print('FORCE_NODE = \n', FORCE_NODE)
    # # print('FORCE_VAL = \n', FORCE_VAL)

    #print('FLUX_q_ELE = \n', FLUX_q_ELE)
    #print('FLUX_q_VAL = \n', FLUX_q_VAL)
    #print('FLUX_c_ELE = \n', FLUX_c_ELE)
    #print('FLUX_c_VAL = \n', FLUX_c_VAL)

    return(N_NODE, N_ELEM, NNODE_ELE, ShapeOrder, Ng, N_PRE_T, 
           N_FLUX_q, N_FLUX_c,COORDS, ELEM_NODE, ELEM_PROPS, ELEM_BLOAD, 
           T_NODE, T_VAL,FLUX_q_ELE, FLUX_q_VAL,FLUX_c_ELE, FLUX_c_VAL)


# Shape functions for isoparametric Q4 element
def shape2d_Q4(r):
    r1=r[0]
    r2=r[1]

    Nhat =  0.25*np.array([(1-r1)*(1-r2), (1+r1)*(1-r2), (1+r1)*(1+r2), (1-r1)*(1+r2) ])
    DNhat = 0.25*np.array([[-(1-r2), (1-r2), (1+r2), -(1+r2)] , 
                           [-(1-r1), -(1+r1) ,(1+r1) ,(1-r1)]])

    return Nhat, DNhat

def EQ_NUM_LM(N_NODE, N_ELEM, ELEM_NODE, T_NODE, N_PRE_T):
    EQ_NUM = np.zeros(N_NODE, dtype=int)

    # Negative for prescribed nodes
    for i in range(1, N_PRE_T + 1):
        EQ_NUM[T_NODE[i - 1] - 1] = -i
        
    # Positive for free nodes
    row = 0
    for i in range(N_NODE):
        if EQ_NUM[i] == 0:
            row += 1
            EQ_NUM[i] = row
            
    LM = np.zeros((N_ELEM, 4), dtype=int)
    for e in range(N_ELEM):
        for a in range(4): 
            node = ELEM_NODE[e, a]
            LM[e, a] = EQ_NUM[node - 1]
            
    N_DOF = row
    TOTAL_DOF = N_NODE
    #print('LM : ', LM)
    return EQ_NUM, LM, N_DOF, TOTAL_DOF

def Assemble_Tp(N_PRE_T, T_NODE, T_VAL, EQ_NUM):
    Tp = np.zeros(N_PRE_T)

    for i in range(N_PRE_T):
        node = T_NODE[i]
        val = T_VAL[i]
        row = -EQ_NUM[node - 1] - 1
        Tp[row] = val

    return Tp

def Kke(k, L1, L2):
        
    K1 = ((k * L2) / (6 * L1)) * np.array([[ 2, -2, -1,  1],
                                           [-2,  2,  1, -1],
                                           [-1,  1,  2, -2],
                                           [ 1, -1, -2,  2]])

    K2 = ((k * L1) / (6 * L2)) * np.array([[ 2,  1, -1, -2],
                                           [ 1,  2, -2, -1],
                                           [-1, -2,  2,  1],
                                           [-2, -1,  1,  2]])

    Ke = K1 + K2
    
    return Ke

def assemble_PQ(N_ELEM, ELEM_NODE, COORDS, ELEM_PROPS, ELEM_BLOAD, LM, N_TOTAL_DOF):

    P = np.zeros(N_TOTAL_DOF)

    for e in range(N_ELEM):
        nodes = ELEM_NODE[e] - 1
        xe = COORDS[nodes]
        q = ELEM_BLOAD[e]
        t = ELEM_PROPS[e, 1]

        gp = [-1/np.sqrt(3), 1/np.sqrt(3)]
        w = [1, 1]
        Pe = np.zeros((4, 1))
        
        for i in range(2):
            for j in range(2):
                r = np.array([gp[i], gp[j]])
                Nhat, DNhat = shape2d_Q4(r)
                J = DNhat @ xe
                detJ = np.linalg.det(J)
                Pe += Nhat[:, None] * q * t * detJ * w[i] * w[j]
        #Pe = np.array([[1],[1],[1],[1]]) * (q * t * area / 4)  # equally distributed to 4 nodes

        for a in range(4):
            row = LM[e, a]
            if row > 0:
                P[row - 1] += Pe[a]
                
    #print("P:")
    #print(P)

    return P

def assemble_Pq(N_FLUX_q, FLUX_q_ELE, FLUX_q_VAL, ELEM_NODE, COORDS, LM, N_TOTAL_DOF):

    P = np.zeros(N_TOTAL_DOF)

    # Face mapping based on local node
    face_nodes = {
        1: [0, 1],
        2: [1, 2],
        3: [2, 3],
        4: [3, 0]
    }

    for i in range(N_FLUX_q):
        e, face = FLUX_q_ELE[i]
        q = FLUX_q_VAL[i]
        e -= 1  # 0 index
        n1_local, n2_local = face_nodes[face]

        # Global node IDs
        nodes = ELEM_NODE[e] - 1
        coords = COORDS[nodes]

        # Edge length
        L_edge = np.linalg.norm(coords[n2_local] - coords[n1_local])

        Pe = q * L_edge / 2 * np.array([1.0, 1.0])

        for local_node, p in zip([n1_local, n2_local], Pe):
            P[nodes[local_node]] += p
        
    #print("P:")
    #print(P)

    return P

def Kce(h, L_edge):

    Ke_edge = (h * L_edge / 6.0) * np.array([[2.0, 1.0],
                                             [1.0, 2.0]])

    return Ke_edge

def assemble_K(N_ELEM, ELEM_NODE, COORDS, ELEM_PROPS, FLUX_c_ELE, FLUX_c_VAL, N_TOTAL_DOF):
    
    face_nodes = {1: (0, 1),
                  2: (1, 2),
                  3: (2, 3),
                  4: (3, 0)}
    
    K = np.zeros((N_TOTAL_DOF, N_TOTAL_DOF))
    
    for i in range(N_ELEM):
        nodes = ELEM_NODE[i] - 1  
        coords = COORDS[nodes]
        k = ELEM_PROPS[i, 0]
        
        L1 = np.linalg.norm(coords[1] - coords[0])
        L2 = np.linalg.norm(coords[3] - coords[0])
        Ke = Kke(k, L1, L2)
        
        for (elem_id, edge_id), (h, T_inf) in zip(FLUX_c_ELE, FLUX_c_VAL):
            if elem_id == i + 1:
                n1, n2 = face_nodes[edge_id]
                L_edge = np.linalg.norm(coords[n2] - coords[n1])
                Kc_edge = Kce(h, L_edge)
                # inject into the 4×4 Ke
                Ke[np.ix_([n1, n2], [n1, n2])] += Kc_edge
        
        for a in range(4):
            A = nodes[a]
            for b in range(4):
                B = nodes[b]
                K[A, B] += Ke[a, b]
    
    return K
    
def assemble_Pc(N_FLUX_c, FLUX_c_ELE, FLUX_c_VAL, ELEM_NODE, COORDS, LM, N_TOTAL_DOF, ELEM_PROPS):

    P_conv = np.zeros(N_TOTAL_DOF)

    face_nodes = {
        1: [0, 1],
        2: [1, 2],
        3: [2, 3],
        4: [3, 0]
    }

    for i in range(N_FLUX_c):
        e, face = FLUX_c_ELE[i]
        h, T_inf = FLUX_c_VAL[i]
        e -= 1
        n1_local, n2_local = face_nodes[face]
        nodes = ELEM_NODE[e] - 1
        coords = COORDS[nodes]
        L_edge = np.linalg.norm(coords[n2_local] - coords[n1_local])

        Pe = h * T_inf * L_edge / 2 * np.array([1.0, 1.0])
        
        for local_node, p in zip([n1_local, n2_local], Pe):
            P_conv[nodes[local_node]] += p

    #print('P_conv = ', P_conv)
    return P_conv

def assemble_partition(K, PQ, Pq, Pc, T_NODE, T_VAL, N_TOTAL_DOF, EQ_NUM, N_NODE):

    Kur = K
    Pur = PQ + Pq + Pc
    
    prescribed = [i for i in range(N_NODE) if EQ_NUM[i] < 0]
    free = [i for i in range(N_NODE) if EQ_NUM[i] > 0]

    KPP = Kur[np.ix_(prescribed, prescribed)]
    KPF = Kur[np.ix_(prescribed, free)]
    KFP = Kur[np.ix_(free, prescribed)]
    KFF = Kur[np.ix_(free, free)]

    PP = np.zeros(len(prescribed))
    for i, node in enumerate(prescribed):
        PP[i] = T_VAL[-EQ_NUM[node] - 1]
    PF = Pur[free]
    
    #print("KFF:")
    #print(KFF)
    #print("PF:")
    #print(PF)

    UP = T_VAL
    UF = np.zeros_like(free, dtype=float)

    return KPP, KPF, KFP, KFF, PP, PF, UP, UF, free, prescribed

def solve_and_assemble(KFF, KFP, PF, UP, free, prescribed, N_TOTAL_DOF):

    #print("KFF:")
    #print(KFF)
    #print('det KFF = ', np.linalg.det(KFF))
    
    UF = np.linalg.solve(KFF, PF - KFP @ UP)


    TUR = np.zeros(N_TOTAL_DOF)
    TUR[free] = UF
    TUR[prescribed] = UP
    
    return TUR

# Final solver function used in main.py
def solve_steady(N_NODE, N_ELEM, NNODE_ELE, ShapeOrder, Ng, N_PRE_T, N_FLUX_q, N_FLUX_c,COORDS, ELEM_NODE, ELEM_PROPS, ELEM_BLOAD, T_NODE, T_VAL,FLUX_q_ELE, FLUX_q_VAL,FLUX_c_ELE, FLUX_c_VAL):

    EQ_NUM, LM, N_DOF, N_TOTAL_DOF = EQ_NUM_LM(N_NODE, N_ELEM, ELEM_NODE, T_NODE, N_PRE_T)

    K = assemble_K(N_ELEM, ELEM_NODE, COORDS, ELEM_PROPS, FLUX_c_ELE, FLUX_c_VAL, N_TOTAL_DOF)
    PQ = assemble_PQ(N_ELEM, ELEM_NODE, COORDS, ELEM_PROPS, ELEM_BLOAD, LM, N_TOTAL_DOF)
    Pq = assemble_Pq(N_FLUX_q, FLUX_q_ELE, FLUX_q_VAL, ELEM_NODE, COORDS, LM, N_TOTAL_DOF)
    
    Pc = assemble_Pc(N_FLUX_c, FLUX_c_ELE, FLUX_c_VAL, ELEM_NODE, COORDS, LM, N_TOTAL_DOF, ELEM_PROPS)

    KPP, KPF, KFP, KFF, PP, PF, UP, UF, free, prescribed = assemble_partition(
        K, PQ, Pq, Pc, T_NODE, T_VAL, N_TOTAL_DOF, EQ_NUM, N_NODE)

    TUR = solve_and_assemble(KFF, KFP, PF, UP, free, prescribed, N_TOTAL_DOF)
    
    return TUR, EQ_NUM, LM, K, KPP, KPF, KFP, KFF, PP, PF, UP, UF

    
# Transient
def capacitance(N_ELEM, ELEM_NODE, COORDS, ELEM_PROPS, N_TOTAL_DOF):

    from math import sqrt
    gp = [-1/sqrt(3), 1/sqrt(3)]
    w  = [1.0, 1.0]
    C = np.zeros((N_TOTAL_DOF, N_TOTAL_DOF))

    for e in range(N_ELEM):
        nodes = ELEM_NODE[e] - 1
        coords = COORDS[nodes]
        k, t, rho, c = ELEM_PROPS[e]           
        Ce = np.zeros((4,4))

        for i in range(2):
            for j in range(2):
                r, s = gp[i], gp[j]
                Nhat, _ = shape2d_Q4((r,s))
                J = _ @ coords   
                detJ = np.linalg.det(J)
                Ce += rho*c*t * np.outer(Nhat, Nhat) * detJ * w[i]*w[j]

        for a in range(4):
            A = nodes[a]
            for b in range(4):
                B = nodes[b]
                C[A,B] += Ce[a,b]
    
    diag = np.sum(C, axis=1)
    C_lump = np.diag(diag)
    #print('C_lump = ', C_lump)
    return C_lump

def crit_step(K, C_lump):

    CiK = np.linalg.inv(C_lump) @ K
    eigs = np.linalg.eigvals(CiK)
    λmax = np.max(eigs.real)
    return 2.0 / λmax
    
def general_theta_solver(N_NODE, N_ELEM, ELEM_NODE, COORDS, ELEM_PROPS, 
                         T_NODE, T_VAL, FLUX_c_ELE, FLUX_c_VAL, N_FLUX_c,
                         N_FLUX_q, FLUX_q_ELE, FLUX_q_VAL, N_PRE_T, factor):
    
    EQ_NUM, LM,_,_ = EQ_NUM_LM(N_NODE, N_ELEM, ELEM_NODE, T_NODE, N_PRE_T)
    K = assemble_K(N_ELEM, ELEM_NODE, COORDS, ELEM_PROPS, FLUX_c_ELE, FLUX_c_VAL, N_NODE)
    Pc = assemble_Pc(N_FLUX_c, FLUX_c_ELE, FLUX_c_VAL, ELEM_NODE, COORDS, LM, N_NODE, ELEM_PROPS)
    Pq = assemble_Pq(N_FLUX_q, FLUX_q_ELE, FLUX_q_VAL, ELEM_NODE, COORDS, LM, N_NODE)
    
    C_lump = capacitance(N_ELEM, ELEM_NODE, COORDS, ELEM_PROPS, N_NODE)  
    crit_t = crit_step(K, C_lump)
    dt = crit_t * factor    #number adjusted according to question
    #dt = 50                # use this step for crank nicoleson
    theta = 0               # 0 for forward euler,  0.5 for crank nicoleson
    sol_t = 5000
    
    print('crit_t = ', crit_t)
    
    N = len(EQ_NUM)

    
    Kbar = theta * K + C_lump / dt
    
    M = -(1 - theta) * K + C_lump / dt
    
    Pur = Pq + Pc
    
    free = sorted([i for i in range(N_NODE) if EQ_NUM[i] > 0],
                  key=lambda i: EQ_NUM[i])
    prescribed = sorted([i for i in range(N_NODE) if EQ_NUM[i] < 0],
                  key=lambda i: -EQ_NUM[i])
    
    Kbar_ff = Kbar[np.ix_(free, free)]
    Kbar_fp = Kbar[np.ix_(free, prescribed)]
    M_ff    = M[np.ix_(free, free)]
    M_fp    = M[np.ix_(free, prescribed)]    
    invK_ff = np.linalg.inv(Kbar_ff)
    
    n_steps = int(np.ceil(sol_t / dt))
    times   = np.arange(n_steps+1) * dt
    TUR    = np.zeros((N, n_steps + 1))
    TUR[:, 0] = 140.0
    Tp0 = dict(zip([n-1 for n in T_NODE], T_VAL))

    if prescribed:
        T0_free = Tp0[prescribed[0]]
    else:
        T0_free = 0.0
    for i in prescribed:
        TUR[i, 0] = Tp0[i]
    for i in free:
        TUR[i, 0] = T0_free
        
    for n in range(n_steps):
        Tn = TUR[:, n]
    
        Pf_bar = M_ff.dot(Tn[free]) + M_fp.dot(Tn[prescribed]) + Pur[free]
               
        Tfn1 = invK_ff.dot(Pf_bar - Kbar_fp.dot(Tn[prescribed]))
        
        for i, idx in enumerate(free):
            TUR[idx, n+1] = Tfn1[i]
        for idx in prescribed:
            TUR[idx, n+1] = TUR[idx, n] 

    return times, TUR