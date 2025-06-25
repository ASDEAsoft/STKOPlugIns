'''
This file contains utility functions for parsing and processing sections in the MidasGenImporter plugin.
'''

from MidasGenImporter.document import section
from typing import Tuple, List

class _globals:
    # key: offest_type, value (y_factor, z_factor)
    offset_factors = { 
        "LT": (-1, 1),
        "CT": (0, 1),
        "RT": (1, 1),
        "LC": (-1, 0),
        "CC": (0, 0),
        "RC": (1, 0),
        "LB": (-1, -1),
        "CB": (0, -1),
        "RB": (1, -1),
    }
    # key: offset_type, value (y_factor, z_factor)
    u_movement = { 
        "LT": (1, -1),
        "LC": (1, 0),
        "LB": (1, 1),
        "CB": (0, 1),
        "RB": (-1, 1),
        "RC": (-1, 0),
        "RT": (-1, -1),
        "CT": (0, -1),
        "CC": (0, 0)
    }
    special_factor = {
        "LT": (1, 1),
        "LC": (1, 0),
        "LB": (1, 1),
        "CB": (0, 1),
        "RB": (1, 1),
        "RC": (1, 0),
        "RT": (1, 1),
        "CT": (0, 1),
        "CC": (0, 0)
    }

def _get_offset_anchor(shape: int, shape_data:List[float], offset_type: str) -> Tuple[float, float]:
    # note: now all supported shapes have the same offset anchor
    H, B = shape_data[0], shape_data[1]
    if offset_type == "RT":
        return (B / 2, H / 2)  # Right Top
    elif offset_type == "LT":
        return (-B / 2, H / 2)  # Left Top
    elif offset_type == "CT":
        return (0, H / 2)  # Center Top
    elif offset_type == "RC":
        return (B / 2, 0)  # Right Center
    elif offset_type == "CC":
        return (0, 0)  # Center Center
    elif offset_type == "LC":
        return (-B / 2, 0)  # Left Center
    elif offset_type == "RB":
        return (B / 2, -H / 2)  # Right Bottom
    elif offset_type == "LB":
        return (-B / 2, -H / 2)  # Left Bottom
    elif offset_type == "CB":
        return (0, -H / 2)  # Center Bottom
    else:
        raise ValueError(f"Unsupported offset type: {offset_type}")

def get_section_centroid(shape: int, shape_data:List[float]) -> Tuple[float, float]:
    """
    Returns the centroid coordinates (y, z) of a section based on its shape type.
    Note: returns the centroid wrt the geometric center of the section.
    Args:
        shape (int): The shape type identifier of the section.
        shape_data (List[float]): List containing shape data as per *.MGT file format.
    Returns:
        Tuple[float, float]: The centroid coordinates (y, z) of the section.
    Raises:
        ValueError: If the provided shape type is not supported.
    """
    if shape == section.shape_type.SB:
        return (0.0, 0.0)
    elif shape == section.shape_type.L:
        H, B, tw, tf = shape_data[:4]
        # compute the areas (same as T-section)
        A_flange = B * tf
        A_web = tw * (H - tf)
        A = A_flange + A_web
        # compute the static moments wrt the y-axis (horizontal) and z-axis (vertical)
        # knowing that the web is at left and the flange is at the top
        Sy_flange =  A_flange * (B / 2)
        Sz_flange = A_flange * (H - tf / 2)
        Sy_web = A_web * (tw / 2)
        Sz_web = A_web * ((H - tf) / 2)
        centroid_y = (Sy_flange + Sy_web) / A
        centroid_z = (Sz_flange + Sz_web) / A
        # wrt to the geometric center of the section
        centroid_y -= B / 2
        centroid_z -= H / 2
        return (centroid_y, centroid_z)
    elif shape == section.shape_type.T:
        H, B, tw, tf = shape_data[:4]
        # compute the centroid for a T-section knowing that the 
        # flange is a rectangle B x tf
        # and the web is a rectangle tw x (H - tf)
        # compute the areas
        A_flange = B * tf
        A_web = tw * (H - tf)
        A = A_flange + A_web
        # compute the static moments wrt the y-axis (horizontal) and z-axis (vertical)
        Sy_flange = A_flange * (B / 2)
        Sz_flange = A_flange * (H - tf / 2)
        Sy_web = A_web * (B / 2)
        Sz_web = A_web * ((H - tf) / 2)
        centroid_y = (Sy_flange + Sy_web) / A
        centroid_z = (Sz_flange + Sz_web) / A
        # wrt to the geometric center of the section
        centroid_y -= B / 2
        centroid_z -= H / 2
        return (centroid_y, centroid_z)
    else:
        raise ValueError(f"Unsupported section shape: {shape}")

def get_section_offset(shape : int, shape_data:List[str], offset_data:List[str]) -> Tuple[float, float]:
    """
    Calculates the offset of a section based on its shape, shape data, and offset data.
    Args:
        shape (int): The type of the section shape, typically an enum value from `frame_section.shape_type`.
        shape_data (List[str]): List of strings containing data specific to the section shape.
        offset_data (List[str]): List of 7 strings representing offset parameters:
            - offset_type (str): Type of offset.
            - icent (str): Centroid indicator.
            - iref (str): Reference indicator.
            - iHORZ (str): Horizontal offset indicator.
            - hUSER (str): User-defined horizontal offset.
            - iVERT (str): Vertical offset indicator.
            - vUSER (str): User-defined vertical offset.
    Returns:
        Tuple[float, float]: The calculated (horizontal, vertical) offset for the section.
    Raises:
        ValueError: If `offset_data` does not contain exactly 7 elements or if the section shape is unsupported.
    """
    
    # get centroid
    centroid_y, centroid_z = get_section_centroid(shape, shape_data)
    center_y, center_z = centroid_y, centroid_z

    # check offset data
    if len(offset_data) != 7:
        raise ValueError(f"Invalid offset data length: {len(offset_data)}. Expected 7 elements.")
    offset_type = offset_data[0]
    icent = float(offset_data[1])
    offset_reference = float(offset_data[2])
    horizontal_offset_presence = float(offset_data[3])
    horizontal_offset = float(offset_data[4])
    vertical_offset_presence = float(offset_data[5])
    vertical_offset = float(offset_data[6])

    if horizontal_offset_presence == 0 and vertical_offset_presence == 0: # caso 1
        ay, az =  _get_offset_anchor(shape, shape_data, offset_type)
        if icent == 1:
            ay = ay - center_y
            az = az - center_z 
        else:
            ay = ay - center_y * _globals.special_factor[offset_type][0]
            az = az - center_z * _globals.special_factor[offset_type][1]
    else:
        if horizontal_offset_presence == 0 and vertical_offset_presence == 1: 
            if offset_reference == 0: # caso 2
                ay, _ = _get_offset_anchor(shape, shape_data, offset_type)
                if icent == 1:
                    ay = ay - center_y
                    az = vertical_offset * _globals.offset_factors[offset_type][1]
                else:
                    ay = ay - center_y * _globals.special_factor[offset_type][0]
                    az = vertical_offset * _globals.offset_factors[offset_type][1]
            else: # caso 3
                ay, az = _get_offset_anchor(shape, shape_data, offset_type)
                if icent == 1:
                    ay = ay - center_y
                    az = az + vertical_offset * _globals.u_movement[offset_type][1] - center_z * _globals.special_factor[offset_type][1]
                else:
                    ay = ay - center_y * _globals.special_factor[offset_type][0]
                    az = az + vertical_offset * _globals.u_movement[offset_type][1] - center_z * _globals.special_factor[offset_type][1]
        
        elif horizontal_offset_presence == 1 and vertical_offset_presence == 0: 
            if offset_reference == 0: # caso 4
                _, az = _get_offset_anchor(shape, shape_data, offset_type)
                if icent == 1:
                    ay = horizontal_offset * _globals.offset_factors[offset_type][0]
                    az = az - center_z
                else:
                    ay = horizontal_offset * _globals.offset_factors[offset_type][0]
                    az = az - center_z * _globals.special_factor[offset_type][1]
            else: # caso 5
                ay, az = _get_offset_anchor(shape, shape_data, offset_type)
                if icent == 1:
                    ay = ay + horizontal_offset * _globals.u_movement[offset_type][0] - center_y * _globals.special_factor[offset_type][0]
                    az = az - center_z
                else:
                    ay = ay + horizontal_offset * _globals.u_movement[offset_type][0] - center_y * _globals.special_factor[offset_type][0]
                    az = az - center_z * _globals.special_factor[offset_type][1]
                    
        elif horizontal_offset_presence == 1 and vertical_offset_presence == 1: 
            if offset_reference == 0: # caso 6
                ay, az = _get_offset_anchor(shape, shape_data, offset_type)
                if icent == 1:
                    ay = horizontal_offset * _globals.offset_factors[offset_type][0]
                    az = vertical_offset * _globals.offset_factors[offset_type][1]
                else:
                    ay = horizontal_offset * _globals.offset_factors[offset_type][0]
                    az = vertical_offset * _globals.offset_factors[offset_type][1]
            else: # caso 7
                ay, az = _get_offset_anchor(shape, shape_data, offset_type)
                if icent == 1:
                    ay = horizontal_offset * _globals.u_movement[offset_type][0] - center_y * _globals.special_factor[offset_type][0]
                    az = vertical_offset * _globals.u_movement[offset_type][1] - center_z * _globals.special_factor[offset_type][1]
                else:
                    ay = ay + horizontal_offset * _globals.u_movement[offset_type][0] - center_y * _globals.special_factor[offset_type][0]
                    az = az + vertical_offset * _globals.u_movement[offset_type][1] - center_z * _globals.special_factor[offset_type][1]

    return (ay, az)

def test_T_section():
    from matplotlib import pyplot as plt
    H = 400.0
    B = 200.0
    tw = 10.0
    tf = 80.0
    shape_data = [H, B, tw, tf]
    centroid_y, centroid_z = get_section_centroid(section.shape_type.T, shape_data)
    # create y,z points for the T-section wrt the center
    Y = [-tw/2, tw/2, tw/2, B/2, B/2, -B/2, -B/2, -tw/2, -tw/2]
    Z = [-H/2, -H/2, H/2-tf, H/2-tf, H/2, H/2, H/2-tf, H/2-tf, -H/2]
    # plot the T-section
    plt.figure(figsize=(6, 6))
    plt.plot(Y, Z, marker='o')
    plt.plot([centroid_y], [centroid_z], marker='x', color='red', markersize=10, label='Centroid')
    plt.axis('equal')
    plt.show()

def test_L_section():
    from matplotlib import pyplot as plt
    H = 400.0
    B = 200.0
    tw = 10.0
    tf = 10.0
    shape_data = [H, B, tw, tf]
    centroid_y, centroid_z = get_section_centroid(section.shape_type.L, shape_data)
    # create y,z points for the T-section wrt the center
    Y = [-B/2, -B/2+tw, -B/2+tw, B/2, B/2, -B/2, -B/2]
    Z = [-H/2, -H/2, H/2-tf, H/2-tf, H/2, H/2, -H/2]
    # plot the T-section
    plt.figure(figsize=(6, 6))
    plt.plot(Y, Z, marker='o')
    plt.plot([centroid_y], [centroid_z], marker='x', color='red', markersize=10, label='Centroid')
    plt.axis('equal')
    plt.show()

def test_section(shape:int, shape_info:List[float], offset_data:List[str]):
    from matplotlib import pyplot as plt
    centroid_y, centroid_z = get_section_centroid(shape, shape_info)
    offset_y, offset_z = get_section_offset(shape, shape_info, offset_data)
    # Visualization can be added here if needed
    H, B, tw, tf = shape_info[:4]
    if shape == section.shape_type.T:
        Y = [-tw/2, tw/2, tw/2, B/2, B/2, -B/2, -B/2, -tw/2, -tw/2]
        Z = [-H/2, -H/2, H/2-tf, H/2-tf, H/2, H/2, H/2-tf, H/2-tf, -H/2]
    elif shape == section.shape_type.L:
        Y = [-B/2, -B/2+tw, -B/2+tw, B/2, B/2, -B/2, -B/2]
        Z = [-H/2, -H/2, H/2-tf, H/2-tf, H/2, H/2, -H/2]
    else:
        return None
    plt.figure(figsize=(6, 6))
    plt.plot(Y, Z, marker='o')
    plt.plot([centroid_y], [centroid_z], marker='x', color='red', markersize=10, label='Centroid')
    plt.plot([offset_y], [offset_z], marker='o', color='green', markersize=10, label='Offset')
    plt.axis('equal')
    plt.show()
    print(offset_y, offset_z)
