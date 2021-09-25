# -*- coding: utf-8 -*-
# ChemTools is a collection of interpretive chemical tools for
# analyzing outputs of the quantum chemistry calculations.
#
# Copyright (C) 2016-2019 The ChemTools Development Team
#
# This file is part of ChemTools.
#
# ChemTools is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# ChemTools is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
"""Output Module.

This module contains functions for generating scripts for visualizing
purposes using ChimeraX.
"""
import os
import numpy as np

__all__ = ['print_cx_script_nci', 'print_cx_script_isosurface',
           'print_cx_script_multiple_cube']


def _cx_script_start():
    """Generate part of the beginning part of the CX script.

    Returns
    -------
    CX script responsible for beginning of VMD script
    """
    return ('#!/usr/local/bin/chimerax\n'
            '# CX script written by save_state $Revision: 1.41 $\n'
            # '# VMD version: 1.8.6\n'
            # 'set viewplist\n'
            # 'set fixedlist\n'
            '#\n'
            '# Display settings\n'
            '# Publication preset includes white background and black silhouettes\n'
            'windowsize 700 500\n'
            'camera ortho\n'
            'preset pub\n'
            # 'color Element {C} gray\n'
            # 'color Element {Cl} green\n'
            # 'axes location Off\n'
            # 'light 2 on\n'
            # 'light 3 on\n'
            '#\n')


def _cx_script_molecule(representation, *mol_files):
    """Generate part of the CX script that loads the molecule information.

    Parameters
    ----------
    representation: str
        Representation of a molecule. Choose either 'ball' or 'stick'.
    mol_files : str
        Names of the input files that represent the molecule
        .mol2, .pdb and .cube files are supported
        Cube files can correspond to the density, reduced density gradient, isosurface, color of
        isosurface, etc

    Returns
    -------
    Part of the CX script that constructs the molecule

    Raises
    ------
    ValueError
        If no files are given
    TypeError
        If unsupported file type (i.e. not mol2, pdb or cube)
    """
    output = '# load new molecule\n'
    if len(mol_files) == 0:
        raise ValueError('Need at least one molecule file')
    for i, mol in enumerate(mol_files):
        ext = os.path.splitext(mol)[1]
        if ext == '.mol2':
            file_type = 'mol2'
        elif ext == '.pdb':
            file_type = 'pdb'
        elif ext == '.cube':
            file_type = 'cube'
        else:
            raise TypeError('Unsupported file type, {0}'.format(ext))

        output += 'open {0} format {1}\n'.format(mol, file_type)
    output += ('#\n'
               '# representation of the atoms\n')
    if representation.lower() == 'ball':
        output += ('style ball\n'
                'size ballScale 0.3\n')
    elif representation.lower() == 'stick':
        output += ('style stick\n'
                'size stickRadius 0.2\n')
    else:
        raise ValueError('Argument representation is not recognized.')
    output += ('color #1 byelement\n'
               '#\n')
    return output

def _cx_script_isosurface(isosurf=0.5, index=1, style_type='surface', scalemin=-0.05, scalemax=0.05, colorscheme='blue', negative=False):
    """Generate part of the CX script that configures the isosurface.

    Parameters
    ----------
    isosurf : float
        Iso-surface value to plot.
    index : int, optional
        Index of the file that contains the iso-surface data (in the order loaded by
        `_cx_script_molecule`)
    style_type : str, optional
        Option that controls how the iso-surface will be drawn
        One of 'surface', 'mesh' or 'solid'.
    scalemin : float, optional
        Smallest value to color on the iso-surface.
    scalemax : float, optional
        Largest value to color on the iso-surface.
    colorscheme : str, optional
        Color scheme used in CX script for coloring the iso-surface.
        The iso-surface can be colored with just one color. For this option, a string is needed
        specifying the color (X11 color names from the CSS3 specification).

        If a cube file is used to color the isosurface, the str must be one of the following options corresponding to a color gradient
            ============  =====================================
            Options       Description
            ============  =====================================
            'rainbow'     small=blue, middle=green, large=red
            'redblue'     small=red, middle=white, large=blue
            'bluered'     small=blue, middle=white, large=red
            'cyanmaroon'  small=cyan, middle=white, large=maroon
            'grayscale'   small=black, large=white
            'rgb'*        small=red, middle=green, large=blue
            ===========   =====================================
            * Note that `rgb` is not one of the predefined palettes in CX.
        Check the program or website for more details on colors and sequence of colors (or palettes) specifications.
    negative : bool, optional
        Determines if you want to plot the negative of the iso-surface as well.

    Returns
    -------
    Part of the ChimeraX script that constructs the isosurface

    Raises
    ------
    TypeError
        If `isosurf` is not a float
        If `index` is not an integer
        If `style_type` is not one of 'surface', 'mesh' or 'solid'
        If `scalemin` is not float
        If `scalemax` is not float
        If `colorscheme` is not a string or one of 'rainbow', 'redblue', 'bluered', 'cyanmaroon', 'grayscale', 'rgb' when `index` is 2.

    """
    if not isinstance(isosurf, float):
        raise TypeError('`isosurf` must be a float')

    if not isinstance(index, int):
        raise TypeError('`index` must be an integer')

    if style_type not in ['surface', 'mesh', 'solid']:
        raise TypeError('Unsupported `style_type`. Must be one of {0}, {1}  or {2}'
                        ''.format('surface', 'mesh', 'solid'))

    if not isinstance(scalemin, float):
        raise TypeError('`scalemin` must be a float')
    if not isinstance(scalemax, float):
        raise TypeError('`scalemax` must be a float')

    palettes = ['rainbow', 'redblue', 'bluered', 'cyanmaroon', 'grayscale', 'rgb']
    if not isinstance(colorscheme, str) and not hasattr(colorscheme, '__iter__'):
        raise TypeError('Unsupported colorscheme, {0}'.format(colorscheme))
    if index==2 and colorscheme not in palettes:
        raise TypeError('Unsupported colorscheme. Must be one of {0}, {1}, {2}, {3}, {4} or {5}'
                        ''.format('rainbow', 'redblue', 'bluered', 'cyanmaroon', 'grayscale', 'rgb'))

    # set color for positive and negative iso surfaces
    if not negative and isinstance(colorscheme, (str)) and (colorscheme not in palettes):
        pos_color, neg_color = colorscheme, colorscheme
    elif negative and hasattr(colorscheme, '__iter__') and len(colorscheme) == 2:
        pos_color, neg_color = colorscheme
    else:
        pos_color, neg_color = ['blue', 'red']

    output = '# add representation of the surface\n'
    if negative:
        output += ('volume #{index} level -{isosurf:.5f} color {neg_color}'
        ' level {isosurf:.5f} color {pos_color} style {style}\n'.format(index=index, isosurf=isosurf, pos_color=pos_color, neg_color=neg_color, style=style_type))
    elif not negative and index == 1:
        output += 'volume #{index} level {isosurf:.5f} style {style}\n'.format(isosurf=isosurf, index=index, style=style_type)
        output += 'volume #{index} color {color}\n'.format(index=index, color=pos_color)
    else:
        output += 'volume #{index} level {isosurf:.5f} style {style}\n'.format(isosurf=isosurf, index=index, style=style_type)
        output += 'color sample #{index1} map #{index2} offset {offset}'.format(index1=index, index2=index-1, offset=0)
        if colorscheme.lower() == 'rgb':
            output += ' palette {min},#ff0000:0,#00ff00:{max},#0000ff\n'.format(min=scalemin, max=scalemax)
        else:
            output += (' palette {scheme} range {min},{max}\n'.format(min=scalemin, max=scalemax))
    output += '#\n'
    return output


def print_cx_script_nci(scriptfile, densfile, rdgfile, isosurf=0.5, denscut=0.05):
    r"""Generate CX (ChimeraX) script for visualizing NCI isosurfaces.

    Visualizes NCI (non-covalent interactions) isosurfaces subject to the constraint of
    density(r) < denscut, i.e. low-density, and colored based on the
    sign(:math:`\lambda_2`) :math:`\rho`.

    Parameters
    ----------
    scriptfile : str
        Name of VMD script file to generate.
    densfile : str
        Name of density cube file.
    rdgfile : str
        Name of reduced density gradient cube file.
    isosurf : float, optional
        Reduced density gradient iso-surface to visualize.
    denscut : float, optional
        Density cutoff used in creating reduced density gradient cube file.
        Similar to NCIPlot program, reduced density gradient of points with
        density > denscut will be set to 100.0 to display reduced density gradient
        iso-surface subject to the constraint of low density.

    Note
    ----
    The script is the same as the one generated by NCIPlot software version 1.0.
    """
    output = _cx_script_start()
    output += _cx_script_molecule('ball', densfile, rdgfile)
    output += _cx_script_isosurface(isosurf=isosurf, scalemin=-denscut, scalemax=denscut, index=1, colorscheme='rainbow')
    with open(scriptfile, 'w') as f:
        f.write(output)

def print_cx_script_isosurface(scriptfile, isofile, colorfile=None, isosurf=0.5,
                                scalemin=-0.05, scalemax=0.05, colorscheme=None, negative=False,
                                representation='ball'):
    """Generate CX (ChimeraX) script for visualizing the isosurface.

    Visualize isosurface based on one cube file when coloring by the value of another cube file on the isosurface.

    Parameters
    ----------
    scriptfile : str
        Name of CX script file to generate.
    isofile : str
        Name of cube file used in CX script for visualizing the isosurface.
    colorfile : str, optional
        Name of cube file used in CX script for coloring the isosurface.
        If None, the isofile is used for coloring.
    isosurf : float, optional
        The value of the isosurface to visualize used in CX script.
    scalemin : float, optional
        Smallest value to color on the iso-surface used in CX script.
    scalemax : float, optional
        Largest value to color on the iso-surface used in CX script.
    colorscheme : str, optional
        Color scheme used in CX script for coloring the iso-surface.
        The iso-surface can be colored with just one color. For this option, a string is needed
        specifying the color (X11 color names from the CSS3 specification).

        If a cube file is used to color the isosurface, the str must be one of the following options corresponding to a color gradient
            ============  =====================================
            Options       Description
            ============  =====================================
            'rainbow'     small=blue, middle=green, large=red
            'redblue'     small=red, middle=white, large=blue
            'bluered'     small=blue, middle=white, large=red
            'cyanmaroon'  small=cyan, middle=white, large=maroon
            'grayscale'   small=black, large=white
            'rgb'*        small=red, middle=green, large=blue
            ===========   =====================================
            * Note that `rgb` is not one of the predefined palettes in CX.
        Check the program or website for more details on colors and sequence of colors (or palettes) specifications.
    negative : bool, optional
        Determines if you want to plot the negative of the iso-surface as well.
    representation: str, optional
        Representation of a molecule.

    """
    # set default color schemes
    if colorscheme is None:
        if colorfile is None:
            colorscheme = 'blue'
        else:
            colorscheme = 'rgb'

    output = _cx_script_start() # write script header
    if colorfile is not None: # load surface files and possibly also coordinates
        output += _cx_script_molecule(representation, colorfile, isofile)
        file_index = 2
    else:
        output += _cx_script_molecule(representation, isofile)
        file_index = 1

    output += _cx_script_isosurface(isosurf=isosurf, index=file_index,
                                     scalemin=scalemin, scalemax=scalemax, colorscheme=colorscheme, negative=negative)

    with open(scriptfile, 'w') as f:
        f.write(output)


def print_cx_script_multiple_cube(scriptfile, cubes, isosurfs=None,
                                   scalemin=-0.05, scalemax=0.05, colors=None,
                                   representation='ball'):
    """Generate CX (ChimeraX) script for visualizing multiple cube files.

    Visualize multiple cube files simultaneously where data from each cube file is colored
    differently.

    Parameters
    ----------
    scriptfile : str
        Name of CX script file to generate.
    cubes : list of str
        Names of cube files to plot
    isosurfs : float, list of float
        Isovalue at which the plot (iso-surface) is generated
        If a float is given, then this is the value of iso-surface for all cube files
        Default value is 0.5 for all iso-surfaces
    scalemin : float
        Smallest value to color on the iso-surface used in VMD script.
        Default is -0.05
    scalemax : float
        Largest value to color on the iso-surface used in VMD script.
        Default is 0.05
    colors : list of str
        Colors of each cube file data
        See CX program or manual for details.

    Raises
    ------
    TypeError
        If cube files are not provided as a list or tuple
        If colors are not provided as a list or tuple of the same length as the cube files
    ValueError
        If any of the cube files cannot be found
        If any of the colors are not a string
    NotImplementedError
        If no colors are specified
    """
    if not isinstance(cubes, (list, tuple)):
        raise TypeError('The cube files must be given as a list or tuple')
    if not all(os.path.isfile(cube) for cube in cubes):
        raise ValueError('Cannot find at least one of the cube files')

    if isosurfs is None:
        isosurfs = [0.5 for i in cubes]
    elif isinstance(isosurfs, float):
        isosurfs = [isosurfs for i in cubes]
    if not (isinstance(isosurfs, (list, tuple)) and len(isosurfs) == len(cubes)):
        raise TypeError('The isosurfs must be provided as a list or tuple of same length as the '
                        'number of cube files')
    elif not all(isinstance(isosurf, float) for isosurf in isosurfs):
        raise TypeError('Each iso-surface value must be a float')

    if colors is None:
        # colors = range(len(cubes))
        raise NotImplementedError('A list of colors must be specified')
    elif not (isinstance(colors, (list, tuple)) and len(colors) == len(cubes)):
        raise TypeError('The colors must be provided as a list or tuple of the same length as the '
                        'number of cube files')
    elif not all(isinstance(color, str) for color in colors):
        raise ValueError('Each color must be given as a string')

    output = _cx_script_start()
    output += _cx_script_molecule(representation, *cubes)
    for i, (isosurf, color) in enumerate(zip(isosurfs, colors)):
        output += _cx_script_isosurface(isosurf=isosurf, index=i, scalemin=scalemin,
                                         scalemax=scalemax, colorscheme=color)

    with open(scriptfile, 'w') as f:
        f.write(output)
