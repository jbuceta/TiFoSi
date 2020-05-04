#!/usr/bin/env python
# coding: utf-8

# In[26]:


#%pylab inline

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.patches import Patch
import sys
import argparse
import re
import os

#################
# Some default settings

matplotlib.rcParams["figure.figsize"] = [4*6.4, 4*4.8]

#https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html

parser = argparse.ArgumentParser(description='Plot frames from TiFoSi simulation.')
parser.add_argument('--cmap', dest='cmap', default="terrain", help='Color map for protein levels.', choices = plt.colormaps())
parser.add_argument('--cmapt', dest='cmapType', default="tab20", help='Color map for cell type.', choices = plt.colormaps())
parser.add_argument('--oformat', dest='oformat', default="png", help='Output images file format.', choices = ["png", "eps", "svg", "jpg"])
parser.add_argument('--dpi', dest='dpi', default=300, help='Output images resolution.', type=int)
parser.add_argument('--title', dest='title', default="plot", help='String prepended to the file name.', type=str)
parser.add_argument('--nolegend', dest='colorbar', default=True, action='store_false', help='Remove legend in plots.')
parser.add_argument('--elw', dest='edgelinewidth', default=2, help='Edge line width.', type=int)
parser.add_argument('--elc', dest='edgecolor', default="#FFFFFF", help='Edge line color.', type=str)

parser.add_argument('--nopframes', dest='movieprotein', default=True, action='store_false', help='Do not compute protein levels frames.')
parser.add_argument('--notframes', dest='moviecelltypes', default=True, action='store_false', help='Do not compute cell type frames.')

parser.add_argument('--ifolder', dest='ifolder', default="./", help='Folder containing input files.', required=True)
parser.add_argument('--ofolder', dest='ofolder', default="./", help='Folder to save output files.', required=True)

globals().update(parser.parse_args().__dict__)

cmap = matplotlib.cm.get_cmap(cmap)
cmapType = matplotlib.cm.get_cmap(cmapType)
if not re.search(r'^#(?:[0-9a-fA-F]{3}){1,2}$', edgecolor):
    sys.exit("Error: The --elc parameter [ %s ] is not a valid color. The color has to be specified in by using a hex color code (e.g. '#FFFFFF' for white or '#FF0000' for red)." % edgecolor)
if not os.path.isdir(ifolder):
    sys.exit("Error: The folder [ %s ] does not exist." % ifolder)
if not os.path.isdir(ofolder):
    sys.exit("Error: The folder [ %s ] does not exist." % ofolder)

#################
# Auxiliar function

def updateBBox(bbox1, bbox2):
    return {"xmin": min(bbox1["xmin"], bbox2["xmin"]),
            "ymin": min(bbox1["ymin"], bbox2["ymin"]),
            "xmax": max(bbox1["xmax"], bbox2["xmax"]),
            "ymax": max(bbox1["ymax"], bbox2["ymax"])}

def updatePBBox(bbox1, bbox2):
    result = {}
    for p in PNAMES:
        result[p] = ( min(bbox1[p][0], bbox2[p][0]), max(bbox1[p][1], bbox2[p][1]) )
    return result

#################
# Cell and Frame definitions

class Cell:
    def __init__(self, data):
        self.id = data[0]
        self.type = int(data[1])
        self.nvertexs = int(data[2])
        self.nproteins = int(data[4])
        self.proteins = {}
        for ip, p in enumerate(PNAMES):
            self.proteins[p] = float(data[5 + ip])
        
        shift = 5 + len(PNAMES) + 2 + self.nvertexs
        self.vertexs = []
        for iv in range(self.nvertexs):
            self.vertexs.append( (float(data[shift + iv*2]), float(data[shift + iv*2 + 1])) )

    def getCellBoundingBox(self):
        bbox = None
        for v in self.vertexs:
            if bbox == None:
                bbox = {"xmin": v[0], "ymin": v[1], "xmax": v[0], "ymax": v[1]}
            bbox = {"xmin": min(v[0], bbox["xmin"]),
                    "ymin": min(v[1], bbox["ymin"]),
                    "xmax": max(v[0], bbox["xmax"]),
                    "ymax": max(v[1], bbox["ymax"])}
        return bbox
            
class Frame:
    def __init__(self, data):
        if len(data) != 2:
            sys.exit("Error: Malformed file.")
        
        self.ncells = int(data[0])
        self.cells = []
        
        self.bbox = None
        self.proteinBBox = None
    
        self.ctypes = None
    
    def addCell(self, cell):
        self.cells.append(cell)
        
    def getFrameBoundingBox(self):
        self.bbox = None
        for c in self.cells:
            if self.bbox == None:
                self.bbox = c.getCellBoundingBox()
            self.bbox = updateBBox(self.bbox, c.getCellBoundingBox())
            
    def getProteinBoundingBox(self):
        self.proteinBBox = {}
        for p in PNAMES:
            temp = []
            for c in self.cells:
                temp.append(c.proteins[p])
            self.proteinBBox[p] = (min(temp), max(temp))
    
    def getCellTypes(self):
        self.ctypes = []
        for c in self.cells:
            self.ctypes.append(c.type)
        self.ctypes = set(self.ctypes)

#################
# Read data

print("Loading data...")

try:
    fileproteinnames = "protein_order.dat"
    f = open(os.path.join(ifolder, fileproteinnames))
    PNAMES = f.readline().split()
    f.close()
except:
    sys.exit("Error: Unable to find [ %s ]" % (ifolder + fileproteinnames))


plotBBox = None
plotPBBox = None
plotCTypes = set([])

filecelulas = "dcells.dat"
try:
    f = open(os.path.join(ifolder, filecelulas))
except:
    sys.exit("Error: Unable to find [ %s ]" % (ifolder + filecelulas))

frames = []
startFrame = True
frame = None
for line in f:
    data = line.split()
    if startFrame == True:
        frame = Frame(data) #Start new Frame
        if frame.ncells != 0: #Is frame empty? Otherwise start reading cells.
            startFrame = False        
        continue

    frame.addCell(Cell(data))

    if frame.ncells == len(frame.cells):
        frame.getFrameBoundingBox()
        frame.getProteinBoundingBox()
        frame.getCellTypes()
        if plotBBox == None:
            plotBBox = frame.bbox
            plotPBBox = frame.proteinBBox
        plotBBox = updateBBox(plotBBox, frame.bbox)
        plotPBBox = updatePBBox(plotPBBox, frame.proteinBBox)
        plotCTypes = plotCTypes.union(frame.ctypes)

        frames.append(frame)
        startFrame = True


plt.switch_backend('Agg')

if movieprotein == True:
    print("Plotting protein frames...")
    
    for iframe in range(len(frames)):
        for pname in PNAMES:

            norm = matplotlib.colors.Normalize(vmin=plotPBBox[pname][0],vmax=plotPBBox[pname][1])

            fig, ax = plt.subplots()

            patches = list()
            for ic in range(frames[iframe].ncells):
                cell = frames[iframe].cells[ic]
                polygon = Polygon(
                    cell.vertexs,
                    closed = True,
                    facecolor = cmap(norm(cell.proteins[pname])),
                    linewidth = edgelinewidth,
                    edgecolor = edgecolor)
                patches.append(polygon)


            p = PatchCollection(patches, match_original=True, alpha=0.4)
            ax.add_collection(p)

            if colorbar == True:
                fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)

            plt.xlim(plotBBox["xmin"], plotBBox["xmax"])
            plt.ylim(plotBBox["ymin"], plotBBox["ymax"])
            plt.axis('off')

            fn = ofolder + '%s.%s.%s.%s' % (title, iframe, pname, oformat)
            plt.savefig(fn, bbox_inches='tight', dpi = dpi)
            plt.close()


if moviecelltypes == True:
    print("Plotting cell type frames...")
    
    plotCTypes = list(plotCTypes)
    
    colorStep = 1
    if cmapType.N < len(plotCTypes):
        sys.exit("There is not enough colors in the selected color scheme to represent all cell types.")
    if cmapType.N > 20:
        colorStep = int(cmapType.N/len(plotCTypes))
    
    typeColor = {}
    legendElements = list()
    for idx, elem in enumerate(plotCTypes):
        typeColor[elem] = cmapType(idx*colorStep)
        legendElements.append( Patch(edgecolor=typeColor[elem], label= "Type: %s" % elem, facecolor=typeColor[elem]) )


    for iframe, frame in enumerate(frames[-1:]):
        fig, ax = plt.subplots()

        patches = list()
        for ic in range(frame.ncells):
            cell = frame.cells[ic]
            polygon = Polygon(
                cell.vertexs,
                closed = True,
                facecolor = typeColor[cell.type],
                linewidth = edgelinewidth,
                edgecolor = edgecolor)
            patches.append(polygon)


        p = PatchCollection(patches, match_original=True)
        ax.add_collection(p)

        if colorbar == True:
            ax.legend(handles=legendElements, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=22)

        plt.xlim(plotBBox["xmin"], plotBBox["xmax"])
        plt.ylim(plotBBox["ymin"], plotBBox["ymax"])
        plt.axis('off')

        fn = ofolder + '%s.%s.%s.%s' % (title, iframe, "type", oformat)
        plt.savefig(fn, bbox_inches='tight', dpi = dpi)
        plt.close()

