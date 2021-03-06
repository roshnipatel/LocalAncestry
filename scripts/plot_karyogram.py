"""
Plots karyograms for specified individual.

Modified lightly from Alicia Martin (armartin via GitHub).
"""

import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pylab
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.collections as mcol
import os


def splitstr(option, opt, value, parser):
    return(setattr(parser.values, option.dest, value.split(',')))

parser = argparse.ArgumentParser()

parser.add_argument('--bed_path', help='All bed files for a particular individual (i.e. both haplotypes on all chromosomes)', required=True, nargs='+')
parser.add_argument('--ind', help='Individual sample ID to process.', required=True)
parser.add_argument('--chrX', help='include chrX?', default=False, action="store_true")
parser.add_argument('--centromeres', help='Filepath for information on centromere position.', required=True)
parser.add_argument('--pop_order', default='AFR,EUR,NAT', help='Comma-delimited list of population labels in the order of rfmix populations (1 first, 2 second, and so on). Used in bed files and karyogram labels.')
parser.add_argument('--colors', default=None)
parser.add_argument('--out')

args = parser.parse_args()

def plot_rects(anc, chr, start, stop, hap, pop_order, colors, chrX):
    centro_coords = centromeres[str(chr)]
    if len(centro_coords) == 2: # Acrocentric chromosome
        mask = [
        (centro_coords[0]+2,chr-0.4), # Add +/- 2 at the end of either end
        (centro_coords[1]-2,chr-0.4),
        (centro_coords[1]+2,chr),
        (centro_coords[1]-2,chr+0.4),
        (centro_coords[0]+2,chr+0.4),
        (centro_coords[0]-2,chr),
        (centro_coords[0]+2,chr-0.4)
        ]

        mask_codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        ]
        clip_mask = Path(vertices=mask, codes=mask_codes)

    else: # Need to write more complicated clipping mask with centromere masked out
        mask = [
        (centro_coords[0]+2,chr-0.4), # Add +/- 2 at the end of either end
        (centro_coords[1]-2,chr-0.4),
        (centro_coords[1]+2,chr+0.4),
        (centro_coords[2]-2,chr+0.4),
        (centro_coords[2]+2,chr),
        (centro_coords[2]-2,chr-0.4),
        (centro_coords[1]+2,chr-0.4),
        (centro_coords[1]-2,chr+0.4),
        (centro_coords[0]+2,chr+0.4),
        (centro_coords[0]-2,chr),
        (centro_coords[0]+2,chr-0.4)
        ]

        mask_codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        ]
        clip_mask = Path(vertices=mask, codes=mask_codes)

    if hap == 'A': # bed_a ancestry goes on top
        verts = [
            (float(start), chr), # Left, bottom
            (float(start), chr + 0.4), # Left, top
            (float(stop), chr + 0.4), # Right, top
            (float(stop), chr), # Right, bottom
            (0, 0), # Ignored
        ]
    else: # bed_b ancestry goes on bottom
        verts = [
            (float(start), chr - 0.4), # Left, bottom
            (float(start), chr), # Left, top
            (float(stop), chr), # Right, top
            (float(stop), chr - 0.4), # Right, bottom
            (0, 0), # Ignored
        ]

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]

    clip_path = Path(verts, codes)
    if anc in pop_order:
        col=mcol.PathCollection([clip_path],facecolor=colors[pop_order.index(anc)], linewidths=0)
    else:
        col=mcol.PathCollection([clip_path],facecolor=colors[-1], linewidths=0)
    if 'clip_mask' in locals():
        col.set_clip_path(clip_mask, ax.transData)
    ax.add_collection(col)


# Read in bed files and get individual name
pop_order = args.pop_order.split(',')
ind = args.ind

# Define plotting space
fig = plt.figure()
ax = fig.add_subplot(111)
ax.set_xlim(-5,300)
chrX = args.chrX
if chrX:
    ax.set_ylim(24,0)
else:
    ax.set_ylim(23,0)
plt.xlabel('Genetic position (cM)')
plt.ylabel('Chromosome')
plt.title(ind)
if args.chrX:
    plt.yticks(range(1,24))
    yticks = range(1,23)
    yticks.append('X')
    ax.set_yticklabels(yticks)
else:
    plt.yticks(range(1,23))

# Define colors
def hex_to_rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

colors = []
color_list = args.colors.split(',')
[colors.append(x) for x in color_list]

# Define centromeres
centro = open(args.centromeres)
centromeres = {}
first_line = True
for line in centro:
    if first_line:
        first_line = False
        continue
    line = line.strip().split('\t')
    if chrX and line[0] == 'X':
        line[0] = '23'
    centromeres[line[0]] = [float(n) for n in line[5:]]

# Plot rectangles
for filepath in args.bed_path:
    bed = open(filepath)

    # Extract chromosome and haplotype information from filepath
    s = filepath.strip().split('.')
    chrom = s[1][3:]
    hapl = s[3]

    first_line = True
    for line in bed:
        line = line.strip().split('\t')
        if first_line:
            plot_rects(line[3], int(line[0]), 0, line[2], hapl, pop_order, colors, chrX)
            first_line = False
        try:
            plot_rects(line[3], int(line[0]), line[1], line[2], hapl, pop_order, colors, chrX)
        except ValueError: # Flexibility for chrX
            plot_rects(line[3], 23, line[1], line[2], hapl, pop_order, colors, chrX)
        except IndexError:
            plot_rects(line[3], int(line[0]), line[1], centromeres[str(chrom)][-1], hapl, pop_order, colors, chrX)

# Write a legend
p = []
for i in range(len(pop_order)):
    p.append(plt.Rectangle((0, 0), 1, 1, color=colors[i]))
p.append(plt.Rectangle((0, 0), 1, 1, color='k'))
labs = list(pop_order)
labs.append('UNK')
leg = ax.legend(p, labs, loc=4, fancybox=True)
leg.get_frame().set_alpha(0)

# Get rid of annoying plot features
spines_to_remove = ['top', 'right']
for spine in spines_to_remove:
    ax.spines[spine].set_visible(False)
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

fig.savefig(args.out)
