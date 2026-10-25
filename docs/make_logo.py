"""Draws the occulta logo (docs/logo.svg, docs/logo.png, docs/logo_wordmark.png): a star, the companion that
occults it, and the flat-bottomed dip that the occultation leaves in the light curve."""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch
from matplotlib.colors import LinearSegmentedColormap

NAVY, LINE = '#0B1B3A', '#EAF0FA'


def draw(ax, wordmark=False):
    ax.set_xlim(0, 512); ax.set_ylim(512, 0); ax.set_aspect('equal'); ax.axis('off')
    ax.add_patch(FancyBboxPatch((0, 0), 512, 512, boxstyle='round,pad=0,rounding_size=96', fc=NAVY, ec='none'))
    cx, cy, R = 256, 214, 118
    # glow
    for r, a in zip(np.linspace(170, R, 14), np.linspace(0.0, 0.18, 14)):
        ax.add_patch(Circle((cx, cy), r, fc='#F7C25A', ec='none', alpha=a * 0.35))
    # limb-darkened star from a radial gradient
    y, x = np.mgrid[-1:1:400j, -1:1:400j]
    r = np.hypot(x + 0.16, y + 0.2)
    cmap = LinearSegmentedColormap.from_list('star', ['#FFF4D6', '#F7C25A', '#E8743B'])
    img = ax.imshow(np.clip(r, 0, 1.25), cmap=cmap, extent=(cx - R, cx + R, cy + R, cy - R), zorder=2)
    img.set_clip_path(Circle((cx, cy), R, transform=ax.transData))
    # the occulting companion
    ax.add_patch(Circle((330, 222), 54, fc=NAVY, ec='#FFE7AE', lw=2.5, zorder=3))
    # the light curve and its dip
    ax.plot([72, 196, 222, 290, 316, 440], [430, 430, 470, 470, 430, 430], color=LINE, lw=7, solid_joinstyle='round',
            solid_capstyle='round', zorder=4)


fig = plt.figure(figsize=(5.12, 5.12), dpi=100)
ax = fig.add_axes([0, 0, 1, 1]); draw(ax)
fig.savefig('docs/logo.png', dpi=100, transparent=True); fig.savefig('docs/logo.svg', transparent=True); plt.close(fig)

fig = plt.figure(figsize=(12.8, 4.0), dpi=100)
ax = fig.add_axes([0.0, 0.0, 0.3125, 1.0]); draw(ax)
ax2 = fig.add_axes([0.33, 0.0, 0.67, 1.0]); ax2.axis('off')
ax2.text(0.0, 0.5, 'occulta', fontsize=96, family='serif', color=NAVY, va='center', fontweight='bold')
fig.savefig('docs/logo_wordmark.png', dpi=100, transparent=True); plt.close(fig)
print('ok')
