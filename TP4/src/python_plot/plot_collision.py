import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

data = np.loadtxt("frames/frames.txt")
frame_ids = np.unique(data[:, 0]).astype(int)

fig, ax = plt.subplots(figsize=(10, 4))

scat_carre = ax.scatter([], [], s=4, c="red", label="Carré")
scat_rect = ax.scatter([], [], s=4, c="blue", label="Rectangle")

ax.set_xlim(0, 250)
ax.set_ylim(0, 120)
ax.set_aspect("equal")
ax.legend()

def init():
    scat_carre.set_offsets(np.empty((0, 2)))
    scat_rect.set_offsets(np.empty((0, 2)))
    return scat_carre, scat_rect

def update(frame_index):
    fid = frame_ids[frame_index]
    frame_data = data[data[:, 0] == fid]

    carre = frame_data[frame_data[:, 4] == 1]
    rect = frame_data[frame_data[:, 4] == 2]

    scat_carre.set_offsets(carre[:, 1:3] if len(carre) else np.empty((0, 2)))
    scat_rect.set_offsets(rect[:, 1:3] if len(rect) else np.empty((0, 2)))

    ax.set_title(f"Frame {fid}")
    return scat_carre, scat_rect

ani = FuncAnimation(fig, update, frames=len(frame_ids), init_func=init, interval=40, blit=True)
plt.show()