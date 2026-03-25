import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

dtype = np.dtype([
    ('frame', np.int32),
    ('x', np.float64),
    ('y', np.float64),
    ('type', np.int32),
])

data = np.fromfile("frames/frames.bin", dtype=dtype)
frame_ids = np.unique(data['frame'])

fig, ax = plt.subplots()
scat_red = ax.scatter([], [], s=2, c="red", label="Carré")
scat_blue = ax.scatter([], [], s=2, c="blue", label="Rectangle")

ax.set_xlim(0, 250)
ax.set_ylim(0, 140)
ax.set_aspect("equal")
ax.legend()

def init():
    scat_red.set_offsets(np.empty((0, 2)))
    scat_blue.set_offsets(np.empty((0, 2)))
    return scat_red, scat_blue

def update(i):
    fid = frame_ids[i]
    frame_data = data[data['frame'] == fid]

    red = frame_data[frame_data['type'] == 1]
    blue = frame_data[frame_data['type'] == 2]

    scat_red.set_offsets(np.c_[red['x'], red['y']] if len(red) else np.empty((0, 2)))
    scat_blue.set_offsets(np.c_[blue['x'], blue['y']] if len(blue) else np.empty((0, 2)))

    ax.set_title(f"Frame {fid}")
    return scat_red, scat_blue

ani = FuncAnimation(fig, update, frames=len(frame_ids), init_func=init, interval=20, blit=True)
plt.show()