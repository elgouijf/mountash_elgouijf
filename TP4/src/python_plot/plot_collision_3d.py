import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Lecture du fichier
df = pd.read_csv(
    "frames_3d/frames.txt",
    sep=r"\s+",
    header=None,
    names=["frame", "x", "y", "z", "type"]
)

frames = sorted(df["frame"].unique())

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")

def update(frame_id):
    ax.clear()

    data = df[df["frame"] == frame_id]

    data1 = data[data["type"] == 1]
    data2 = data[data["type"] == 2]

    ax.scatter(data1["x"], data1["y"], data1["z"], s=8, label="Type 1")
    ax.scatter(data2["x"], data2["y"], data2["z"], s=8, label="Type 2")

    ax.set_title(f"Frame {frame_id}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    ax.set_xlim(df["x"].min(), df["x"].max())
    ax.set_ylim(df["y"].min(), df["y"].max())
    ax.set_zlim(df["z"].min(), df["z"].max())

    ax.legend()

ani = FuncAnimation(fig, update, frames=frames, interval=80, repeat=True)

plt.show()