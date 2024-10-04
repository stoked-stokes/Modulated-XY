import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
from matplotlib import animation
import imageio

L = 30
EQ_NO = 100 * L
MEAS_NO = 1000

hsv = cm.get_cmap('hsv', 12)

def angle_to_color(angle):
    if angle < 0:
        angle += 2 * np.pi
    
    return hsv((angle / (2*np.pi)))

# def animate(i):

#     U = np.cos(data[EQ_NO+1+i])
#     V = np.sin(data[EQ_NO+1+i])
#     print(np.array([angle_to_color(i) for i in data[EQ_NO+1+i].flatten()]))
    
#     plt.title(f'frame {1+i}')
#     qu.set_UVC(U, V, C=np.array([angle_to_color(i) for i in data[EQ_NO+1+i].flatten()]))

run = 21
temp = 0.1
filename = f"./Data/run{run}/config-{temp}.h5"

with h5py.File(filename, "r") as file:

    data = dict()
    for i in range(0, EQ_NO+MEAS_NO):
        matrix = file["matrix"][i, :, :]
        data[i] = np.transpose(matrix)

    file.close()

x = np.arange(1, L+1, 1)
y = np.arange(1, L+1, 1)
X, Y = np.meshgrid(x, y)

def create_frame(i):
    U = np.cos(data[i])
    V = np.sin(data[i])
    colorscheme = np.array([angle_to_color(j) for j in data[i].flatten()])

    fig = plt.figure()
    plt.quiver(X, Y, U, V, units='xy', scale=1.1, pivot='mid', color=colorscheme)
    plt.xlim(-1, L+2)
    plt.ylim(-1, L+2)
    plt.title(f'frame {i}')

    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig(f'./Animation/run{run}/{temp}/img_{i}.png', 
                transparent = False,  
                facecolor = 'white'
               )
    # plt.show()
    plt.close()

a = EQ_NO+1
b = EQ_NO+MEAS_NO

for k in range(a, b):
    create_frame(k)

frames = []
for k in range(a, b):
    image = imageio.imread(f'./Animation/run{run}/{temp}/img_{k}.png') # throws the error "module 'imageio' has no attribute 'v2'"
    frames.append(image)

imageio.mimsave(f'XY-{run}@{temp}.gif', # output gif
                frames,          # array of input frames
                fps = 0.7)         # optional: frames per second



# array = np.array([])
# for i in range((100*L)+1, ITERS+1):
#     array = np.append(array, data[i])

# c = np.array(list(map(ang_to_color, array.flatten())))

# anim = animation.FuncAnimation(fig, animate, frames=MEAS_NO, interval=400, repeat=False)