from mpl_toolkits.basemap import Basemap
import datetime
import matplotlib.pyplot as plt
from matplotlib import animation


fig = plt.figure()
ax = plt.axes()
line, = ax.plot([], [], lw=0, c='r', alpha=0.7)
m = Basemap(projection='robin', lon_0=0, resolution='c')
lats = []
longs = []
dates = []


def load_papers(filename):
    papers = []
    with open(filename) as f:
        for line in f:
            papers.append(Paper(line))
    return papers


def animate(i):
    x, y = m(longs, lats)
    line.set_data(x, y)
    ax.set_title(str(i+1960))
    ax.plot(x, y, lw=0, c='r', alpha=0.7)
    return line,


def init():
    m.drawcoastlines()
    # draw parallels and meridians.

    m.drawmapboundary(fill_color='black')  # fill to edge
    m.drawcountries()
    m.fillcontinents(color='white', lake_color='white', zorder=0)


class Paper:
    def __init__(self, line):
        words = line.strip().split(",")
        self.title = words[0]
        if words[1] == 'None':
            self.date = None
        else:
            self.date = datetime.datetime.strptime(words[1], "%Y-%m-%d")
        self.coords = []
        for institution, country, lat, long in zip(*[iter(words[2:])]*4):
            self.coords.append((float(lat), float(long)))


if __name__ == "__main__":
    papers = load_papers("Output.csv")
    for paper in papers:
        if paper.date is not None:
            for coord in paper.coords:
                dates.append(paper.date)
                lats.append(coord[0])
                longs.append(coord[1])
    init()
    anim = animation.FuncAnimation(fig, animate, frames=56, interval=60, blit=True)
    anim.save('basic_animation.mp4', fps=30)
    plt.show()