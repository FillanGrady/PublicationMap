from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np
import datetime
import matplotlib.animation as manimation


def load_papers(filename):
    papers = []
    with open(filename) as f:
        for line in f:
            papers.append(Paper(line))
    return papers

class Paper:
    def __init__(self, line):
        words = line.strip().split(",")
        self.title = words[0]
        if words[2] == 'None':
            self.date = None
        else:
            self.date = datetime.datetime.strptime(words[2], "%Y-%m-%d")
        self.coords = []
        for institution, country, lat, long in zip(*[iter(words[3:])]*4):
            try:
                self.coords.append((float(lat), float(long)))
            except ValueError:
                self.coords.append((0, 0))


def draw_globe(satellite_lat, satellite_long, lats, longs, date, writer):
    if date > datetime.date.today():
        date = datetime.date.today()
    satellite_lat %= 360
    satellite_long %= 360
    map = Basemap(projection='ortho',lat_0=satellite_lat,lon_0=satellite_long,resolution='c')
    # draw coastlines, country boundaries, fill continents.
    map.drawcoastlines(linewidth=0.25, zorder=4, color='w')
    map.drawcountries(linewidth=0.25, zorder=2, color='w')
    # draw the edge of the map projection region (the projection limb)
    map.drawmapboundary(fill_color='black', zorder=0)
    # draw lat/lon grid lines every 30 degrees.
    map.drawmeridians(np.arange(0,360,30), color='w')
    map.drawparallels(np.arange(-90,90,30), color='w')
    lats = np.array([lats])
    lons = np.array([longs])
    # compute native map projection coordinates of lat/lon grid.
    x, y = map(lons, lats)
    # contour data over the map.
    cs = map.scatter(x,y, zorder=4)
    plt.title(str(date))
    writer.grab_frame()
    plt.clf()


if __name__ == "__main__":
    lats = []
    longs = []
    dates = []
    papers = load_papers("Output.csv")
    for paper in papers:
        if paper.date is not None:
            for coord in paper.coords:
                dates.append(paper.date)
                lats.append(coord[0])
                longs.append(coord[1])
    lats = np.array(lats)
    longs = np.array(longs)
    dates = np.array(dates)
    satellite_long = 0
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='Glioblastoma Multiforme', artist='MSTP Group 3',
                    comment='')
    writer = FFMpegWriter(fps=25, metadata=metadata, bitrate=4e3)
    fig = plt.figure()

    with writer.saving(fig, "SpinningGlobe.mp4", 100):
        for year in range(2015, 2018, 1):
            for month in range(1, 12, 1):
                for day in range(1, 28, 1):
                    current_date = datetime.datetime(year=year, month=month, day=day)
                    after_longs = longs[dates < current_date]
                    after_lats = lats[dates < current_date]
                    draw_globe(20, satellite_long, after_lats, after_longs, current_date.date(), writer)
                    satellite_long -= 1.5
                    print "%s-%s-%s" % (year, month, day)
