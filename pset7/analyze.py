from __future__ import print_function

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import openmc


class Viewer(object):
    def __init__(self, sources):
        # Initialize state.
        self.sources = sources
        self.current_frame = 0
        self.points = None
        self.frame_marker = None
        self.t_trace = None

        # Build the figure and axis.
        self.fig = plt.figure()
        #self.ax1 = self.fig.add_subplot(121, aspect='equal')
        #self.ax2 = self.fig.add_subplot(122)
        self.ax1 = self.fig.add_axes((0.07, 0.1, 0.38, 0.8), aspect='equal')
        self.ax2 = self.fig.add_axes((0.55, 0.55, 0.38, 0.38))
        self.ax3 = self.fig.add_axes((0.55, 0.07, 0.38, 0.38))
        self.ax1.set_xlim((-150, 150))
        self.ax1.set_ylim((-150, 150))
        self.ax1.autoscale(enable=False)
        self.ax3.set_xscale('log')
        self.ax3.set_yscale('log')
        self.fig.canvas.mpl_connect('key_press_event', self.on_press)
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)

        # Initialize the plot.
        self.ax2.plot([len(s) for s in self.sources])
        self.t_min = 1e-10
        self.t_max = max(np.max(s['t']) for s in self.sources)
        self.t_bins = np.logspace(np.log10(self.t_min), np.log10(self.t_max),
                                  1000)
        self.update_plot()

    def on_press(self, event):
        # Set the frame number.
        changed = False
        if event.key == 'down' and self.current_frame > 0:
            self.current_frame -= 1
            changed = True
        elif event.key == 'up' and self.current_frame < len(self.sources)-1:
            self.current_frame += 1
            changed = True

        # Update the plot.
        if changed:
            self.update_plot()

    def on_click(self, event):
        if event.inaxes == self.ax2:
            # Figure out which frame the user picked.
            x = int(round(event.xdata))

            # Make sure it's bounded within our data.
            x = min(x, len(self.sources) - 1)
            x = max(x, 0)

            # Update the plot.
            self.current_frame = x
            self.update_plot()

    def update_plot(self):
        # Remove any actors that are about to be replaced.
        if self.points is not None: self.points.remove()
        if self.frame_marker is not None: self.frame_marker.remove()
        if self.t_trace is not None: self.t_trace.remove()

        # Plot the source sites.
        x = self.sources[self.current_frame]['xyz'][:, 0]
        y = self.sources[self.current_frame]['xyz'][:, 1]
        self.points, = self.ax1.plot(x, y, linestyle='', marker='.', c='k')

        # Plot the time distribution.
        #self.ax3.hist(self.sources[self.current_frame]['t'], bins=self.t_bins)
        hist, edges = np.histogram(self.sources[self.current_frame]['t'],
                                   bins=self.t_bins)
        self.t_trace, = self.ax3.step(self.t_bins, stepify(hist), where='pre',
                                      c='k')

        # Indicate the current frame.
        self.frame_marker = self.ax2.axvline(self.current_frame, c='r')

        self.fig.canvas.draw()
        print(self.current_frame)


def stepify(arr):
    """Repeat the first value in an array for step-plotting purposes."""
    out = np.zeros(len(arr) + 1)
    out[0] = arr[0]
    out[1:] = arr
    return out


if __name__ == '__main__':
    # This will plot the absorption rate over the long time scale from 0 to
    # 1000 s.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sp = openmc.StatePoint('statepoint.1000.h5')
    t_bins = sp.tallies[2].filters[0].bins
    rate = stepify(sp.tallies[2].mean[:, 0, 0])
    ax.step(t_bins, rate, where='pre')
    ax.set_xlabel('Time [s]')
    ax.set_xscale('log')
    ax.set_yscale('log')

    # This will plot the absorption rate over the short time scale to 1500
    # mean generation times.
    fig = plt.figure()
    ax = fig.add_subplot(111)
    sp = openmc.StatePoint('statepoint.1000.h5')
    t_bins = sp.tallies[3].filters[0].bins / 1.31e-5
    rate = stepify(sp.tallies[3].mean[:, 0, 0])
    ax.step(t_bins, rate, where='pre')
    ax.set_xlabel('Time / Lambda')

    # This will compute the prompt neutron lifetime.  (It also computes an
    # uncertainty which is clearly a massive over-estimate.)
    sp = openmc.StatePoint('statepoint.100.h5')
    absr_t = sp.tallies[1].mean[0, 0, 0]
    absr_t_sd = sp.tallies[1].std_dev[0, 0, 0]
    absr = sp.tallies[1].mean[0, 0, 1]
    absr_sd = sp.tallies[1].std_dev[0, 0, 1]
    pnl = absr_t / absr
    pnl_sd = pnl * np.sqrt(absr_t_sd**2 / absr_t**2 + absr_sd**2 / absr_sd**2)
    print('Prompt neutron lifetime = {:.2e}'.format(pnl))
    plt.show()
    exit()

    # If you output batch-by-batch statepoints this will open a GUI for viewing
    # the data.  Click on the upper right to pick a generation to view.  Use the
    # up and down arrow keys to increment or decrement the plotted batch by one.
    # (Make sure to comment out the exit() above.)
    sources = []
    for i in range(1000):
        fname = '/tmp/sp_scratch/statepoint.{:04d}.h5'.format(i+1)
        sp = openmc.StatePoint(fname)
        sources.append(sp.source)

    viewer = Viewer(sources)

    plt.show()
