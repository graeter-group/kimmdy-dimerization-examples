import matplotlib.pyplot as plt
import seaborn as sns
from PIL import Image

# Nature guide
single_column_mm = 88
double_column_mm = 180


def mm_to_inch(x_mm):
    return x_mm * 0.0393700787


single_column, double_column = mm_to_inch(single_column_mm), mm_to_inch(double_column_mm)

# --- Set here the color by category ---
# The colors are defined below
category_color = {
    "experiment": "HITS_CYAN",
    "experiment_light": "HITS_CYAN_LIGHT",
    "kimmdy": "HITS_MAGENTA",
    "kimmdy_light": "HITS_MAGENTA_LIGHT",
}
# --- ------------------------------ ---

sns.set_context("paper")  # options: paper, talk, poster
sns.set_theme(style="ticks")

plot_colors = {
    # HITS Blues
    "HITS_DARKBLUE": "#003063",
    "HITS_SIGNALBLUE": "#004f9f",

    "HITS_CYAN": "#0088c2",
    "HITS_CYAN_LIGHT": "#55b4dc",

    "HITS_GREEN": "#019050",
    "HITS_GREEN_LIGHT": "#89b77a",

    "HITS_MAGENTA": "#c3006b",
    "HITS_MAGENTA_LIGHT": "#da7da6",

    "HITS_YELLOW": "#ffcc00",
    "HITS_YELLOW_LIGHT": "#ffe07d",

    "MPI_GREEN": "#006c66",
    "MPI_GREEN_SECONDARY": "#055",

    "MPG_grey_dark": "#777777"  # https://www.mpikg.mpg.de/6339023/Logo-Guide-Print-_-Max-Planck-Gesellschaft.pdf
}
# Better would be |= but do this to be compatible with >=3.5
plot_colors = {**plot_colors, **{category: plot_colors[color_ref] for category, color_ref in category_color.items()}}


default_plot_config = {
    # Font sizes (good for figsize ~ (10, 10); adjust as needed for smaller/larger figures)
    'xtick.labelsize': 7,            # Font size for x-axis tick labels
    'ytick.labelsize': 7,            # Font size for y-axis tick labels
    'axes.titlesize': 7,             # Font size for the axes (subplot) title
    'axes.labelsize': 7,             # Font size for axis labels (x and y labels)
    'legend.fontsize': 7,            # Font size for legend text
    'figure.titlesize': 7,           # Font size for the figure-level suptitle

    # Line and border thickness
    'lines.linewidth': 1.5,          # Width of plotted lines
    'axes.linewidth': 0.8,           # Width of axes border/spines
    'xtick.major.width': 0.8,        # Width of major x-axis ticks
    'ytick.major.width': 0.8,        # Width of major y-axis ticks
    'xtick.major.size': 3,
    'ytick.major.size': 3,

    # Tick visibility
    'xtick.top': True,                # Show ticks on top of x-axis
    'ytick.right': True,              # Show ticks on right of y-axis
    'xtick.bottom': True,             # Show ticks on bottom of x-axis
    'ytick.left': True,               # Show ticks on left of y-axis

    # Axes spines (borders)
    'axes.spines.top': True,          # Show top spine
    'axes.spines.right': True,        # Show right spine
    'axes.spines.bottom': True,       # Show bottom spine
    'axes.spines.left': True,         # Show left spine

    # Tick direction (in, out, inout)
    'xtick.direction': 'inout',       # Ticks on x-axis point both in and out
    'ytick.direction': 'inout',       # Ticks on y-axis point both in and out

    # Marker and font settings
    'lines.markersize': 2,           # Size of markers on lines
    'font.sans-serif': ['Roboto'],    # Set sans-serif font; list format ensures fallback options
    'axes.labelcolor': 'black',       # Color for axis labels

    # Colors for axes and ticks
    'axes.edgecolor': 'black',        # Color of the axes border/spines
    'xtick.color': 'black',           # Color of x-axis tick labels and ticks
    'ytick.color': 'black',           # Color of y-axis tick labels and ticks

    # Background transparency (RGBA, where A=0 means fully transparent)
    # 'savefig.facecolor': (0.0, 0.0, 0.0, 0.0),  # Transparent background when saving the figure
    # 'figure.facecolor': (0.0, 0.0, 0.0, 0.0),   # Transparent background of the figure itself
    # 'axes.facecolor': (0.0, 0.0, 0.0, 0.0),     # Transparent background inside each axes

    # Legend appearance
    'legend.frameon': False,          # Don't draw a frame (box) around the legend
    'legend.framealpha': 0,           # Fully transparent legend background (has no effect if frame is off)

    # Saving quality
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',

    # Grid style
    'grid.color': '#EAEAEA',
    'grid.linestyle':  '-',
    'grid.linewidth': 0.7,

    # Nature guide
    'pdf.fonttype': 42,

    "figure.figsize": [single_column, single_column/1.33333]

}


def apply_plot_config(plot_config):
    """Update matplotlib plotting settings"""
    plt.rcParams.update(plot_config)


def init_roboto_font(path_to_roboto_folder=None):
    """Add Roboto font option to matplotlib, if not yet available"""
    import matplotlib.font_manager as fm
    from pathlib import Path

    cwd = Path(__file__).parent
    if not path_to_roboto_folder:
        path_to_roboto_folder = cwd / "assets" / "roboto"

    font_dirs = [path_to_roboto_folder]
    font_files = fm.findSystemFonts(fontpaths=font_dirs)

    for font_file in font_files:
        fm.fontManager.addfont(font_file)


def rgb2hex(r, g, b):
    """Helper function RGB to HEX color conversion"""
    return f"#{r:02x}{g:02x}{b:02x}"


def auto_init():
    apply_plot_config(plot_config=default_plot_config)  # Apply custom matplotlib style
    init_roboto_font()
    return plot_colors


def convert_to_rgb(outputfile):
    img = Image.open(outputfile)
    print(f"File {outputfile} mode BEFORE: {img.mode}")

    print("Convert to RGB")
    img = img.convert("RGB")
    img.save(outputfile)

    img = Image.open(outputfile)
    print(f"File {outputfile} mode AFTER: {img.mode}")

    mode_after = img.mode
    return mode_after


if __name__ == "__main__":
    # Show all colors. Modified script from https://matplotlib.org/stable/gallery/color/named_colors.html
    from matplotlib.patches import Rectangle
    import math

    def plot_color_panel(file_name, colors):

        ncols = 3

        cell_width = 312
        cell_height = 22
        swatch_width = 48
        margin = 12

        names = colors.keys()

        n = len(names)
        nrows = math.ceil(n / ncols)

        width = cell_width * ncols + 2 * margin
        height = cell_height * nrows + 2 * margin
        dpi = 72

        fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
        fig.subplots_adjust(margin / width, margin / height,
                            (width - margin) / width, (height - margin) / height)
        ax.set_xlim(0, cell_width * ncols)
        ax.set_ylim(cell_height * (nrows - 0.5), -cell_height / 2.)
        ax.yaxis.set_visible(False)
        ax.xaxis.set_visible(False)
        ax.set_axis_off()

        for i, name in enumerate(reversed(names)):
            row = i % nrows
            col = i // nrows
            y = row * cell_height

            swatch_start_x = cell_width * col
            text_pos_x = cell_width * col + swatch_width + 7

            ax.text(text_pos_x, y, name, fontsize=14,
                    horizontalalignment='left',
                    verticalalignment='center')
            ax.add_patch(
                Rectangle(xy=(swatch_start_x, y - 9), width=swatch_width,
                          height=18, facecolor=colors[name], edgecolor='0.7')
            )

        plt.savefig(file_name, bbox_inches='tight')

    plot_color_panel("color_categories.png", {category: plot_colors[color_ref] for category, color_ref in category_color.items()})
    plot_color_panel("color_palette.png", {category: color for category, color in plot_colors.items()
                                           if category not in category_color.keys()})