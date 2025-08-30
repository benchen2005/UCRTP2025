import subprocess

cmd = [
    "python",
    "/full/path/to/src/rmats2sashimiplot/rmats2sashimiplot.py",
    "--b1", "/path/to/control1.bam,/path/to/control2.bam",
    "--b2", "/path/to/kd1.bam,/path/to/kd2.bam",
    "-e", "/path/to/events_to_plot.txt",
    "--l1", "Control",
    "--l2", "KD",
    "--exon_s", "1",
    "--intron_s", "5",
    "-o", "/path/to/output_dir"
]

subprocess.run(cmd)
