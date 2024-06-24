from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

setup(
    name='analysispkg',
    version='0.0.1.dev1',
    description='Python-based analysis package for OPP FA12.',
    long_description=readme,
    author='Stephanne Taylor, Hauke Blanken, Maxim Krassovski, Michael Dunphy',
    packages=find_packages(exclude='bin,config'),
    entry_points={
        'console_scripts': ['scan.py=analysispkg.scr_scan:main',
                            'extract.py=analysispkg.scr_extract:main',
                            'analyze.py=analysispkg.scr_analyze:main',
                            'scores.py=analysispkg.scr_scores:main',
                            'plots.py=analysispkg.scr_plots:main',
                            'score_comparisons.py=analysispkg.scr_score_comparisons:main',
                            'plot_comparisons.py=analysispkg.scr_plot_comparisons:main',
                            ],
    },
)
