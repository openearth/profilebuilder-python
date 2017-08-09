from setuptools import setup, find_packages

setup(
    name='profilebuilder',
    version='0.1',
    description='Toolbox for creating flume profiles',
    author = 'Bas Hoonhout',
    author_email='bas.hoonhout@deltares.nl',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'docopt',
        'six',
    ],
    tests_require=[
        'nose'
    ],
    test_suite='nose.collector',
    entry_points={'console_scripts': [
    ]},
)
