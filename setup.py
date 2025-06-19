from setuptools import setup, find_packages

setup(
    name='ljts_simulator',
    version='0.1',
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "numpy",
        "matplotlib"
    ],
    entry_points={
        "console_scripts": [
            "ljts-run=main:main"
        ]
    }
)
