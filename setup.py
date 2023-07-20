import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cryoorigami", # Replace with your own username
    version="1.0.0",
    author="Tural Aksel",
    author_email="turalaksel@gmail.com",
    description="Python toolset for cryo-EM data analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/douglaslab/cryoorigami",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=['matplotlib==3.1.2',
                      'mrcfile==1.1.0',
                      'numba==0.42.0',
                      'numpy==1.17.4',
                      'packaging==19.2',
                      'pandas==0.25.3',
                      'parso==0.5.1',
                      'pexpect==4.7.0',
                      'pickleshare==0.7.5',
                      'prompt-toolkit==3.0.2',
                      'ptyprocess==0.6.0',
                      'pyFFTW==0.11.1',
                      'Pygments==2.15.0',
                      'pyparsing==2.4.5',
                      'python-dateutil==2.8.1',
                      'pytz==2019.3',
                      'PyYAML==5.4',
                      'scipy==1.2.0',
                      'sip==5.0.0',
                      'six==1.13.0',
                      'termcolor==1.1.0',
                      'toml==0.10.0',
                      'traitlets==4.3.3',
                      'wcwidth==0.1.7'],
    python_requires='>=3.5',
)
