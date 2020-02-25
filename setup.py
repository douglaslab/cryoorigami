import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cryoorigami-turalaksel", # Replace with your own username
    version="0.0.1",
    author="Tural Aksel",
    author_email="turalaksel@gmail.com",
    description="Set of python scripts for cryo-EM data analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/douglaslab/pyOrigamiEM",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = ['mrcfile==1.1.0', 
                        'numpy==1.14.0',
                        'scipy==1.2.0',
                        'pandas==0.22.0']
    python_requires='>=3.5',
)