import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mt-dynamics", # Replace with your own username
    version="0.0.1",
    author="Florian Huber and Maurits Kok",
    author_email="f.huber@esciencecenter.nl",
    description="MT dynamics simulation",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/florian-huber/mt-dynamics",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)