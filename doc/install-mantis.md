# MANTIS Installation

Install MANTIS from source.

```bash
cd ~/opt
wget https://github.com/OSU-SRLab/MANTIS/archive/refs/tags/v1.0.5.tar.gz
tar xzf v1.0.5.tar.gz
rm v1.0.5.tar.gz

export PATH=$PATH:$HOME/opt/MANTIS-1.0.5  # in .bashrc
```

Fix Python bug manually due to incompatibility with Python 3.12. Comment out the following line in the `mantis.py` file:

```python
# required_modules_present(['Pysam', 'NumPy'])
```

Fix the `default.py` file:

```python
with open(filepath, 'Ur') as f:  # line 44

# change to
with open(filepath, 'r') as f:
```

Compile the `RepeatFinder` in MANTIS.

```bash
cd ~/opt/MANTIS-1.0.5/tools
make
mv RepeatFinder ../
```

*In fact, MANTIS is so old that there are many other deprecated syntax needed to be fixed. I have my own fork of it.*
