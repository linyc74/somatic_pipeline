Stack Overflow page: [x86_64-conda_cos6-linux-gnu-gcc: not found](https://stackoverflow.com/questions/55227065/x86-64-conda-cos6-linux-gnu-gcc-not-found)

Perl uses `cpan` to build packages.
The perl in conda has been problematic with hard-coded paths, which results in build errors like:

```
/bin/sh: 1: /tmp/build/80754af9/perl_1527832170752/_build_env/bin/x86_64-conda_cos6-linux-gnu-gcc: not found
```

To resolve this build issue of `cpan` (and possibly `cpanm`), a very specific version of perl is needed. The version that works is `5.26.2=h470a237_0`.

```bash
conda install -c conda-forge perl=5.26.2=h470a237_0
```

Too bad, yikes!
