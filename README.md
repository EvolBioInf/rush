# [`rush`](https://owncloud.gwdg.de/index.php/s/ipd9GfIvaEKNgx4)
## Authors
Bernhard Haubold, Linda Krause, Thomas Horn, Beatriz Vieira Mourato, Peter Pfaffelhuber
## Description
Recombination detection Using SHustrings as further described [here](https://doi.org/10.1093/bioinformatics/btt550).
## Dependencies
`rush` depends on the Gnu scientific library and the divsufsort
library. On Debian-like systems these can be installed with
```
sudo apt install libgsl-dev libdivsufsort-dev
```
On macOS we installed the divsufsort library as explained in its
[repo](https://github.com/y-256/libdivsufsort).
## Compile
Compile the sources using 
```
make
```
the executable `rush` is now in the directory `bin` and can be tested
with
```
make test
```
## Documentation
The command `make doc` generates the manual [`doc/rushDoc.pdf`](https://owncloud.gwdg.de/index.php/s/ipd9GfIvaEKNgx4).
## License
GNU General Public License, https://www.gnu.org/licenses/gpl.html
