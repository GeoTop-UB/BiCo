# BiCo<br><sub><sup>Sage Script for Computing Invariants of Bigraded Complexes</sup></sub>

This Sage script is aimed at computing invariants of Bigraded Complexes. Its main goals are the computation of (unchecked ones are upcoming features):

- [x] Dolbeault, anti-Dolbeault, Bott-Chern and Aeppli cohomologies, 
- [x] Higher pluripotential operations of arity 3,
- [ ] Zigzag decomposition,
- [ ] Frölicher spectral sequence, and
- [ ] LaTeX formatting of the outputs.

It has been funded by the project **Europa Excelencia "Homotopical Invariants of Almost Complex Manifolds" (EUR2023-143450), AEI, Spain**.

## Set up

The only prerequisites for using this Sage script are to have both [Python](https://www.python.org/) and [Sage](https://www.sagemath.org/) installed. Then, one must download the file [bigraded_complexes.py.sage](https://github.com/GeoTop-UB/BiCo/blob/main/bigraded_complexes.py.sage) and paste it in the source folder of the desired project. Afterwards, one has to run/type the following command in the Sage environment/file:

```sage
attach("bigraded_complexes.py.sage")
```

## Example

The following is an example of usage of the Sage script. It loads the predefined Bigraded Complex associated to the Kodaira-Thurston manifold and displays its Aeppli cohomology in ascii art form.

*Input:*

```sage
attach("bigraded_complexes.py.sage")
KT = BidifferentialBigradedCommutativeAlgebraExample.KodairaThurston()
KT._ascii_art_aeppli_cohomology()
```

*Output:*

```txt
Aeppli cohomology:

2  |  [abar*bbar]     |  [b*abar*bbar]                 |  [a*b*abar*bbar]
---+------------------+--------------------------------+-------------------
1  |  [abar], [bbar]  |  [b*abar], [a*bbar], [b*bbar]  |  [a*b*bbar]     
---+------------------+--------------------------------+-------------------
0  |  [1]             |  [a], [b]                      |  [a*b]          
---+------------------+--------------------------------+-------------------
   |  0               |  1                             |  2              
```

## Citation

```bibtex
@misc{BiCoSage,
      title  = {{B}i{C}o: {S}age {S}cript for {C}omputing {I}nvariants of {B}igraded {C}omplexes}, 
      author = {Roger Garrido-Vilallave},
      year   = {2024},
      url    = {https://github.com/GeoTop-UB/BiCo}
}
```
