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

## Examples

### Example 1

In this example, the predefined Bigraded Complex associated to the Kodaira-Thurston manifold is loaded, and its Aeppli cohomology is displayed in ascii art form.

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

### Example 2

In this example, an instance of the ```BidifferentialBigradedCommutativeAlgebra``` class representing the Iwasawa manifold is constructed. The product $\mu_3(ab\overline{a}, \overline{b}, \overline{a}) = c \overline{ac}$ is computed, where $\mu_3$ denotes the pluripotential arity-3 product (for some choice of homotopy transfer data).

*Input:*

```sage
attach("bigraded_complexes.py.sage")
QQi = QuadraticField(-1, 'I')
lie_algebra = LieAlgebra(QQi, 'p,ip,q,iq,z,iz', {
            ('p','q'): {'z':1},
            ('p', 'iq'): {'iz':1},
            ('ip','q'): {'iz':1},
            ('ip', 'iq'): {'z':-1}
        })
acs = Matrix(QQi,6,[
            [0,1,0,0,0,0],
            [-1,0,0,0,0,0],
            [0,0,0,1,0,0],
            [0,0,-1,0,0,0],
            [0,0,0,0,0,1],
            [0,0,0,0,-1,0]
        ])
names = ['a','b','c','abar','bbar','cbar']
Iwasawa = BidifferentialBigradedCommutativeAlgebra.from_nilmanifold(lie_algebra, acs, names, normalization_coefficients=[1/2,1,1,1/2,1,1])
Iwasawa.algebra().inject_variables()
values = Iwasawa.operation(arity=3, operation_bidegree=(-1,-1))
print(values[(a*b*abar, bbar, abar)])
```

It should be noted that the normalization coefficients are purely cosmetic. These are introduced to match the generators used in the literature for the bigraded algebra of the Iwasawa manifold.

*Output:*

```txt
Defining a, b, c, abar, bbar, cbar
c*abar*cbar          
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
